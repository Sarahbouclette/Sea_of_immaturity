rm(list=ls())

##-----------------Loading packages-------------------
pkgs <- c("here", "dplyr", "ggplot2", "funbiogeo", "missForest", "slam",
          "pbmcapply", "patchwork", "ggplot2", "maditr","rnaturalearth",
          "tibble", "stringr", "hrbrthemes", "randomForest", "ranger", "caret", "ggpubr", "pheatmap")
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

library(funbiogeo)
library(ggplot2)

##-------------loading data and functions-------------
# species traits 
load(file = here::here("Documents", "Sea_of_immaturity","outputs", "MF2_species_traits_imputed.Rdata"))
load(file = here::here("Documents", "Sea_of_immaturity","data", "derived_data", "01_d_species_traits_final.Rdata"))

##-----------1.Preparing data and selection of explicative variables------------

#Add the phylogeny predictive variable doing an MCA

var_phylo <- c("species_name", "Class", "Order", "Family", "Genus") 

phylogeny <- species_traits_final |> 
  dplyr::select(all_of(var_phylo)) |> 
  tibble::column_to_rownames(var="species_name") |> 
  as.matrix()

dimensionality <- FactoMineR::MCA(phylogeny, ncp = 4, graph=F) 

phylo_space <- dimensionality[["ind"]][["coord"]]
colnames(phylo_space)<- gsub(" ", "", colnames(phylo_space))

traits_and_phylo <- cbind(phylo_space, species_traits_final)
Dims <- colnames(traits_and_phylo[,1:4])

# Add the already known length at maturity to the data inferred with Miss forest

traits_and_phylo <- traits_and_phylo %>%
  select(species_name, LengthMatMin_unsexed, all_of(Dims))

#species_traits_imputed <- species_traits_imputed %>%
  #rename(species_name = species)

data_rf <- species_traits_imputed %>%
  left_join(traits_and_phylo, by = "species_name")

# Liste des prédicteurs potentiels basés sur des connaissances écologiques
selected_traits <- c("Length", "Troph", "TempPrefMean", "Climate", 'K', "Vulnerability_fishing", "IUCN_inferred_Loiseau23", "TempPrefM","TempPrefMin", "LengthMax_unsexed", "Status", "DemersPelag")

# Filtrer les données sans valeurs manquantes pour l’entraînement
data_rf_known <- data_rf |>
  tidyr::drop_na(all_of(selected_traits), LengthMatMin_unsexed, Dims) %>%  # Supprime les lignes avec NA uniquement dans ces colonnes
  select(all_of(selected_traits), LengthMatMin_unsexed, Dims, all_of(var_phylo)) #Species with known LengthMat and known predictors
data_rf_unknown <- data_rf[is.na(data_rf$LengthMatMin_unsexed), ] #Species with unknown maturity Length 
data_rf_unknown <- data_rf_unknown %>% 
  tidyr::drop_na(all_of(selected_traits), Dims) %>% #but known predictors
  select(all_of(selected_traits), Dims, all_of(var_phylo)) 

##-----------2. Evaluation des variables predictrices------------------------

# Sélection des variables optimales avec RFE (Recursive Feature Elimination)
cv_control <- rfeControl(functions = rfFuncs, method = "cv", number = 5)

rfe_result <- rfe(data_rf_known %>% select(-LengthMatMin_unsexed, -all_of(var_phylo)), data_rf_known$LengthMatMin_unsexed,
                  sizes = c(3, 5, 6), rfeControl = cv_control)

best_vars <- predictors(rfe_result)  # Variables retenues par RFE
print(best_vars)

#best_vars <- c("Length", "Troph", "TempPrefMean", "Climate", 'K', "Vulnerability_fishing", "IUCN_inferred_Loiseau23", "TempPrefM","TempPrefMin", "LengthMax_unsexed", "Status", "DemersPelag")

##---------------3.optimization Random Forest-------------
tune_grid <- expand.grid(mtry = c(2, 5, 10),
                         splitrule = c("variance",NULL),
                         min.node.size = c(1, 5, 10))

cv_control <- trainControl(method = "cv", number = 5)

bestTune<-data.frame(num_trees_values = c(100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050), MSE = NA, mtry = NA, splitrule=NA, min.node.size = NA)

for (n in bestTune$num_trees_values) {
  rf_tuned <- train(
    LengthMatMin_unsexed ~ ., 
    data = data_rf_known %>% select(all_of(best_vars), LengthMatMin_unsexed),
    method = "ranger", 
    trControl = cv_control, 
    tuneGrid = tune_grid,
    num.trees = n
  )
  
  bestTune[bestTune$num_trees_values==n,]$MSE <- rf_tuned$finalModel$prediction.error
  bestTune[bestTune$num_trees_values==n,]$mtry <- rf_tuned$bestTune$mtry
  bestTune[bestTune$num_trees_values==n,]$splitrule <- rf_tuned$bestTune$splitrule
  bestTune[bestTune$num_trees_values==n,]$min.node.size <- rf_tuned$bestTune$min.node.size
}


print(rf_tuned$bestTune)  # Meilleurs paramètres trouvés

# Erreur en fonction du nombre d'arbres
plot(bestTune$num_trees_values, bestTune$MSE,
     type='l',
     main = "Erreur en fonction du nombre d'arbres", 
     xlab = "Nombre d'arbres", 
     ylab = "Erreur")

##----------------4. Training Random Forest-------------------------

set.seed(123)  # Reproductibilité
num_folds <- 5  # Nombre de folds

folds <- createFolds(data_rf_known$LengthMatMin_unsexed, k = num_folds, list = TRUE)

all_predictions <- list()  # Stocker les prédictions par fold
cv_metrics <- data.frame(Fold = integer(), R2 = numeric(), Type = character())  # Stocke les R²

for (i in seq_along(folds)) {
  train_idx <- unlist(folds[-i])  # Les autres folds servent pour le train
  test_idx <- folds[[i]]  # Le fold i sert pour le test
  
  train_fold <- data_rf_known[train_idx, ]
  test_fold <- data_rf_known[test_idx, ]
  
  rf_model <- ranger(LengthMatMin_unsexed ~ ., 
                     data = train_fold %>% select(all_of(best_vars), LengthMatMin_unsexed),
                     num.trees = 500, importance = "permutation",
                     mtry = 2, splitrule = "variance", min.node.size = 10)
  
  # Prédictions
  pred_train <- predict(rf_model, train_fold %>% select(all_of(best_vars)))$predictions
  pred_test <- predict(rf_model, test_fold %>% select(all_of(best_vars)))$predictions
  
  # Calcul du R² pour Train
  ss_res_train <- sum((train_fold$LengthMatMin_unsexed - pred_train)^2)
  ss_tot_train <- sum((train_fold$LengthMatMin_unsexed - mean(train_fold$LengthMatMin_unsexed))^2)
  r2_train <- 1 - (ss_res_train / ss_tot_train)
  
  # Calcul du R² pour Test
  ss_res_test <- sum((test_fold$LengthMatMin_unsexed - pred_test)^2)
  ss_tot_test <- sum((test_fold$LengthMatMin_unsexed - mean(test_fold$LengthMatMin_unsexed))^2)
  r2_test <- 1 - (ss_res_test / ss_tot_test)
  
  # Stockage des résultats
  cv_metrics <- rbind(cv_metrics, 
                      data.frame(Fold = i, R2 = r2_train, Type = "Train"),
                      data.frame(Fold = i, R2 = r2_test, Type = "Test"))
  
  # Stocker les prédictions
  all_predictions[[i]] <- list(
    train = data.frame(Observed = train_fold$LengthMatMin_unsexed, Predicted = pred_train, Fold = i, Set = "Train"),
    test = data.frame(Observed = test_fold$LengthMatMin_unsexed, Predicted = pred_test, Fold = i, Set = "Test")
  )
}

# Fusionner les prédictions
df_predictions <- do.call(rbind, lapply(all_predictions, function(x) rbind(x$train, x$test)))

# Transformer les noms des folds en facteur pour l'affichage ggplot
cv_metrics$Fold <- as.factor(cv_metrics$Fold)

# 🔹 Création du graphique
ggplot(cv_metrics, aes(x = Fold, y = R2, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "R² par fold (Train vs Test)", x = "Fold", y = "R²") +
  scale_fill_manual(values = c("steelblue", "darkorange"))

ggsave(file = here("Documents", "Sea_of_immaturity", "figures", "RF2_rsquared_by_fold.png"), width = 8, height = 6)

df_test <- df_predictions %>% filter(Set == "Test")
df_train <- df_predictions %>% filter(Set == "Train")
r2_test_tot <- mean(cv_metrics[cv_metrics$Type=="Test",]$R2)
r2_train_tot <- mean(cv_metrics[cv_metrics$Type=="Train",]$R2)

# Visualisation pred train
ggplot(df_train, aes(x = Observed, y = Predicted)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Comparaison Observed vs Predicted on the Train set", 
       x = "Observed values", y = "Predicted values") + 
  annotate("text", x = min(df_train$Observed), y = max(df_train$Predicted), 
           label = paste("R² =", round(r2_train_tot, 3)), 
           hjust = 0, size = 5, color = "blue") 

# Sauvegarde du graphique
ggsave(file = here("Documents", "Sea_of_immaturity", "figures", "RF2_regression_CV_train.png"), width = 8, height = 6)

# Visualisation pred test
ggplot(df_test, aes(x = Observed, y = Predicted)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Comparaison Observed vs Predicted on the Test set", 
       x = "Observed values", y = "Predicted values") + 
  annotate("text", x = min(df_test$Observed), y = max(df_test$Predicted), 
           label = paste("R² =", round(r2_test_tot, 3)), 
           hjust = 0, size = 5, color = "blue") 
# Sauvegarde du graphique
ggsave(file = here("Documents", "Sea_of_immaturity", "figures", "RF2_regression_CV_test.png"), width = 8, height = 6)


##------------5. Prediction of the unknown maturity length----------------------

pred_unknown <- predict(rf_model, data_rf_unknown %>% select(all_of(best_vars)))$predictions
pref_unkown_filtered <- pmin(pred_unknown,  data_rf_unknown$Length)

##-------------6. Integration of the predictions into the initial dataset-------

data_rf_unknown$LengthMatMin_unsexed <- pred_unknown

##-------------7. Plot results preedictions--------------------------

# Extraire l'importance des variables
importance_df <- as.data.frame(rf_model$variable.importance) %>%
  rownames_to_column(var = "Variable") %>%
  rename( importance = `rf_model$variable.importance`)%>%
  arrange(desc(importance))  # Trier par importance décroissante

ggplot(importance_df, aes(x = reorder(Variable, importance), y = importance)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Importance des variables", x = "Variable", y = "Importance") 

ggsave(plot= last_plot(), width = 8, height = 8, 
       file= here::here("Documents", "Sea_of_immaturity", "figures","2_RF_importance_variables.png"))


plot_predictions <- ggplot(data_rf_unknown, aes(x=LengthMatMin_unsexed)) + geom_histogram(color='red') + ggtitle('Distribution of the maturity length predicted')
plot_observed_data <- ggplot(data_rf_known, aes(x=LengthMatMin_unsexed))+ geom_histogram(color='orange') + ggtitle('Distribution of the maturity length already observed')
plot_lengths <- ggplot(data_rf, aes(x=Length)) + geom_histogram(color='skyblue') + ggtitle('Distribution of the lengths')
plot_lengthsMax <-  ggplot(data_rf, aes(x=LengthMax_unsexed)) + geom_histogram(color='blue') + ggtitle('Distribution of the maximum lengths')

# Distribution lengths
ggpubr::ggarrange(plot_predictions, plot_observed_data, plot_lengths, plot_lengthsMax)
summary(data_rf_known$LengthMatMin_unsexed)
summary(data_rf_unknown$LengthMatMin_unsexed)
summary(data_rf$Length)

ggsave(plot= last_plot(), width = 8, height = 8, 
       file= here::here("Documents", "Sea_of_immaturity", "figures","2_RF_predictions_distribution.png"))

##-----------8. Save------------------------------------------------------

##-----------9. Search for causality-----------------------------------

library(ape)
library(pheatmap)
library(dplyr)
library(vegan)

## Length at maturity, phylogeny heritage ? ##

##Similarity based on maturity lengths 
pred_matrix <- as.matrix(pred_unknown)
dist_pred <- as.matrix (dist(pred_matrix))
rownames(dist_pred) <- data_rf_unknown$species_name
colnames(dist_pred) <- data_rf_unknown$species_name
sim_pred <- 1 / (1 + dist_pred)

pheatmap(sim_pred, 
         clustering_method = "ward.D2", 
         main = "Clustering basé sur Random Forest",
         labels_row = rownames(sim_matrix),
         labels_col = colnames(sim_matrix), 
         fontsize_row = 6, fontsize_col = 6)


## Similarity based on the phylogeny
taxo_df <- data_rf_unknown %>%
  select(all_of(var_phylo)) %>%
  na.omit()

taxo_df$species_name <- as.factor(taxo_df$species_name)

tree <- ape::as.phylo(~Class/Order/Family/Genus/species_name, data = taxo_df)

# Ajouter les noms des taxons supérieurs
tree$node.label <- unique(c(taxo_df$Class, taxo_df$Order, taxo_df$Family, taxo_df$Genus))

# Affichage de l'arbre phylogénétique
plot(tree, type = "cladogram", cex = 0.5, no.margin = TRUE, label.offset = 0.5)
nodelabels(tree$node.label, frame = "none", cex = 0.7, col = "red")

# Matrice de distances phylogénétiques
phylo_dist <- cophenetic(tree)  

# Matrice de similarité phylogénétique
phylo_sim <- 1 / (1 + phylo_dist)

# 🔹 S'assurer que les matrices ont les mêmes espèces et le même ordre
common_species <- intersect(rownames(sim_pred), rownames(phylo_sim))

sim_pred_filtered <- sim_pred[common_species, common_species]
phylo_sim_filtered <- phylo_sim[common_species, common_species]

# Vérification des dimensions
print(dim(sim_pred_filtered))
print(dim(phylo_sim_filtered))

## Comparaison des matrices avec un test de Mantel
mantel_test_sim <- vegan::mantel(sim_matrix_filtered, phylo_sim_filtered, method = "spearman")

# Affichage des résultats
print(mantel_test_sim)

# result
# r=0.1793
#significance = 0.001

#ccl: la prédiction de la taille à maturité avec le RF ne semble pas bien reproduire la structure phylogénétique. 


## Correlation between LengthMatMin and Length ? ##








