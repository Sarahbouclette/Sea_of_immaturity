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
# species traits #
load(file = here::here("Documents", "Sea_of_immaturity","outputs", "species_traits_imputed.Rdata"))
load(file = here::here("Documents", "Sea_of_immaturity","data", "derived_data", "01_d_species_traits_final.Rdata"))

# Functions to filter and clean the traits

filter_na_columns <- function(data, na_threshold = 0.1) {
  # Calcul du pourcentage de NA par colonne
  na_percent <- colMeans(is.na(data))
  
  # Sélectionner les colonnes avec un taux de NA inférieur au seuil
  filtered_data <- data[, na_percent < na_threshold]
  
  return(filtered_data)
}

replace_empty_null_with_na <- function(data) {
  data[is.null(data) | data == "" | trimws(data) == ""] <- NA
  return(data)
}

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

species_traits_imputed <- species_traits_imputed %>%
  rename(species_name = species)

data_rf <- species_traits_imputed %>%
  left_join(traits_and_phylo, by = "species_name")

# Liste des prédicteurs potentiels basés sur des connaissances écologiques
selected_traits <- c("Length", "Troph", "TempPrefMean", "Climate", "trophic_guild", 'K', "Vulnerability_fishing", "IUCN_inferred_Loiseau23")

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

##--------3.Division en set d'entraînement (80%) et set de test (20%)--------

set.seed(123) #Reproductibility
train_index <- createDataPartition(data_rf_known$LengthMatMin_unsexed, p = 0.8, list = FALSE)
train_data <- data_rf_known[train_index, ]
test_data <- data_rf_known[-train_index, ]

##---------------4.Training and optimization Random Forest-------------
tune_grid <- expand.grid(mtry = c(2, 5, 10),
                         splitrule = c("variance",NULL),
                         min.node.size = c(1, 5, 10))

custom_metric <- function(data, lev = NULL, model = NULL) {
  predictions <- data$pred
  predictions <- pmin(predictions, data_rf_known$Length)  # Contraindre les valeurs
  
  RMSE <- sqrt(mean((data$obs - predictions)^2))  # Calcul du RMSE
  return(c(RMSE = RMSE))
}

cv_control <- trainControl(method = "cv", number = 5, summaryFunction = custom_metric)

bestTune<-data.frame(num_trees_values = c(100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050), MSE = NA, mtry = NA, splitrule=NA, min.node.size = NA)

for (n in bestTune$num_trees_values) {
  rf_tuned <- train(
    LengthMatMin_unsexed ~ ., 
    data = train_data %>% select(all_of(best_vars), LengthMatMin_unsexed),
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


rf_model <- ranger(LengthMatMin_unsexed ~ ., data = train_data %>% select(all_of(best_vars), LengthMatMin_unsexed),
                   num.trees = 800, importance = "permutation", mtry=10, splitrule = "variance", min.node.size = 1)

##-------------5. Évaluation du modèle sur le set de test--------------------

pred_test <- predict(rf_model, test_data %>% select(all_of(best_vars)))$predictions
rmse_test <- sqrt(mean((pred_test - test_data$LengthMatMin_unsexed)^2))
cat("Test RMSE :", rmse_test, "\n")

##-------------6. Training the model on all the known data----------------------

rf_final <- ranger(LengthMatMin_unsexed ~ ., data = data_rf_known %>% select(all_of(best_vars), LengthMatMin_unsexed),
                   num.trees = 800, importance = "permutation", mtry = 10, splitrule = "variance", min.node.size = 1)

##------------7. Prediction of the unkown maturity length----------------------

pred_unknown <- predict(rf_final, data_rf_unknown %>% select(all_of(best_vars)))$predictions

##-------------8. Integration of the predictions into the initial dataset-------

data_rf_unknown$LengthMatMin_unsexed <- pred_unknown

##-------------9. Plot final performances and results--------------------------

# Extraire l'importance des variables
importance_df <- as.data.frame(rf_model$variable.importance) %>%
  rownames_to_column(var = "Variable") %>%
  rename( importance = `rf_model$variable.importance`)%>%
  arrange(desc(importance))  # Trier par importance décroissante

ggplot(importance_df, aes(x = reorder(Variable, importance), y = importance)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Importance des variables", x = "Variable", y = "Importance") +
  theme_minimal()

# Erreur en fonction du nombre d'arbres
plot(bestTune$num_trees_values, bestTune$MSE,
     type='l',
     main = "Erreur en fonction du nombre d'arbres", 
     xlab = "Nombre d'arbres", 
     ylab = "Erreur")

# Prédictions vs observations
plot(test_data$LengthMatMin_unsexed, pred_test, 
     xlab = "Valeurs réelles", 
     ylab = "Valeurs prédites", 
     main = "Prédictions vs Réel")
abline(0, 1, col = "red", lwd = 2)

# Distribution des erreurs
errors <- test_data$LengthMatMin_unsexed - pred_test
hist(errors, breaks = 20, main = "Distribution des erreurs", xlab = "Erreur (Réel - Prédit)")


plot_predictions <- ggplot(data_rf_unknown, aes(x=LengthMatMin_unsexed)) + geom_histogram(color='red') + ggtitle('Distribution of the maturity length predicted')
plot_observed_data <- ggplot(data_rf_known, aes(x=LengthMatMin_unsexed))+ geom_histogram(color='orange') + ggtitle('Distribution of the maturity length already observed')
plot_lengths <- ggplot(data_rf, aes(x=Length)) + geom_histogram(color='skyblue') + ggtitle('Distribution of the lengths')

# Distribution lengths
ggarrange(plot_predictions, plot_observed_data, plot_lengths)
summary(data_rf_known$LengthMatMin_unsexed)
summary(data_rf_unknown$LengthMatMin_unsexed)
summary(data_rf$Length)

##-----------10. Save------------------------------------------------------

##-----------11. Search for causality-----------------------------------

# Phylogeny clusters ?

## Similarity based on maturity lengths
pred_matrix <- as.matrix(pred_unknown)
proximity_approx <- as.matrix(dist(pred_matrix))
rownames(proximity_approx) <- data_rf_unknown$species_name
colnames(proximity_approx) <- data_rf_unknown$species_name
sim_matrix <- 1 - as.matrix(proximity_approx)

pheatmap(sim_matrix, 
         clustering_method = "ward.D2", 
         main = "Clustering basé sur Random Forest",
         labels_row = rownames(sim_matrix),
         labels_col = colnames(sim_matrix),fontsize_row = 6, fontsize_col = 6)

hc_rf <- hclust(as.dist(sim_matrix), method = "ward.D2")

# Découper en 5 clusters (modifiable selon ton besoin)
clusters <- cutree(hc_rf, k = 5)

# Associer chaque espèce à son cluster
species_clusters <- data.frame(Species = rownames(sim_matrix), Cluster = clusters)
print(species_clusters)

## Similarity based on the phylogeny 

taxo_df <- data_rf_unknown %>%
  select(all_of(var_phylo))

hc <- hclust(dist(1:nrow(taxo_df)), method = "average")  # Clustering basé sur l'ordre des données
phylo_tree <- as.phylo(hc)
plot(phylo_tree, type="cladogram")










