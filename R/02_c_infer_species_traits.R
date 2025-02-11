################################################################################
##
##  
##
## 1c_infer_species_traits.R
##
## 23/10/2022
##
## Ulysse Flandrin
##
###############################################################################"
#----------------- cleaning memory -----------------
rm(list=ls())

##-----------------Loading packages-------------------
pkgs <- c("here", "dplyr", "ggplot2", "funbiogeo", "missForest", "slam",
           "pbmcapply", "patchwork", "ggplot2", "maditr","rnaturalearth",
           "tibble", "stringr", "hrbrthemes")
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

library(funbiogeo)
library(ggplot2)
library(caret)

##-------------loading data and functions-------------
# species traits #
load(file = here::here("Documents", "Sea_of_immaturity","data", "derived_data", "01_d_species_traits_final.Rdata"))

# functions #
source(here::here("Documents", "Sea_of_immaturity","R","02_a_MissForest_functions.R"))
source(here::here("Documents", "Sea_of_immaturity","R","02_b_evaluation_prediction_model.R"))

filter_na_columns <- function(data, na_threshold = 0.3) {
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


##-------------Check for traits distribution-------------
species_traits_final$spec_code <- as.character(species_traits_final$spec_code)
#species_traits_final<-species_traits_final %>%
  #select(-Spawning, -LengthMatMin_female, -LengthMatMin_male, -LengthMatMin_unsexed, -n_T, -LengthFemale, -WeightMax, -FamCode, -GenCode, -IUCN, -Emblematic, -Abyssopelagic, -Phylum, -ParentalCareQ, -range_n_cells_01, -LTypeMaxM) # To simplfy

species_traits_filtered <- filter_na_columns(species_traits_final)
species_traits_filtered <- replace_empty_null_with_na(species_traits_filtered)
species_traits_filtered <- species_traits_filtered %>%
  select(-Phylum, -FamCode, -GenCode, -Emblematic, -LTypeMaxM, -MGillnets, -geographic_range_Albouy19, -UsedforAquaculture, -Salinity, 
         -Area, -SouthernLat, -WesternLat, -NorthernLat, -EasternLat,
         -Aquarium, -Abyssopelagic, -AverageDepth, -ClimVuln_SSP585, -range_n_cells_01)

species_traits_filtered[] <- lapply(species_traits_filtered, function(x) {
  if (is.character(x)) return(factor(x))
  if (is.integer(x)) return(as.numeric(x))
  return(x)
})
# Sélection des colonnes numériques pour la corrélation
numeric_columns <- sapply(species_traits_final, is.numeric)
corr_matrix <- cor(species_traits_final[, numeric_columns], use = 'complete.obs')

# Identifier les colonnes fortement corrélées (cutoff 0.9)
high_corr <- findCorrelation(corr_matrix, cutoff = 0.9)
traits_correlated <- species_traits_final[, numeric_columns]
traits_correlated <- traits_correlated[,high_corr]
traits_correlated <- colnames(traits_correlated)
species_traits_final <- species_traits_final |>
  select(-AverageDepth, -body_depth_ratio, -ClimVuln_SSP126) #Chosen from the traits_correlated list



numeric_cols <- colnames(species_traits_filtered)[sapply(species_traits_filtered, class) == "numeric"]
df <- tidyr::pivot_longer(species_traits_filtered, cols = all_of(numeric_cols) ,
                          names_to = "trait", values_to = "value") |> 
  dplyr::filter(!is.na(value))

ggplot(data = df)+
  aes(x=value, group=trait, fill=trait) +
  geom_histogram(aes(y = ..density.., fill = trait), bins = 20,
                 color = "grey40", alpha = 0.2) +
  facet_wrap(~trait, scales = "free") +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    axis.ticks.x=element_blank()
  )

ggsave(plot=last_plot(), here::here("Documents", "Sea_of_immaturity", "outputs", "Distribution_variables_beforeMF3.png"))

#DEPTH MIN AND DEPTH MAX + GEOGRAPHIC RANGE ARE HIGLY RIGHT-SKEWED -> LOG TRANSFORMATION
cols <- c("depth_min", "depth_max", "a", "b", "Length", "K", "LengthMax_unsexed")

if (all(cols %in% colnames(species_traits_filtered))) {
  species_traits_filtered <- species_traits_filtered %>%
    mutate(across(all_of(cols), ~ log10(. + 1), .names = "log_{.col}")) %>%
    select(-all_of(cols))  # Supprime les anciennes colonnes
}



################################################################################"
##
##       #### Integrate taxonomy with Principal Component Analyses ####
##
###############################################################################"

##-------------Choose the number of dimensions-------------

phylogeny <- species_traits_filtered |> 
  dplyr::select(species_name, Class, Order, Family, Genus) |> 
  tibble::column_to_rownames(var="species_name") |> 
  as.matrix()

#Integrate phylogeny with MCA
dimensionality <- FactoMineR::MCA(phylogeny,ncp=10, graph=F) 
factoextra::fviz_eig(dimensionality, choice = "variance",
                     addlabels = F, ylim = c(0, 0.75), barcolor = "darkred",
                     barfill = "darkred", ncp=550)


# Choose the number of dimensions
dimensions <- c(0,2,5,10,15,20,25,30,35,40,45,50,55,60,75,100,150)
## Run missforest for each dimension of phylogeny
test_dimensionality <- lapply(dimensions, FUN = function(dim){
  #dim=0
  cat(dim, "dimensions:", "\n")
  
  if( dim == 0 ){
    Dims <- paste0("Dim", 1:dim)
    
    #No phylogeny
    traits_and_phylo <- species_traits_filtered
    #prep data
    data_for_missforest <- preping_data(species_traits_df =traits_and_phylo)
    
    data_to_infer <- data_for_missforest[[1]]
    traits_data_factors <- data_for_missforest[[2]]
    traits_data_num <- data_for_missforest[[3]]
    factor_length <- data_for_missforest[[4]]
    
  }else{
    Dims <- paste0("Dim", 1:dim)
    
    dimensionality <- FactoMineR::MCA(phylogeny,ncp=dim, graph=F) 
    
    phylo_space <- dimensionality[["ind"]][["coord"]]
    colnames(phylo_space)<- gsub(" ", "", colnames(phylo_space))
    
    traits_and_phylo <- cbind(phylo_space, species_traits_filtered) 
    rownames(traits_and_phylo) <- NULL
    
    data_for_missforest <- preping_data(species_traits_df = traits_and_phylo)
    
    data_to_infer <- data_for_missforest[[1]]
    traits_data_factors <- data_for_missforest[[2]]
    traits_data_num <- data_for_missforest[[3]]
    factor_length <- data_for_missforest[[4]]
    
    
  } # data prepped
  
  
  # Run Missforest
  model_eval_missforest_MCA <- fct_missforest_evaluation_MCA(
    data_to_infer, traits_data_factors = traits_data_factors,
    traits_data_num = traits_data_num, factor_length = factor_length,
    model_iteration = 5, #100
    maxiter_mf=5, #10
    ntree_mf=100,  #100
    prop_NA = 0.2)
  
  # Vérifier les résultats
  if (is.null(model_eval_missforest_MCA) || length(model_eval_missforest_MCA) == 0) {
    warning("Aucun résultat retourné pour dim = ", dim)
    return(NULL)
  } else {
    return(model_eval_missforest_MCA)
  }

  # Extract results
 #results <- extract_model_perf(raw_result = test_dimensionality)
  
  #traits_performance <- results[[1]] |> 
    #dplyr::mutate(dimensions = dim)
  
  #return(traits_performance)
}) #END of lapply on dimension number

# #time of computation:
# time <- c("0 dim 04:46", "2 dim 04:31", "5 dim 03:50", "10 dim 05:44", "25 dim 04:52",
#           "50 dim 06:16", "75 dim 07:35", "100 dim 10:05", "200 dim 13:05")

results <- data.frame()

for (i in 1:18) {
  traits_performance <- extract_model_perf(test_dimensionality[[i]])[[1]]
  traits_performance$iteration <- i 
  results <- bind_rows(results, traits_performance)  # Ajout correct des résultats
}


# Ajout de la colonne dimensions en répétant la liste sur les groupes de 24 variables
results_dimensionality <- results |> 
  dplyr::mutate(iteration = rep(dimensions, each = 20)[seq_len(nrow(results))]) |>
  rename(dimensions = iteration)

save(results_dimensionality, file = here::here("Documents", "Sea_of_immaturity", "outputs","choose_dimensionality_phylogeny_mF2.Rdata"))

resumed_data <- results_dimensionality |>
  dplyr::group_by(variable, dimensions) |>
  dplyr::summarise(median_estimate = median(estimate))|>
  tidyr::drop_na(median_estimate) 

resumed_data$dimensions <- as.factor(resumed_data$dimensions)

## Plot performance against the number of dimensions
ggplot(resumed_data) +
  aes(x= dimensions, y= median_estimate, fill= dimensions)+
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width=0.15, alpha=0.2)+
  xlab("") + ylab("Assessement quality (median of R-squared and accuracy among traits)") +
  theme_minimal() +
  scale_fill_manual(values = grDevices::colorRampPalette(
    RColorBrewer::brewer.pal(9, "OrRd"))(length(unique(resumed_data$dimensions)))) +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 10))
  
ggsave(filename = here::here("Documents","Sea_of_immaturity","figures", "02_c_phylogeny_importance_in_MF2_perf.jpg"),
       plot = last_plot(), width = 10, height =8 )


## Plot performance of prediction for each dimensionality ##
dim = 100
res <- results_dimensionality |> dplyr::filter(dimensions == dim)
# Estimate distributions (boxplot)
boxplot_missforest <- estimates_boxplot(df_estimates = res)
boxplot_missforest
ggsave(filename= here::here("Documents","Sea_of_immaturity","figures", paste0("02_c_Missforest3_performance_boxplot_traits_phylogeny_",dim,".jpg")),
       boxplot_missforest, width = 12, height =8 )

# # Estimate distribution (Histograms)
# hist_missforest <- estimates_histogramm(data = res)
# hist_missforest




##-------------Assess error in MissForest-------------

## number of dimensions chosen 
dim = 100
##

Dims <- paste0("Dim", 1:dim)

## Preping data
phylogeny <- species_traits_filtered |> 
  dplyr::select(species_name, Class, Order, Family, Genus) |> 
  tibble::column_to_rownames(var="species_name") |> 
  as.matrix()

dimensionality <- FactoMineR::MCA(phylogeny,ncp=dim, graph=F) 

phylo_space <- dimensionality[["ind"]][["coord"]]
colnames(phylo_space)<- gsub(" ", "", colnames(phylo_space))

traits_and_phylo <- cbind(phylo_space, species_traits_filtered) 
rownames(traits_and_phylo) <- NULL

data_for_missforest <- preping_data(species_traits_df = traits_and_phylo)

data_to_infer <- data_for_missforest[[1]]
traits_data_factors <- data_for_missforest[[2]]
traits_data_num <- data_for_missforest[[3]]
factor_length <- data_for_missforest[[4]]


## Run Missforest evaluation
model_eval_missforest_MCA <- fct_missforest_evaluation_MCA(
  data_to_infer, traits_data_factors , traits_data_num, factor_length,
  model_iteration = 10, #100
  maxiter_mf=10, #10
  ntree_mf=100,  #100
  prop_NA = 0.2)


save(model_eval_missforest_MCA, file = here::here("Documents","Sea_of_immaturity","outputs", "predictive_model_eval_species_traits3.Rdata"))
# load(file = here::here("outputs", "predictive_model_eval_species_traits.Rdata"))


## Extract results
results <- extract_model_perf(raw_result = model_eval_missforest_MCA)
traits_performance <- results[[1]]
order_performance <- results[[2]]
raw_factors_perf <- results[[3]]
raw_numeric_perf <- results[[4]]


## Plot performance of prediction
# Estimate distributions (boxplot)
boxplot_missforest <- estimates_boxplot(df_estimates = traits_performance)
boxplot_missforest
ggsave(filename = here::here("Documents", "Sea_of_immaturity", "figures", "3_Missforest_final_performance_boxplot_traits.jpg"),
       boxplot_missforest, width = 12, height =8 )

# Estimate distributions (Histograms)
hist_missforest <- estimates_histogramm(data = traits_performance)
hist_missforest
ggsave(filename = here::here("Documents", "Sea_of_immaturity","figures", "3_Missforest_performance_distribution.png"),
       hist_missforest, width = 22, height =14 )


##------------- Infer data on chosen variables-------------

#Choose variable to infer
# var_to_infer <- c("Length", "K", "trophic_guild", "IUCN_inferred_Loiseau23", 
#                   "ClimVuln_SSP585", "Vulnerability", "body_depth_ratio", "body_width_ratio",
#                   "Schooling", "log(depth_min)", "log(depth_max)")
var_to_infer <- c("K", "Length", "Vulnerability_fishing", "Troph", "trophic_guild", "TempPrefM", "TempPrefMean", "TempPrefMin", "IUCN_inferred_Loiseau23", "DemersPelag", "Climate", "Status", "LengthMax_unsexed")

#Number of dimensions chosen 
dim = 100

#Preping data
phylogeny <- species_traits_filtered |> 
  dplyr::select(species_name, Class, Order, Family, Genus) |> 
  tibble::column_to_rownames(var="species_name") |> 
  as.matrix()

dimensionality <- FactoMineR::MCA(phylogeny,ncp=dim, graph=F) 

phylo_space <- dimensionality[["ind"]][["coord"]]
colnames(phylo_space)<- gsub(" ", "", colnames(phylo_space))

traits_and_phylo <- cbind(phylo_space, species_traits_filtered)
rownames(traits_and_phylo) <- NULL

data_for_missforest <- preping_data(species_traits_df = traits_and_phylo)

data_to_infer <- data_for_missforest[[1]] |> dplyr::select(-Class, -Order, -Family, -Genus)
traits_data_factors <- data_for_missforest[[2]]
traits_data_num <- data_for_missforest[[3]]
factor_length <- data_for_missforest[[4]]


#Run several missforest, and keep the most frequent inference if most models 
# converge to the same result
inferred_data <- 
  missforest_applied(
    data_to_infer, 
    factor_length,
    traits_data_factors,
    traits_data_num,
    var_to_infer,
    confidence_threshold = 0.7,
    model_iteration = 10,
    maxiter_mf=10, #missforest parameter
    ntree_mf=100) #missforest parameter


inferred_data_factors <- inferred_data[[1]] %>%
  pivot_wider(names_from = variable, values_from = imputed)
inferred_data_num <- inferred_data[[2]] %>%
  pivot_wider(names_from = variable, values_from = imputed)

inferred_species_traits <- inferred_data_factors  |> 
  dplyr::left_join(inferred_data_num, by = "species_name")

inferred_species_traits <- inferred_data  |> 
  dplyr::left_join( 
    dplyr::select(species_traits_filtered, 
                  species_name)) 

inferred_data <- inferred_data %>%
  rename(species_name = species)
  
#  dplyr::rename(species = "species_name")
# Explore data
fb_plot_species_traits_completeness(inferred_data)
fb_plot_number_species_by_trait(inferred_data)
ggsave(plot= last_plot(), width = 8, height = 8, 
       file= here::here("Documents", "Sea_of_immaturity", "figures","MF3_percent_species_INFERRED.png"))

# #check trophic guild inference
# df <- data.frame(troph_fishbase  = data_to_infer$Troph[which(is.na(data_to_infer$trophic_guild))],
#                  trophic_guild_inferred = inferred_data$trophic_guild[which(is.na(data_to_infer$trophic_guild))])
# 
# ggplot(data=df, aes(x=troph_fishbase, group=trophic_guild_inferred, fill=trophic_guild_inferred)) +
#   geom_histogram(aes(y = ..density..), bins = 20, color = "grey40", fill ="white") +
#   geom_density(aes(fill = trophic_guild_inferred), alpha = 0.2) +
#   hrbrthemes::theme_ipsum() +
#   facet_wrap(~trophic_guild_inferred, scales = "free_y") +
#   theme(
#     legend.position="none",
#     panel.spacing = unit(0.1, "lines"),
#     axis.ticks.x=element_blank()
#   )
# ggsave(filename = here::here("figures", "1_trophic_guild_inference_VS_troph.jpg"),
#        width = 12, height = 8)



##-------------save data-------------

save(inferred_data, file= here::here("Documents", "Sea_of_immaturity", "outputs", "MF3_inferred_species_traits.Rdata"))
# load(file= here::here("outputs", "species_traits_inferred.Rdata"))
write.csv(inferred_species_traits, file= here::here("outputs", "species_traits_inferred.csv"))

##------------Compare data before and after impute--------

species_traits_imputed <- species_traits_filtered %>%
  left_join(inferred_data, by = "species_name") %>%
  mutate(across(matches("\\.x$"), ~ coalesce(.x, get(sub(".x$", ".y", cur_column()))))) %>%
  select(-ends_with(".y")) %>%
  rename_with(~ sub("\\.x$", "", .), ends_with(".x"))


species_traits_imputed <- species_traits_imputed |>
  dplyr::rename(species = "species_name")

fb_plot_number_species_by_trait(species_traits_imputed, threshold_species_proportion = 1)

save(species_traits_imputed, file= here::here("Documents", "Sea_of_immaturity", "outputs", "MF3_species_traits_imputed.Rdata"))
ggsave(plot= last_plot(), width = 8, height = 8, 
       file= here::here("Documents", "Sea_of_immaturity", "figures","1_percent_species_per_trait_afterinference.png"))
