################################################################################
##
##  Takes the Species x traits matrices from fishbase and other sources,
##   check all the species names, and merge the tables.
##
## 1c_merge_species_traits.R
##
## 23/10/2022
##
## Ulysse Flandrin
##
################################################################################
#----------------- cleaning memory -----------------
rm(list=ls())

##-----------------Loading packages-------------------
pkgs <- c("here", "dplyr", "gtools")
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
# ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

# remotes::install_github("FRBCesab/funbiogeo")
##-------------loading data and functions-------------

source(here::here("Documents", "Sea_of_immaturity", "R", "01_a_check_scientific_names.R"))

load(file = here::here("Documents", "Sea_of_immaturity", "data", "derived_data", "01_b_species_traits_fishbase.Rdata"))
load(file = here::here("Documents", "Sea_of_immaturity", "data", "derived_data", "01_c_other_traits_raw.Rdata"))


##-------------Deal with scientific names-------------

species_with_spec_code <- other_traits_raw[!is.na(other_traits_raw$spec_code),]
species_to_find <- other_traits_raw[is.na(other_traits_raw$spec_code),] |> 
  dplyr::select(-c(spec_code, worms_id, fishbase_name))

# /!\ long time to run, and large RAM is needed, reduce mc.core if necessary.
#  put mc_cores = 1 if problem with fishbase
subset1 <- code_sp_check(species_to_find[c(1:5000),], mc_cores = 12) 
subset2 <- code_sp_check(species_to_find[c(5001:10000),], mc_cores = 12) 
subset3 <- code_sp_check(species_to_find[c(10001:nrow(species_to_find)),], mc_cores = 12) 

all_subset <- rbind(subset1, subset2, subset3)
summary(all_subset$check) #should be only 1
NA_fishbase <- all_subset[is.na(all_subset$fishbase_name),] #Mostly non fishes -> OK

other_traits_large_matrix <- all_subset |> 
  dplyr::select(-check) |> 
  gtools::smartbind(species_with_spec_code) 

save(other_traits_large_matrix, file = here::here("data", "fishbase_data", "V0_other_traits_large_matrix.Rdata") )



## Complete data for the same species (same fishbase name)
other_traits_with_SpecCode <- other_traits_large_matrix |> 
  dplyr::group_by(fishbase_name) |> 
  tidyr::fill(tidyr::everything(), .direction = 'updown') |> #fill identical lines
  # dplyr::mutate(across(.cols = where(is.numeric), .fns = mean, .names = "{.col}")) |> #some rows are became duplicates -> means them.
  dplyr::ungroup() |> 
  dplyr::select( -range_n_cells_005, -worms_id) |> #-ClimVuln_SSP126,
  dplyr::distinct() |> 
  dplyr::filter(!is.na(fishbase_name))

save(other_traits_with_SpecCode, file = here::here("data", "fishbase_data", "V0_other_traits_with_spec_code.Rdata") )
# load( file =here::here("data", "fishbase_data", "V0_other_traits_with_spec_code.Rdata") )



## Check duplicates in names
dup_names <- other_traits_with_SpecCode[duplicated(other_traits_with_SpecCode$fishbase_name) |
               duplicated(other_traits_with_SpecCode$fishbase_name, fromLast =T),]
#SOME DUPLICATED SPECIES HAVE SLIGHTLY DIFFERENT VALUES IN TRAITS

other_traits <- other_traits_with_SpecCode |> 
  dplyr::group_by(fishbase_name) |> 
  # Mean the duplicated numeric traits
  dplyr::mutate(across(.cols = where(is.numeric), .fns = ~mean(., na.rm = TRUE), .names = "{.col}")) |> # Mean numeric variables
 
  #Arrange IUCN categories: keep only the highest IUCN status
  dplyr::mutate(IUCN_inferred_Loiseau23 = 
                  ifelse(IUCN_inferred_Loiseau23 == "No Status", 0,
                         ifelse(IUCN_inferred_Loiseau23 == "Non Threatened", 1,
                                ifelse(IUCN_inferred_Loiseau23 == "Threatened", 2, NA)))) |>
  # dplyr::mutate(IUCN_inferred_Loiseau23 = max(IUCN_inferred_Loiseau23, na.rm = TRUE)) |> 
  dplyr::mutate(IUCN_inferred_Loiseau23 = if(all(is.na(IUCN_inferred_Loiseau23))){
    NA} else max(IUCN_inferred_Loiseau23, na.rm = TRUE)) |> 
  dplyr::mutate(IUCN_inferred_Loiseau23 = 
                  ifelse(IUCN_inferred_Loiseau23 == 0, "No Status",
                         ifelse(IUCN_inferred_Loiseau23 == 1, "Non Threatened",
                                ifelse(IUCN_inferred_Loiseau23 == 2, "Threatened", NA)))) |>
  

    dplyr::distinct(across(-species_name), .keep_all = TRUE)


other_traits[duplicated(other_traits$fishbase_name) | duplicated(other_traits$fishbase_name, fromLast =T),]
#Still problem with diet
other_traits$trophic_guild[other_traits$fishbase_name == "Brachygenys chrysargyreum"] <- "microinvertivore"


## Final data for other traits
other_traits <- other_traits |> 
  dplyr::mutate_all(~ifelse(. == "NaN", NA, .))|> #Replace NaN by NA
  dplyr::select(-species_name) |> 
  dplyr::distinct()
  
                                                        
  


##-------------Merge datasets-------------
species_traits_final <- species_traits_fishbase |> 
  dplyr::left_join(other_traits) |> 
  dplyr::select(-UsedforAquaculture, -head_depth_ratio, -lower_jaw_ratio,
                -caudalpeduncle_depth_ratio, -caudalpeduncle_width_ratio,
                -mouth_width_ratio,
                -range_n_cells_01)


##-------------Check and observe data-------------
## IUCN data: IUCN categories used to complete inference by Loiseau et al. (2023)
## mainly for Elasmobranch (that are not taken in Loiseau's study)
iucn <- species_traits_final |> 
  tibble::column_to_rownames("species_name") |> 
  dplyr::select(Class, fishbase_name,
                iucn_inferred = IUCN_inferred_Loiseau23, 
                iucn_redlist = IUCN_category) 

missing_values <- dplyr::filter(iucn, is.na(iucn$iucn_inferred)) |>
  dplyr::mutate(iucn_redlist = dplyr::recode(iucn_redlist,  
                                             "NE" = "No Status",
                                             "DD" = "No Status",
                                             "LC" = "Non Threatened",
                                             "NT" = "Non Threatened",
                                             "VU" = "Threatened",
                                             "EN" = "Threatened",
                                             "CR" = "Threatened")) 

iucn[rownames(missing_values), "iucn_inferred"] <- missing_values$iucn_redlist
iucn$iucn_inferred[iucn$iucn_inferred == "No Status"] <- NA

species_traits_final[, "IUCN_inferred_Loiseau23"] <- iucn[species_traits_final$species_name, 
                                                          "iucn_inferred"]



#Check names
species_traits_final[is.na(species_traits_final$fishbase_name),] #no NA ?
species_traits_final[which(grepl("Genus Species",species_traits_final$fishbase_name)),] #OK
summary(species_traits_final$spec_code) #ok


length(unique(species_traits_final$fishbase_name))# 632 -> 6 DUPLICATED NAMES IN RLS
dup_names<- species_traits_final[duplicated(species_traits_final$fishbase_name) |
                            duplicated(species_traits_final$fishbase_name, fromLast =T),]
# FEW NAs IN DUPLICATED NAMES ARE REMAINING IN FISHBASE DATA -> FILL THESE USELESS GAPS:
dup_names <- dup_names |> 
  dplyr::group_by(fishbase_name) |> 
  tidyr::fill(IUCN_category:last_col(), .direction = "downup") |> 
  dplyr::distinct(across(-species_name), .keep_all = TRUE)

#Same on the full dataframe:
species_traits_final <- species_traits_final |> 
  dplyr::group_by(fishbase_name) |> 
  tidyr::fill(IUCN_category:last_col(), .direction = "downup") |> 
  dplyr::select(-trophic_guild)


## Explore data ##
library(funbiogeo)
library(ggplot2)
species_traits <- dplyr::rename(species_traits_final,species = species_name)

fb_plot_species_traits_completeness(species_traits)
ggsave(plot= last_plot(), file= here::here("figures", "traits_completedness.png"), width = 15, height = 7)
fb_plot_number_species_by_trait(species_traits, threshold_species_proportion = 1)
ggsave(plot= last_plot(), file= here::here("figures", "percent_species_per_traits.png"), width = 8, height = 8)
fb_plot_trait_correlation(species_traits|>dplyr::select(species, where(is.numeric)))
ggsave(plot= last_plot(), file= here::here("figures", "traits_correlogram.png"), width = 10, height = 10)


##-------------Save final data-------------
save(species_traits_final, file = here::here("data", "derived_data", "1c_species_traits_final.Rdata"))
# load( file = here::here("data", "derived_data", "1c_species_traits_final.Rdata"))
save(species_traits_final, file = here::here("outputs", "species_traits_BEFORE_MF.Rdata"))
