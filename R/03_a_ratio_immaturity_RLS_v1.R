pkgs <- c("here", "dplyr", "ggplot2", "funbiogeo", "missForest", "slam",
          "pbmcapply", "patchwork", "ggplot2", "maditr","rnaturalearth",
          "tibble", "stringr", "hrbrthemes", "randomForest", "ranger", "caret", "ggpubr", "pheatmap", "tiblr")
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

library(funbiogeo)
library(ggplot2)

##----------load data-----------------
RLS_raw_data <- read.csv(here::here("Documents", "Sea_of_immaturity", "data", "raw_data", "RLS_raw.csv"))
RLS_Med_data <- RLS_raw_data %>%
  filter(country =="France")
RLS_FrenchPolynesia_data <- RLS_raw_data %>%
  filter(country =="French Polynesia")

load(here::here("Documents", "Sea_of_immaturity", "outputs", "RF1_log_LengthMat_predictions.Rdata"))

##----------Determination ratios v1 Med---------------

RLS_Med_data_merged <- left_join(RLS_Med_data, LengthMatMin_u_predicted, by = "species_name")

# Déterminer si chaque individu est mature
RLS_Med_data_merged <- RLS_Med_data_merged %>%
  mutate(is_mature = size_class >= LengthMatMin_unsexed)

# Calcul du ratio de poissons matures par site et par date du survey
maturity_ratios_Med <- RLS_Med_data_merged %>%
  group_by(site_code, survey_date, location, site_name, latitude, longitude) %>%
  summarise(maturity_ratio = mean(is_mature, na.rm = TRUE)) %>%
  ungroup()

# Affichage du résultat
print(maturity_ratios_Med)

##---------Determination ratios v1 French Polynesia---------------------

RLS_FrenchPolynesia_data<- left_join(RLS_FrenchPolynesia_data, LengthMatMin_u_predicted, by = "species_name")

# Déterminer si chaque individu est mature
RLS_FrenchPolynesia_data <- RLS_FrenchPolynesia_data%>%
  mutate(is_mature = size_class >= LengthMatMin_unsexed)

# Calcul du ratio de poissons matures par site et par date du survey
maturity_ratios_FP <- RLS_FrenchPolynesia_data %>%
  group_by(site_code, survey_date, location, site_name, latitude, longitude) %>%
  summarise(maturity_ratio = mean(is_mature, na.rm = TRUE)) %>%
  ungroup()

# Affichage du résultat
print(maturity_ratios_FP)

## --------------Save maturity ratios ---------------------------------
maturity_ratios <- rbind(maturity_ratios_FP, maturity_ratios_Med)

write.csv(maturity_ratios, file = here::here("Documents", "Sea_of_immaturity", "outputs", "Maturity_ratios_v1.csv"), row.names = FALSE)
