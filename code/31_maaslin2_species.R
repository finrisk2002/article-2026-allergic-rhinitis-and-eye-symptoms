
library(dplyr)
library(vegan)
library(scater)
library(ggplot2)
library(gridExtra)
library(ggplot2)
library(tidyr)
library(SummarizedExperiment)
library(Maaslin2)

tse <- readRDS("input/tse_greengenes2_MGS_allergy_withFunct.rds")

altExp(tse, "prevalent_species") <- agglomerateByPrevalence(tse, rank = "species",
                                                          prevalence = 1/100,
                                                          detection =0.1/100,
                                                          assay.type="relabundance",
                                                          include_lowest=FALSE)

# Prepare the data and metadata for MaAsLin2
spec_df <- as.data.frame(t(assay(altExp(tse, "prevalent_species"), 
                                 "counts")))  # Using count data
meta_spec <- as.data.frame(colData(altExp(tse, "prevalent_species")))

### MaasLin2
variables <- c("AR", "AES", "Allergies")
metadata_cat= c("EAST","MEN","CURR_SMOKE")
metadata = c("BL_AGE","BMI")

# Function to run MaAsLin2 for each variable using count data
run_maaslin2 <- function(variable, data_spec) {
  # Create the formula
  formula <- paste(paste(metadata, collapse = " + "), "+", paste(metadata_cat, 
                                               collapse = " + "), "+", variable)
  meta_spec[metadata_cat] <- lapply(meta_spec[metadata_cat], as.factor)
  meta_spec[[variable]] <- as.factor(meta_spec[[variable]])
  # Run MaAsLin2
  fit_data_spec <- Maaslin2(
    input_data = data_spec,
    input_metadata = meta_spec,
    output = paste0("output/species/maaslin2_output_spec_", 
                    variable),
    fixed_effects = c(metadata, metadata_cat, variable),
    normalization = "TSS",   # Total Sum Scaling normalization
    transform = "LOG",       # Log transformation
    analysis_method = "LM",
    correction="holm",
    max_significance = 0.05
  )
  return(list(species = fit_data_spec))
}

# # Run MaAsLin2 for each variable
results_maaslin2 <- lapply(variables, function(var) run_maaslin2(var, spec_df))

# # Run retrieve result for each variable
# sig_maaslin2_species <- lapply(variables, function(var) load_maaslin2_results_species(var))
# ar <- read.table(sig_maaslin2_species[[1]], header = TRUE, sep = ",")
# aes <- read.table(sig_maaslin2_species[[2]], header = TRUE, sep = ",")
# both <- read.table(sig_maaslin2_species[[3]], header = TRUE, sep = ",")
