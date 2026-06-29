
library(dplyr)
library(vegan)
library(scater)
library(ggplot2)
library(gridExtra)
library(ggplot2)
library(tidyr)
library(TreeSummarizedExperiment)
library(SummarizedExperiment)
library(Maaslin2)

# Extract the tree experiment data
function_tse <- readRDS("input/tse_greengenes2_MGS_allergy_withFunct.rds")
main_colData <- colData(function_tse)

#Prevalent pathway for DAA
altExp(function_tse, "prevalent_pathway") <- agglomerateByPrevalence(altExp(function_tse, 
                                                                            "Pathways"),
                                                            assay.type="relabundance",
                                                            prevalence = 10/100,
                                                            include_lowest=TRUE)

prevalent_pathway <- altExp(function_tse, "prevalent_pathway")
colData(prevalent_pathway) <- main_colData
prevalent_KO <- altExp(function_tse, "prevalent_KO")
colData(prevalent_KO) <- main_colData

# Prepare the data and metadata for maaslin2
pathway_df <- as.data.frame(t(assay(prevalent_pathway, "rel_ab")))
ko_df <- as.data.frame(t(assay(prevalent_KO, "relabundance"))) 
meta_funct<- as.data.frame(colData(function_tse))

# maaslin requires short names
shorten_colnames <- function(df) {
  names(df) <- sapply(names(df), function(name) {
    # Split the name by ":" and take the first part
    shortened_name <- strsplit(name, ":")[[1]][1]
    return(shortened_name)
  })
  return(df)
}

# Apply to your dataframe
pathway_df <- shorten_colnames(pathway_df)
ko_df <- shorten_colnames(ko_df)

### MaasLin2
variables <- c("AR", "AES", "Allergies")
metadata_cat <- c("EAST", "MEN","CURR_SMOKE")
metadata <- c("BL_AGE","BMI")

# Function to run MaAsLin2 for each variable using count data
run_maaslin2 <- function(variable, data_pathway, data_KO) {
  # Create the formula
  formula <- paste(paste(metadata, collapse = " + "), "+", 
                   paste(metadata_cat, collapse = " + "), "+", variable)
  meta_funct[metadata_cat] <- lapply(meta_funct[metadata_cat], as.factor)
  meta_funct[[variable]] <- as.factor(meta_funct[[variable]])
  # Run MaAsLin2
  fit_data_pathway <- Maaslin2(
    input_data = pathway_df,
    input_metadata = meta_funct,
    min_prevalence = 0,
    output = paste0("output/function_cpm/maaslin2_output_pathway_", variable),
    fixed_effects = c(metadata, metadata_cat, variable),
    normalization = "NONE", # cpm
    transform = "LOG",
    analysis_method = "LM",
    correction="holm",
    max_significance = 0.05
  )
  fit_data_KO <- Maaslin2(
    input_data = ko_df,
    input_metadata = meta_funct,
    min_prevalence = 0,
    output = paste0("output/function_cpm/maaslin2_output_KO_", variable),
    fixed_effects = c(metadata, metadata_cat, variable),
    normalization = "NONE", # cpm
    transform = "LOG",
    analysis_method = "LM",
    correction="holm",
    max_significance = 0.05
  )
  
  return(list(pathway = fit_data_pathway, KO = fit_data_KO))
}

# Run MaAsLin2 for each variable
results_maaslin2 <- lapply(variables, function(var) run_maaslin2(var, pathway_df, ko_df))


# Run retrieve result for each variable
# sig_maaslin2_pathway <- lapply(variables, function(var) load_maaslin2_results_pathway(var))
# sig_maaslin2_ko <- lapply(variables, function(var) load_maaslin2_results_ko(var))
# ar <- read.table(sig_maaslin2_ko[[1]], header = TRUE, sep = ",")
# aes <- read.table(sig_maaslin2_ko[[2]], header = TRUE, sep = ",")
# both <- read.table(sig_maaslin2_ko[[3]], header = TRUE, sep = ",")
