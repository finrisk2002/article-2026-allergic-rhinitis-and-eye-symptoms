
library(dplyr)
library(vegan)
library(scater)
library(ggplot2)
library(gridExtra)
library(ggplot2)
library(tidyr)
library(Maaslin2)
library(TreeSummarizedExperiment)


# Gets a subset of object that includes prevalent taxa, pathway and KO
# Extract the tree experiment data
tse <- readRDS("input/tse_greengenes2_MGS_allergy_withFunct.rds")

variables <- c("AR", "AES", "Allergies")
metadata_cat= c("MEN","CURR_SMOKE")
metadata = c("BL_AGE","BMI")

main_colData <- colData(tse)

# Prevalent pathways for DAA
altExp(function_tse, "prevalent_pathway") <- agglomerateByPrevalence(altExp(function_tse, 
                                                                            "Pathways"),
                                                                     assay.type="relabundance",
                                                                     prevalence = 10/100,
                                                                     include_lowest=TRUE)
prevalent_pathway <- altExp(tse, "prevalent_pathway")
colData(prevalent_pathway) <- main_colData
prevalent_KO <- altExp(tse, "prevalent_KO")
colData(prevalent_KO) <- main_colData

prevalent_species <- altExp(tse, "prevalent_species")

# Prepare the data and metadata for MaAsLin2
pathway_df <- as.data.frame(t(assay(prevalent_pathway, "rel_ab")))
ko_df <- as.data.frame(t(assay(prevalent_KO, "relabundance"))) 
species_df <- as.data.frame(t(assay(prevalent_species, "counts")))

meta_funct <- as.data.frame(colData(tse))

# maaslin required short name for analysis
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


# Function to subset data and run maaslin2
run_maaslin2_by_area <- function(variable, data_pathway, data_KO, data_species, meta_funct) {
  
  areas <- c("East" = 1, "West" = 0)
  results <- list()
  meta_funct[c("EAST", metadata_cat)] <- lapply(meta_funct[c("EAST",
                                                             metadata_cat)], 
                                                as.factor)
  
  for (area in names(areas)) {
    area_value <- areas[[area]]
    # Subset metadata
    meta_area <- meta_funct[meta_funct$EAST == area_value, ]
    # Subset data based on samples in metadata
    samples <- rownames(meta_area)

    pathway_area <- data_pathway[samples, ,drop = FALSE]
    ko_area <- data_KO[samples, ,drop = FALSE]
    species_area <- data_species[samples, ,drop = FALSE]
    
    # Convert categorical metadata to factors
    meta_area[metadata_cat] <- lapply(meta_area[metadata_cat], as.factor)
    meta_area[[variable]] <- as.factor(meta_area[[variable]])
    
    # Create the formula
    formula <- paste(paste(metadata, collapse = " + "), "+", 
                     paste(metadata_cat, collapse = " + "), "+", variable)
    
    # Run MaAsLin2 for each dataset
    fit_pathway <- Maaslin2(
      input_data = pathway_area,
      input_metadata = meta_area,
      min_prevalence = 0,
      output = paste0("output/area/", area, "_maaslin2_pathway_", 
                      variable),
      fixed_effects = c(metadata, metadata_cat, variable),
      normalization = "NONE",
      transform = "LOG",
      analysis_method = "LM",
      correction = "holm",
      max_significance = 0.05
    )
    
    fit_KO <- Maaslin2(
      input_data = ko_area,
      input_metadata = meta_area,
      min_prevalence = 0,
      output = paste0("output/area/", area, "_maaslin2_KO_",
                      variable),
      fixed_effects = c(metadata, metadata_cat, variable),
      normalization = "NONE",
      transform = "LOG",
      analysis_method = "LM",
      correction = "holm",
      max_significance = 0.05
    )

    fit_species <- Maaslin2(
      input_data = species_area,
      input_metadata = meta_area,
      output = paste0("output/area/", area, "_maaslin2_species_",
                      variable),
      fixed_effects = c(metadata, metadata_cat, variable),
      normalization = "TSS",   # Total Sum Scaling normalization
      transform = "LOG",       # Log transformation
      analysis_method = "LM",
      correction="holm",
      max_significance = 0.05
    )
    # Store results
    results[[area]] <- list(pathway = fit_pathway, KO = fit_KO, species = fit_species)
  }

}

# Run MaAsLin2 for each variable in both East and West
results_maaslin2_by_area <- lapply(variables, function(var) 
  run_maaslin2_by_area(var, pathway_df, ko_df, species_df, meta_funct))


# # Function to calculate prevalence of variables in East and West
calculate_prevalence <- function(area_value, area_name, variables) {
  # Subset metadata for the specified area
  colData(prevalent_species_tse)$EAST <- as.factor(colData(prevalent_species_tse)$EAST)
  meta_spec<- as.data.frame(colData(prevalent_species_tse))
  meta_area <- meta_spec %>% filter(EAST == area_value)
  print(unique(meta_spec$EAST))
  print(area_name)
  print(dim(meta_spec))
  print(dim(meta_area))
  
  # Calculate prevalence for each variable
  prevalence_results <- lapply(variables, function(var) {
    prevalence_table <- meta_area %>%
      group_by(!!sym(var)) %>%
      summarise(count = n(), .groups = "drop") %>%
      mutate(prevalence = count / sum(count) * 100)

    # Save results
    write.csv(prevalence_table, file = paste0("output/area/prevalence_", 
                                              area_name, "_", var, ".csv"), 
              row.names = FALSE)
    return(prevalence_table)
  })

  names(prevalence_results) <- variables
  return(prevalence_results)
}

# Run prevalence calculation for East and West
prevalence_east <- calculate_prevalence(1, "East", c("AR", "AES", "Allergies"))
prevalence_west <- calculate_prevalence(0, "West", c("AR", "AES", "Allergies"))

