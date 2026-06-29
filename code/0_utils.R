
# # 1) Maaslin2 results
# species
load_maaslin2_results_species <- function(variable) {
  results <- read.csv(paste0("output/species/maaslin2_output_spec_",
                             variable, "/all_results.tsv"), sep = "\t")
  # Filter results where metadata == variable
  filtered_results <- subset(results, metadata == variable)
  # Sort filtered results based on qval (ascending)
  sorted_results <- filtered_results[order(filtered_results$qval), ]
  output_file <- paste0("output/species/maaslin2_output_spec_",
                        variable, "/all_results_",variable,".csv")
  write.csv(sorted_results, file = output_file, row.names = FALSE)
  return(output_file)
}

# metacyc pathways
load_maaslin2_results_pathway <- function(variable) {
  results <- read.csv(paste0("output/function_cpm/maaslin2_output_pathway_",
                             variable, "/significant_results.tsv"), sep = "\t")
  # Filter results where metadata == variable
  filtered_results <- subset(results, metadata == variable)
  # Sort filtered results based on qval (ascending)
  sorted_results <- filtered_results[order(filtered_results$qval), ]
  output_file <- paste0("output/function_cpm/maaslin2_output_pathway_",
                        variable, "/significant_results_",variable,".csv")
  write.csv(sorted_results, file = output_file, row.names = FALSE)
  return(output_file)
}

# Function to load MaAsLin2 results
load_maaslin2_results_ko <- function(variable) {
  results <- read.csv(paste0("output/function_cpm/maaslin2_output_KO_",
                             variable, "/all_results.tsv"), sep = "\t")
  # Filter results where metadata == variable
  filtered_results <- subset(results, metadata == variable)
  # Sort filtered results based on qval (ascending)
  sorted_results <- filtered_results[order(filtered_results$qval), ]
  output_file <- paste0("output/function_cpm/maaslin2_output_KO_",
                        variable, "/all_results_",variable,".csv")
  write.csv(sorted_results, file = output_file, row.names = FALSE)
  return(output_file)
}

# Regional
load_maaslin2_results <- function(area, analysis, variable, path) {
  results <- read.csv(path, sep = "\t")
  # Filter results where metadata == variable
  filtered_results <- subset(results, metadata == variable)
  # Sort filtered results based on qval (ascending)
  sorted_results <- filtered_results[order(filtered_results$qval), ]
  output_file <- paste0("output/area/", area, "_", 
                        analysis, "_significant_results_",variable,".csv")
  write.csv(sorted_results, file = output_file, row.names = FALSE)
  return(output_file)
}