
library(mia)
library(dplyr)

tse <- readRDS("tse_greengenes2_MGS.rds")
# remove shallow samples
tse <- tse[, colSums(assay(tse,"counts")) >= 50000]
# New variables
# AR vs true control
# AES vs true control
# union vs true control
# AR: KY56ever == 2
# AES: KY57ever == 2
# union: KY57ever == 2 | KY57ever == 2
# control: KY56ever == 1 & KY57ever == 1
colData(tse)$AR <- case_when(tse$KY56ever == 2 ~ 2, 
                             tse$KY56ever == 1 & tse$KY57ever == 1 ~ 1,
                             TRUE ~ NA) 
colData(tse)$AES <- case_when(tse$KY57ever == 2 ~ 2,
                              tse$KY56ever == 1 & tse$KY57ever == 1 ~ 1,
                              TRUE ~ NA)
colData(tse)$Allergies <- case_when(tse$KY57ever == 2 | tse$KY56ever == 2 ~ 2,
                                    tse$KY56ever == 1 & tse$KY57ever == 1 ~ 1,
                                    TRUE ~ NA)
#treeSE prevalent
tse <- transformAssay(tse,  method = "relabundance")
# removing plasmids (already done)
tse <- tse[grep("plasmid", rowData(tse)[,"domain"],
                ignore.case = TRUE, invert = TRUE),]
altExp(tse, "prevalent_species") <- subsetByPrevalentTaxa(tse, rank = "species",
                                                          prevalence = 1/100,
                                                          detection =0.1/100,
                                                          assay.type="relabundance",
                                                          include_lowest=FALSE)
## add alpha diversity in the object
# estimate (observed) richness, faith
# agglomerating to species
tse <- agglomerateByRank(tse, rank="species", na.rm = TRUE)
tse <- mia::addAlpha(tse, assay.type = "counts",
                     index = "observed_richness",
                     name="observed", niter=100)
tse <- mia::addAlpha(tse,assay.type = "counts",
                     index = "shannon_diversity",
                     name="shannon", niter=100)
tse <- mia::addAlpha(tse, assay.type = "counts",
                     index = "faith_diversity",
                     name="faith", niter=100)
## add beta diversity in the object
tse <- runMDS(tse,
  FUN = getDissimilarity,
  method = "unifrac",
  name = "unifrac",
  tree = rowTree(tse),
  ntop = nrow(tse),
  assay.type = "counts",
  niter = 100,
  keep_dist = TRUE,
  pseudocount = TRUE)

# Add metacyc
tse_func <- readRDS("tse_humann3_functional.rds")
# Only keep the functions that are aggregated at the pathway level
tse_func <- tse_func[!grepl("\\|", rownames(tse_func)),]

# Define taxonomic abundances as the main TreeSE object
altExp(tse, "Pathways") <- tse_func[, colnames(tse)]

# Remove the UNMAPPED and UNINTEGRATED PWs
altExp(tse, "Pathways") <- altExp(tse, "Pathways")[setdiff(rownames(altExp(tse, 
                                                                           "Pathways")), 
                                                           c("UNMAPPED", 
                                                             "UNINTEGRATED")), ]
# Add relabundance transformation
altExp(tse, "Pathways") <- transformAssay(altExp(tse, "Pathways"), 
                                          assay.type="rel_ab", 
                                          method="relabundance", 
                                          pseudocount=FALSE,
                                          MARGIN = "samples")
# Prevalent pathways
altExp(tse, "prevalent_pathway") <- agglomerateByPrevalence(altExp(tse, "Pathways"),
                                                              assay.type="relabundance",
                                                              prevalence = 10/100,
                                                              include_lowest=TRUE)
# Add KO functions
ko_df <- read.table("genefamilies_ko_uniref90_with_names_cpm.tsv",
                    header = TRUE, sep = "\t",
                    comment.char = "", row.names = 1, check.names=FALSE, 
                    fill = TRUE)

ko <- TreeSummarizedExperiment(
  assays = list(relabundance = as.matrix(ko_df)))

colnames(ko) <- colnames(ko) %>%
  gsub("^\\s+|\\s+$", "", .) %>%  # Remove leading/trailing spaces
  gsub("\\.R1\\.trimmed\\.filtered_Abundance-RPKs$", "", .) %>%  # Remove the unwanted substring
  gsub("12142.", "", .) %>%
  gsub("\\.", "-", .)%>%
  sapply(function(name) {
    if (grepl("^40", name)) {  # Check if the name starts with '40'
      gsub("-", "", name)  # Remove hyphen
    } else {
      name
    }
  })
rownames(ko) <- sapply(rownames(ko), function(name) {
    # Split the name by ":" and take the first part
    shortened_name <- strsplit(name, ":")[[1]][1] #:
    return(shortened_name)
})
subs <- ko[, colnames(tse)]
altExp(tse, "KO") <- subs
altExp(tse, "KO") <- altExp(tse, "KO")[setdiff(rownames(altExp(tse, "KO")), 
                                               c("UNMAPPED", "UNINTEGRATED")), ]

# Prevalent KOs
altExp(tse, "prevalent_KO") <- agglomerateByPrevalence(altExp(tse, "KO"),
                                                         assay.type="relabundance",
                                                         prevalence = 10/100,
                                                         include_lowest=TRUE)

saveRDS(tse, file = "input/tse_greengenes2_MGS_allergy_withFunct.rds")
