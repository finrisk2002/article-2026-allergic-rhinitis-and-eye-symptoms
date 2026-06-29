
library(tidyr)
library(dplyr)
library(data.table)
library(mia)

tse <- readRDS("input/tse_greengenes2_MGS_allergy_withFunct.rds")
table_data_ori <- as.data.frame(colData(tse))
setDT(table_data_ori) 

var_num = c("BL_AGE", "BMI","SMOKING_PACK_YEARS")
var_char = c("MEN", "CURR_SMOKE","EAST")

output_datadir1 = "output/clin_stats_AR.tsv"
output_datadir2 = "output/clin_stats_AES.tsv"
output_datadir3 = "output/clin_stats_Allergies.tsv"

# 1) Summary table of all allergies
selected_cols <- c(var_num, var_char, "AR")
table_data <- table_data_ori[, .SD, .SDcols = selected_cols]

table1 <- transpose(merge(rbind(table_data[, c(AR = "all", 
                 lapply(.SD, function(x) paste0(round(mean(x, na.rm=T),2), "+-", 
                                                round(sd(x, na.rm=T),2),
                                                "(",round(min(x, na.rm=T),2),"-",
                                                round(max(x, na.rm=T),2),")")) ), 
                 .SDcols=var_num],
                  table_data[, lapply(.SD, function(x) paste0(round(mean(x, na.rm=T),2), "+-", 
                                                round(sd(x, na.rm=T),2),
                                                "(",round(min(x, na.rm=T),2),"-",
                                                round(max(x, na.rm=T),2),")")), 
                             by=AR, .SDcols=var_num]),
                 rbind(table_data[, c(AR = "all", 
                 lapply(.SD, function(x) paste0(table(x)[2], " (", 
                                                round(table(x)[2]/length(x)*100, 1), ")")) ), 
                 .SDcols=var_char],
                 table_data[, lapply(.SD, function(x) paste0(table(x)[2], " (", 
                                                 round(table(x)[2]/length(x)*100, 1), ")")), 
                            by=AR, .SDcols=var_char]), 
                 by = "AR"), keep.names = "col", make.names = "AR")

table1 <- table1[,c("col", "all", "1", "2")]
p_values_cont <- as.vector(as.matrix(table_data[, 
                                                lapply(.SD, function(x) 
                                                  signif(wilcox.test(x ~ AR)$p.value, 
                                                                   digits = 2)), 
                                                .SDcols = var_num]))
for (i in 1:length(var_char)) {
  p_values_cont[length(var_num)+i] <- signif(fisher.test(rbind(table(table_data[table_data$AR == 1,][[var_char[i]]]), 
                                                               table(table_data[table_data$AR == 2, ][[var_char[i]]])))$p.value, 
                                             digits = 2)
}
table1$p_value <- p_values_cont
table1 <- rbind(list("N", nrow(table_data), nrow(table_data[table_data$AR == 1,]), 
                     nrow(table_data[table_data$AR == 2,])), 
                table1, fill=T)
colnames(table1) <- c("Variable", "Total", "No Allergic Rhinitis Symptoms", 
                      "With Allergic Rhinitis Symptoms", "p-value")
write.table(table1, output_datadir1, sep="\t" ,col.names=NA)

selected_cols <- c(var_num, var_char, "AES")
table_data <- table_data_ori[, .SD, .SDcols = selected_cols]
table1 <- transpose(merge(rbind(table_data[, c(AES = "all", 
                                               lapply(.SD, function(x) paste0(round(mean(x, na.rm=T),2), "+-", 
                                                                              round(sd(x, na.rm=T),2),
                                                                              "(",round(min(x, na.rm=T),2),
                                                                              "-",round(max(x, na.rm=T),2),")")) ), 
                                           .SDcols=var_num],
                                table_data[, lapply(.SD, function(x) paste0(round(mean(x, na.rm=T),2), "+-", 
                                                                            round(sd(x, na.rm=T),2),
                                                                            "(",round(min(x, na.rm=T),2),"-",
                                                                            round(max(x, na.rm=T),2),")")), 
                                           by=AES, .SDcols=var_num]),
                          rbind(table_data[, c(AES = "all", lapply(.SD, 
                                                                   function(x) paste0(table(x)[2], " (",
                                                                      round(table(x)[2]/length(x)*100, 1), ")")) ), 
                                           .SDcols=var_char],
                                table_data[, lapply(.SD, function(x) paste0(table(x)[2], " (", 
                                                                      round(table(x)[2]/length(x)*100, 1), ")")), 
                                           by=AES, .SDcols=var_char]), 
                          by = "AES"), keep.names = "col", 
                    make.names = "AES")

table1 <- table1[,c("col", "all", "1", "2")]
p_values_cont <- as.vector(as.matrix(table_data[, lapply(.SD, function(x) signif(wilcox.test(x ~ AES)$p.value, 
                                                                                 digits = 2)), 
                                                .SDcols = var_num]))
for (i in 1:length(var_char)) {
  p_values_cont[length(var_num)+i] <- signif(fisher.test(rbind(table(table_data[table_data$AES == 1,][[var_char[i]]]), 
                                                               table(table_data[table_data$AES == 2, ][[var_char[i]]])))$p.value, 
                                             digits = 2)
}
table1$p_value <- p_values_cont
table1 <- rbind(list("N", nrow(table_data), nrow(table_data[table_data$AES == 1,]), 
                     nrow(table_data[table_data$AES == 2,])), table1, fill=T)
colnames(table1) <- c("Variable", "Total", "No Allergic Eye Symptoms", 
                      "With Allergic Eye Symptoms", "p-value")

#write.csv(table1, output_datadir, row.names=F)
write.table(table1,output_datadir2, sep="\t" ,col.names=NA)


selected_cols <- c(var_num, var_char, "Allergies")
table_data <- table_data_ori[, .SD, .SDcols = selected_cols]

table1 <- transpose(merge(rbind(table_data[, c(Allergies = "all", 
                 lapply(.SD, function(x) paste0(round(mean(x, na.rm=T),2), "+-", 
                                                round(sd(x, na.rm=T),2),
                                                "(",round(min(x, na.rm=T),2),"-",
                                                round(max(x, na.rm=T),2),")")) ), 
                 .SDcols=var_num],
                  table_data[, lapply(.SD, function(x) paste0(round(mean(x, na.rm=T),2), "+-", 
                                                round(sd(x, na.rm=T),2),
                                                "(",round(min(x, na.rm=T),2),"-",
                                                round(max(x, na.rm=T),2),")")), 
                             by=Allergies, .SDcols=var_num]),
                 rbind(table_data[, c(Allergies = "all", 
                 lapply(.SD, function(x) paste0(table(x)[2], " (", 
                                                round(table(x)[2]/length(x)*100, 1), ")")) ), 
                 .SDcols=var_char],
                 table_data[, lapply(.SD, function(x) paste0(table(x)[2], " (", 
                                                 round(table(x)[2]/length(x)*100, 1), ")")), 
                            by=Allergies, .SDcols=var_char]), 
                 by = "Allergies"), keep.names = "col", make.names = "Allergies")

table1 <- table1[,c("col", "all", "1", "2")]
p_values_cont <- as.vector(as.matrix(table_data[, 
                                                lapply(.SD, function(x) signif(wilcox.test(x ~ Allergies)$p.value, 
                                                                               digits = 2)), 
                                                .SDcols = var_num]))
for (i in 1:length(var_char)) {
  p_values_cont[length(var_num)+i] <- signif(fisher.test(rbind(table(table_data[table_data$Allergies == 1,][[var_char[i]]]), 
                                                               table(table_data[table_data$Allergies == 2, ][[var_char[i]]])))$p.value, 
                                             digits = 2)
}
table1$p_value <- p_values_cont
table1 <- rbind(list("N", nrow(table_data), nrow(table_data[table_data$Allergies == 1,]), 
                     nrow(table_data[table_data$Allergies == 2,])), 
                table1, fill=T)
colnames(table1) <- c("Variable", "Total", "No Allergies", 
                      "With Allergies", "p-value")

#write.csv(table1, output_datadir, row.names=F)
write.table(table1,output_datadir3, sep="\t" ,col.names=NA)



# 2) Maaslin2 results tables
maketables_regional <- function(analysis, set) {
  results_aes <- read.csv(paste0("output/area/", set, "_maaslin2_", analysis, 
                                 "_AES/all_results.tsv"), 
                          sep = "\t")
  filtered_aes <- subset(results_aes, metadata == "AES")
  results_ar <- read.csv(paste0("output/area/", set, "_maaslin2_", analysis, 
                                "_AR/all_results.tsv"), 
                         sep = "\t")
  filtered_ar <- subset(results_ar, metadata == "AR")
  results_allergies <- read.csv(paste0("output/area/", set, "_maaslin2_", 
                                       analysis, "_Allergies/all_results.tsv"), 
                                sep = "\t")
  filtered_allergies <- subset(results_allergies, metadata == "Allergies")
  filtered_results <- rbind(filtered_aes, filtered_ar, filtered_allergies)
  # Sort filtered results based on qval (ascending)
  sorted_results <- filtered_results[order(filtered_results$qval), ]
  write.csv(sorted_results, file = paste0("output/results/regional_", set, "_", 
                                          analysis, ".tsv"), row.names = FALSE)
}
maketables_regional("KO", "East")
maketables_regional("KO", "West")
maketables_regional("pathway", "East")
maketables_regional("pathway", "West")
maketables_regional("species", "East")
maketables_regional("species", "West")

maketables_functional <- function(analysis) {
  results_aes <- read.csv(paste0("output/function_cpm/maaslin2_output_"
                                 , analysis, "_AES/all_results.tsv"), 
                          sep = "\t")
  filtered_aes <- subset(results_aes, metadata == "AES")
  results_ar <- read.csv(paste0("output/function_cpm/maaslin2_output_", 
                                analysis, "_AR/all_results.tsv"), 
                         sep = "\t")
  filtered_ar <- subset(results_ar, metadata == "AR")
  results_allergies <- read.csv(paste0("output/function_cpm/maaslin2_output_"
                                       , analysis, "_Allergies/all_results.tsv"), 
                                sep = "\t")
  filtered_allergies <- subset(results_allergies, metadata == "Allergies")
  filtered_results <- rbind(filtered_aes, filtered_ar, filtered_allergies)
  # Sort filtered results based on qval (ascending)
  sorted_results <- filtered_results[order(filtered_results$qval), ]
  write.csv(sorted_results, file = paste0("output/results/funct_",
                                          analysis, ".tsv"), 
            row.names = FALSE)
}
maketables_functional("KO")
maketables_functional("pathway")


maketables_species <- function() {
  results_aes <- read.csv(paste0("output/species/maaslin2_output_spec_AES/all_results.tsv"), 
                          sep = "\t")
  filtered_aes <- subset(results_aes, metadata == "AES")
  results_ar <- read.csv(paste0("output/species/maaslin2_output_spec_AR/all_results.tsv"), 
                         sep = "\t")
  filtered_ar <- subset(results_ar, metadata == "AR")
  results_allergies <- read.csv(paste0("output/species/maaslin2_output_spec_Allergies/all_results.tsv"), 
                                sep = "\t")
  filtered_allergies <- subset(results_allergies, metadata == "Allergies")
  filtered_results <- rbind(filtered_aes, filtered_ar, filtered_allergies)
  # Sort filtered results based on qval (ascending)
  sorted_results <- filtered_results[order(filtered_results$qval), ]
  write.csv(sorted_results, file = paste0("output/results/species.tsv"), 
            row.names = FALSE)
}
maketables_species()

# Save definitions tables
ko_definitions <- rownames(altExp(tse, "prevalent_KO"))
ko_definitions <- sub("\\[EC.*", "", ko_definitions)
ko_definitions <- ko_definitions[startsWith(ko_definitions, "K")]
ko_definitions <- do.call(rbind, 
                          strsplit(as.character(ko_definitions), 
                                   split = ":"))
ko_definitions <- ko_definitions[,-3]
ko_definitions[,2] <- trimws(ko_definitions[,2])
write.csv(ko_definitions, file = "output/results_final/ko_definitions.tsv")

pw_definitions <- rownames(altExp(tse, "Pathways"))
pw_definitions <- do.call(rbind, 
                          strsplit(as.character(pw_definitions), 
                                   split = ":"))
pw_definitions <- pw_definitions[,-3]
write.csv(pw_definitions, file = "output/results_final/pw_definitions.tsv")
