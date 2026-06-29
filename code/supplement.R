
library(mia)
library(ggplot2)
library(dplyr)
library(ggforce)
library(tidyr)

# Supplement

# 1) Volcano plot in the functional predictions analysis
# volcano plot about KO functions

for (i in 1:2) { # KO, pathway
  if (i == 1) {
    df <- read.table("output/results/funct_KO.tsv",
                     header = TRUE, sep = ",")
  } else {
    df <- read.table("output/results/funct_pathway.tsv",
                     header = TRUE, sep = ",")
  }
  df$neglog10pval <- -log10(df$pval)
  df <- df %>%
    mutate(
      significant = case_when(
        qval < 0.05 ~ "qval<0.05",
        TRUE ~ "Not significant"
      )
    )
  p <- ggplot(data=df, aes(x=coef, y=neglog10pval, color = significant)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("qval<0.05" = "deeppink3", 
                                  "Not significant" = "#3C3C3C")) +
    guides(color = guide_legend(override.aes = list(size = 3, shape = 15))) +
    theme_minimal(base_size = 14) +
    ylab("-log10(q)") + xlab("coefficient")
  if (i == 1) {
    p_ko <- p + labs(title = "A", colour = "", subtitle = "KO") +
      theme(plot.title = element_text(size = 24, face = "bold"))
  } else {
    p_pw <- p + labs(title = "B", colour = "", subtitle = "Pathway") +
      theme(plot.title = element_text(size = 24, face = "bold"))
  }
}
fig_v <- grid.arrange(p_ko, p_pw, ncol = 2)
ggsave(fig_v,
       file="output/results_final/fig_volcano.png",
       width = 10, height=4, dpi = 300, device = "png", bg = "white")

# 2) Regional analysis
# Overlap features between all and east or west
## species and KO
for (i in 1:2) {
  if (i == 1) {
    df <- read.table("output/results/species.tsv",
                       header = TRUE, sep = ",") 
    e_df <- read.table("output/results/regional_East_species.tsv",
                         header = TRUE, sep = ",")
    w_df <- read.table("output/results/regional_West_species.tsv",
                         header = TRUE, sep = ",")
  } else {
    df <- read.table("output/results/funct_KO.tsv",
                       header = TRUE, sep = ",") 
    e_df <- read.table("output/results/regional_East_KO.tsv",
                         header = TRUE, sep = ",")
    w_df <- read.table("output/results/regional_West_KO.tsv",
                         header = TRUE, sep = ",")
  }
  # We notice east and west are highly overlapping with the whole dataset but 
  # not that much with each other
  signifs <- df %>% filter(qval < 0.05) %>% select(feature) %>% unique()
  e_sig <- e_df %>% filter(qval < 0.05) %>% select(feature) %>% 
    unique() # two taxa not in signifs
  w_sig <- w_df %>% filter(qval < 0.05) %>% select(feature) %>% 
    unique() # two taxa not in signifs
  
  all_features <- unique(c(signifs$feature, e_sig$feature, w_sig$feature))
  df <- df %>% filter(feature %in% all_features) %>% 
    select(feature, metadata, coef, qval)
  e_df <- e_df %>% filter(feature %in% all_features) %>% 
    select(feature, metadata, coef, qval)
  colnames(e_df) <- c("feature", "metadata","e_coef", "e_qval")
  w_df <- w_df %>% filter(feature %in% all_features) %>% 
    select(feature, metadata, coef, qval)
  colnames(w_df) <- c("feature", "metadata","w_coef", "w_qval")
  
  overlap <- df %>%
    dplyr::left_join(e_df, by = c("feature", "metadata")) %>%
    dplyr::left_join(w_df, by = c("feature", "metadata"))
  
  # check in how many rows these is inconsistent coefficients
  sum(overlap$coef < 0 & overlap$e_coef < 0 & overlap$w_coef < 0) # 90 all negative
  sum(overlap$coef > 0 & overlap$e_coef > 0 & overlap$w_coef > 0) # 62 all positive
  # --> in one taxa there is inconsistent result (nonsignificant)
  # --> all KOs are consistent
  # --> no need to visualize coefficient, all positive
  overlap$feature[!((overlap$coef > 0 & overlap$e_coef > 0 & overlap$w_coef > 0) | 
                      (overlap$coef < 0 & overlap$e_coef < 0 & overlap$w_coef < 0))]
  overlap[c("qval", "e_qval", "w_qval")] <- lapply(overlap[c("qval", "e_qval", 
                                                             "w_qval")], 
                                                   as.numeric)
  
  overlap$cat1 <- ifelse(overlap$qval < 0.05, 1, 0)
  overlap$cat2 <- ifelse(overlap$e_qval < 0.05, 1, 0)
  overlap$cat3 <- ifelse(overlap$w_qval < 0.05, 1, 0)
  
  overlap$pattern <- overlap$cat1 + overlap$cat2*2 + overlap$cat3*3 # to distinguish
  overlap <- overlap[order(overlap$pattern), ]
  overlap$cat1 <- as.factor(overlap$cat1)
  overlap$cat2 <- as.factor(overlap$cat2)
  overlap$cat3 <- as.factor(overlap$cat3)
  
  # plot illustration of overlapping results
  overlap$feature <- gsub("\\.", " ", overlap$feature)
  overlap$feature <- gsub("_", " ", overlap$feature)
  overlap$feature <- gsub("^([^ ]+) ([A-Z]) ", "\\1_\\2 ", overlap$feature)
  
  # make a new dataset only with CA coefficients to show the direcation
  effect <- overlap %>% filter(metadata == "Allergies")
  effect$log2fc <- log2(exp(effect$coef)) #log2FC
  # Plot effect size of CA group in the whole dataset
  effect$metadata[effect$metadata == "Allergies"] <- "Effect"
  
  # Basic plot
  p <- ggplot() + 
    geom_tile(effect, mapping = aes(x = metadata, 
                                    y = feature, fill = log2fc),
              alpha = 0.8, color = "gray30") +
    scale_fill_gradient2(
      low = "#0571b0", mid = "white", high = "maroon4",
      midpoint = 0, 
      name = expression(log[2]~FC),
      position = "top") +
    ggnewscale::new_scale_fill() +
    guides(fill = guide_legend(order = 1)) +
    geom_tile(overlap, mapping = aes(x = metadata, y = reorder(feature, pattern), 
                                     fill = cat1),
              alpha = 0.8, color = "gray30") + 
    scale_fill_discrete(name = "All regions",
                        type = c(
                          "1" = "dodgerblue",
                          "0" = "white"
                        )) +
    ggnewscale::new_scale_fill() +
    guides(fill = guide_legend(order = 2)) +
    geom_tile(overlap, mapping = aes(x = metadata, 
                                     y = reorder(feature, pattern), fill = cat2),
              alpha = 0.5, color = "gray30") +
    scale_fill_discrete(name = "East",
                        type = c(
                          "1" = "deeppink",
                          "0" = "white"
                        )) +
    ggnewscale::new_scale_fill() +
    guides(fill = guide_legend(order = 3)) +
    geom_tile(overlap, mapping = aes(x = metadata, 
                                     y = reorder(feature, pattern), fill = cat3),
              alpha = 0.5, color = "gray30") +
    scale_fill_discrete(name = "West",
                        type = c(
                          "1" = "yellow",
                          "0" = "white"
                        ))
  # Unique additions
  if (i == 1) {
    p1 <- p + labs(title = "A", subtitle = "Significant (q-val<0.05) species", 
                   y = "",
                   x = "") +
      theme_classic(base_size = 13) +
      scale_x_discrete(position = "top", labels = c(AR = "AR",
                                                    AES = "AES",
                                                    Allergies = "CA")) +
      theme(plot.title = element_text(face = "bold", size = 24),
            plot.margin = unit(c(20, 0, 20, 20), "pt"),
            axis.text.y = element_text(face = "italic"))
  } else {
    p2 <- p + labs(title = "B", subtitle = "Significant (q-val<0.05) KOs", 
                   y = "",
                   x = "") +
      theme_classic(base_size = 13) +
      scale_x_discrete(position = "top", labels = c(AR = "AR",
                                                    AES = "AES",
                                                    Allergies = "CA")) +
      theme(plot.title = element_text(face = "bold", size = 24),
            plot.margin = unit(c(20, 0, 20, 20), "pt"))
  }
}
fig_reg <- grid.arrange(p1, p2, ncol = 2, widths = c(0.61, 0.39))
ggsave(fig_reg,
       file="output/results_final/fig_regional.png",
       width = 10, height=12, dpi = 300, device = "png", bg = "white")


# 3) COPD
tse <- readRDS("input/tse_greengenes2_MGS_allergy_withFunct.rds")

colData(tse)$INC_COPD <- ifelse(colData(tse)$COPD == 1 & colData(tse)$COPD_AGEDIFF > 0, 
                                1, 0)
colData(tse)$INC_ASTHMA <- ifelse(colData(tse)$ASTHMA == 1 & colData(tse)$ASTHMA_AGEDIFF > 0, 
                                  1, 0)
sum(colData(tse)$Allergies == 1 & colData(tse)$INC_ASTHMA == 0, na.rm = TRUE)
sum(colData(tse)$Allergies == 1 & colData(tse)$INC_COPD == 0, na.rm = TRUE)

for (i in 1:2) { # asthma and copd
  if (i == 1) {
    df <- data.frame(
      A = c("Yes","Yes","No", "No"),
      B = c("Yes","No","Yes", "No"),
      value = c(314,1991,244,2931))
  } else {
    df <- data.frame(
      A = c("Yes","Yes","No", "No"),
      B = c("Yes","No","Yes", "No"),
      value = c(85,2220,100,3075))
  }
  df <- df %>% mutate(id = interaction(A,B))
  df_long <- df %>%
    pivot_longer(cols = c(A, B), names_to = "x", values_to = "category")
  p <- ggplot(df_long, aes(x = x,
                            id = id,
                            split = category,
                            value = value))
  if (i == 1) {
    p_asthma <- p + 
      scale_x_discrete(labels = c("Baseline allergies", "Future asthma")) +
      xlab("") +
      labs(title = "A") +
      geom_parallel_sets(aes(fill = id), alpha = 0.7) +
      geom_parallel_sets_axes() +
      geom_parallel_sets_labels(angle = 0) + theme_classic(base_size = 12) +
      theme(legend.position = "none", text = element_text(size = 16),
            title = element_text(face = "bold", size = 24))
  } else {
    p_copd <- p + 
      scale_x_discrete(labels = c("Baseline allergies", "Future COPD")) +
      xlab("") +
      labs(title = "B") +
      geom_parallel_sets(aes(fill = id), alpha = 0.7) +
      geom_parallel_sets_axes() +
      geom_parallel_sets_labels(angle = 0) + theme_classic(base_size = 12) +
      theme(legend.position = "none", text = element_text(size = 16),
            title = element_text(face = "bold", size = 24))
  }
}

supp_sankey <- grid.arrange(p_asthma, p_copd, widths = c(0.5, 0.5))

ggsave(supp_sankey,
       file="output/results_final/supp_sankey.png",
       width = 10, height=4, dpi = 300, device = "png", bg = "white")
