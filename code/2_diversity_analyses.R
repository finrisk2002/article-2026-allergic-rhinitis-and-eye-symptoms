
library(mia)
library(dplyr)
library(vegan)
library(scater)
library(gridExtra)
library(ggsignif)
library(ggplot2)
library(patchwork)
library(tidyr)
library(cowplot)
library(miaViz)
library(stats)
library(scales)

# 1) Alpha diversity for all allergy variables (AES, AR, CA) using
# shannon index, observed richness and faith index. Testing between groups using 
# wilcoxon test
# 2) Beta diversity analysis using unifrac distance (PERMANOVA test) + beta dispersion
# 3) Combining the results into the figure 1

tse <- readRDS("input/tse_greengenes2_MGS_allergy_withFunct.rds")

counts <- colSums(assay(tse, "counts"))
richness <- colData(tse)$observed
cor(counts, richness) # ok

variables <- c("AR", "AES", "Allergies")
names <- c("AR", "AES", "CA")
measurements <- c("observed", "shannon") 

# 1) Alpha diversity
# Convert the variable to a factor
# Subset the data based on the variable, remove NA
plot_diversity <- function(tse, variable, diversity) {
  valid_indices <- complete.cases(tse[[variable]])
  tse_sub <- tse[, valid_indices]
  if (variable == "AR") {
    cols <- c("#0571b0", "magenta")
    name <- "AR"
  } else if (variable == "AES") {
    cols <- c("#0571b0", "deeppink")
    name <- "AES"
  } else {
    cols <- c("#0571b0","maroon4")
    name <- "CA"
  }
  if (diversity == "observed") {
    ylab <- "Richness"
    acc <- 1
  } else if (diversity == "shannon") {
    ylab <- "Shannon"
    acc <- 0.01
  } else {
    ylab <- "Faith"
    acc <- 1
  }
  tse_sub[[variable]] <- factor(tse_sub[[variable]],
                                levels=c("1","2"),
                                labels=c("No symptoms", name))
  var_level <- levels(tse_sub[[variable]])
  comb <- combn(var_level, 2, simplify = FALSE)
  # Create the richness plot
  p <- plotColData(
    tse_sub,
    y = diversity,
    x = variable,
    colour_by = variable,
    jitter_type = "jitter",
    point_size = 1.5,
    point_alpha = 0.4,
    show_violin = FALSE
  ) +
    stat_summary(
      fun = mean,
      geom = "text",
      aes(label = paste0("mean=", scales::number(after_stat(y), 
                                                 accuracy = acc,
                                                 big.mark = ",",
                                                 decimal.mark = "."))),
      vjust = 0,
      color = "black",
      size = 3.5) +
    geom_signif(comparisons = comb, textsize = 3.5,
                test = "wilcox.test",
                map_signif_level = function(p) {
                  paste0("p = ", format.pval(p, digits = 2))
                }) + 
    theme_classic() +
    scale_color_manual(values = cols)+
    ylim(0, 1.2*max(colData(tse_sub)[[diversity]])) + 
    labs(y = ylab) +
    theme(text = element_text(size = 10), legend.position = "none",
          axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14))
  if (diversity == "faith") {
    p <- p +
      scale_y_continuous(labels = label_comma(), limits = c(10000, 30000))
  }
  return(p)
}

# Add y-axis to richness and shannon plots
diversity_plots <- list()
shannon_plots <- list()
faith_plots <- list()
diversity_plots[[1]] <- plot_diversity(tse, "AR", "observed") + 
                            labs(title = "B") + 
                            theme(plot.title = element_text(face = "bold", 
                                                            size = 24))
diversity_plots[[2]] <- plot_diversity(tse, "AES", "observed")
diversity_plots[[3]] <- plot_diversity(tse, "Allergies", "observed")

shannon_plots[[1]] <- plot_diversity(tse, "AR", "shannon") +
                            labs(title = "C") +
                            theme(plot.title = element_text(face = "bold", 
                                                            size = 24))
shannon_plots[[2]] <- plot_diversity(tse, "AES", "shannon")
shannon_plots[[3]] <- plot_diversity(tse, "Allergies", "shannon")

faith_plots[[1]] <- plot_diversity(tse, "AR", "faith") +
                            labs(title = "D") +
                            theme(plot.title = element_text(face = "bold", 
                                                            size = 24))
faith_plots[[2]] <- plot_diversity(tse, "AES", "faith")
faith_plots[[3]] <- plot_diversity(tse, "Allergies", "faith")


# 2) Beta diversity

# PCoA
# Add variables for having both symptoms or neither symptoms
df <- as.data.frame(colData(tse))
# Calculate explained variance
e <- attr(reducedDim(tse, "unifrac"), "eig")
rel_eig <- e / sum(e[e > 0])

# Flip the first coordinate
#reducedDim(tse, "PCoA_BC")[, 1] <- -reducedDim(tse, "PCoA_BC")[, 1]
#reducedDim(tse_transform, "PCoA_BC")[, 2] <- -reducedDim(tse_transform, "PCoA_BC")[, 2]

# Extract PCoA coordinates
pcoa_coords <- as.data.frame(reducedDim(tse, "unifrac"))
colnames(pcoa_coords) <- c("PC1", "PC2")
pcoa_coords$Allergies <- colData(tse)$Allergies

values=c("#0571b0", "maroon4")

# Filter out NA values in the variable column
pcoa_coords <- na.omit(pcoa_coords)

pcoa_plot <- ggplot(pcoa_coords, aes(x = PC1, y = PC2, 
                                     color = as.factor(Allergies)))+
  geom_point(alpha=0.5, size=1.5) +
  theme_classic() +
  # stat_ellipse(inherit.aes = F,
  #       aes(color=as.factor(Allergy), x=PC1, y=PC2)) +
  scale_color_manual(labels = c("No symptoms   ", "CA"), values = values)+
  labs(x = paste("PCo1 (", round(100 * rel_eig[1], 1), "%)", sep = ""),
       y = paste("PCo2 (", round(100 * rel_eig[2], 1), "%)", sep = ""),
       title = "A",
       subtitle = "Unifrac") +
  theme(legend.title= element_blank(),
        axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size=12),
        legend.position='bottom',
        plot.title = element_text(face = "bold", size = 24)) +
  guides(
    color = guide_legend(override.aes = list(size = 3.5, alpha = 1)))


# Unifrac dissimilarity PERMANOVA analysis
# Set seed for reproducibility
set.seed(1576)

variable <- "Allergies"
metadata_cat <- c("EAST", "MEN","CURR_SMOKE")
metadata <- c("BL_AGE","BMI")

# Find indices of rows with complete data in the specified metadata columns
valid_indices <- complete.cases(colData(tse)[, c(metadata, metadata_cat, variable)])

# Subset the tse_transform object to only include rows with complete metadata
tse_rda <- tse[, valid_indices]

dist <- attr(reducedDim(tse, "unifrac"), "dist")
dist_sub <- as.matrix(dist)[valid_indices, valid_indices]

adonis_res <- adonis2(dist(dist_sub) ~ Allergies, 
                      data = colData(tse_rda),
                      by = "margin", permutations = 999)
# Beta dispersion
group <- colData(tse)$Allergies
bd <- betadisper(dist, group)
permutest(bd) # significance not explained by homogeneity condition

# 3) Combine the plots
# Arrange plots in a grid and add labels (A and B)
up <- grid.arrange(diversity_plots[[1]], 
                   shannon_plots[[1]], 
                   faith_plots[[1]],
                   ncol = 3,
                   widths = c(0.17, 0.17, 0.21)
)
mid <- grid.arrange(diversity_plots[[2]], shannon_plots[[2]], 
                    faith_plots[[2]],
                    ncol = 3,
                    widths = c(0.17, 0.17, 0.21)
)
bot <- grid.arrange(diversity_plots[[3]], shannon_plots[[3]], 
                    faith_plots[[3]],
                    ncol = 3,
                    widths = c(0.17, 0.17, 0.21)
)
alphs <- grid.arrange(up, mid, bot, ncol = 1, heights = c(0.2, 0.17, 0.17))
beta <- grid.arrange(pcoa_plot, ncol = 1)

fig1 <- grid.arrange(beta, alphs, ncol = 1, heights = c(0.4, 0.6))

ggsave(fig1,
       file="output/results_final/fig1_check.png",
       width = 10, height=10, dpi = 300, device = "png")
