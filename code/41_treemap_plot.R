
library(mia)
library(ggtree)
library(ggnewscale)
library(miaViz)

tse <- readRDS("input/tse_greengenes2_MGS_allergy_withFunct.rds")
variables <- c("AR", "AES", "Allergies")

sig_maaslin2_species <- lapply(variables, function(var) 
  load_maaslin2_results_species(var))
ar <- read.table(sig_maaslin2_species[[1]], header = TRUE, sep = ",")
aes <- read.table(sig_maaslin2_species[[2]], header = TRUE, sep = ",")
both <- read.table(sig_maaslin2_species[[3]], header = TRUE, sep = ",")

ar$feature <- gsub("\\.", " ", ar$feature)
aes$feature <- gsub("\\.", " ", aes$feature)
both$feature <- gsub("\\.", " ", both$feature)

signifs <- ar$feature[ar$qval < 0.05]
signifs <- c(signifs, aes$feature[aes$qval < 0.05])
signifs <- c(signifs, both$feature[both$qval < 0.05])
signifs <- unique(signifs)

rownames(altExp(tse, "prevalent_species")) <- gsub("-", " ", 
                                                   rownames(altExp(tse, 
                                                         "prevalent_species")))
x <- altExp(tse, "prevalent_species")[
  rownames(altExp(tse, "prevalent_species")) %in% signifs, ]

names <- rownames(x)
ar <- ar[match(names, ar$feature), ]
aes <- aes[match(names, aes$feature), ]
both <- both[match(names, both$feature),]

# make non significant white for plotting
ar$coef <- ifelse(ar$qval < 0.05, ar$coef, 0)
aes$coef <- ifelse(aes$qval < 0.05, aes$coef, 0)
both$coef <- ifelse(both$qval < 0.05, both$coef, 0)

# log2FC
rowData(x)$AR <- log2(exp(ar$coef[ar$feature %in% signifs]))
rowData(x)$AES <- log2(exp(aes$coef[aes$feature %in% signifs]))
rowData(x)$Allergies <- log2(exp(both$coef[both$feature %in% signifs]))

rowData(x)$phylum <- gsub("_", " ", rowData(x)$phylum)
#rownames(x) <- gsub("_", " ", rownames(x))

rownames(x) <- sub("_", "<", rownames(x))
rownames(x) <- gsub("_", " ", rownames(x))
rownames(x) <- sub("<", " _", rownames(x))

rownames(x) <- gsub("^([A-Z])[a-z]+\\s+", "\\1. ", rownames(x))
rownames(x) <- sub("s _", "s ", rownames(x))

# Annotation data.frame
annotation_df <- as.data.frame(
  rowData(x)[, c("AR", "AES", "Allergies")]
)
colnames(annotation_df) <- c("AR", "AES", "CA")

# Base tree
p <- plotRowTree(
  x,
  tip.colour.by = "phylum",
  #tip.size.by   = "log_mean",
  show.label    = FALSE,
  point_size = 3,
  line.width = 0.8
) +
  scale_size_continuous(range = c(2, 6))

p$data$label <- gsub("species:", "", p$data$label)
p$data$label <- gsub("-", " ", p$data$label)
#p$data$label <- gsub("_", " ", p$data$label)

p$data$label <- sub("_", "<", p$data$label)
p$data$label <- gsub("_", " ", p$data$label)
p$data$label <- sub("<", " _", p$data$label)

p$data$label <- gsub("^([A-Z])[a-z]+\\s+", "\\1. ", p$data$label)
p$data$label <- sub("s _", "s ", p$data$label)

# add angle manually as it does align with the heatmap
angles <- seq(from = 10, by = 7.234, length.out = 47)
tip_order <- subset(p$data, isTip)$label[
  order(subset(p$data, isTip)$y, decreasing = TRUE)]
inds <- match(tip_order, p$data$label[1:47])
angles <- angles[inds]

tips <- subset(p$data, isTip)
tips$angle_custom <- angles
p$data$angle[p$data$isTip] <- ifelse(p$data$angle[p$data$isTip] < 200, 
                                     p$data$angle[p$data$isTip], 
                                     p$data$angle[p$data$isTip] - 5)

# Full annotated circular tree
p_annotated <- gheatmap(
  p + new_scale_fill(),
  annotation_df,
  #offset = 0,
  width = 0.4,
  colnames = TRUE,
  colnames_angle = 45,
  colnames_offset_y = 0,
  font.size = 3.5
) +
  layout_circular() +
  ggtree::geom_tiplab(
    size   = 3.5,
    align  = TRUE,
    fontface = "italic",
    offset = 200
  
  ) +
  scale_fill_gradient2(
    low = "#0571b0", mid = "white", high = "maroon4",
    midpoint = 0, name = expression(log[2]~FC),
    position = "bottom", 
    guide = guide_colorbar(
      theme = theme(
        legend.text = element_text(face = "plain")))
  ) +
  guides(colour = guide_legend(override.aes = list(size = 5))
  ) +
  theme(
    #legend.text  = element_text(size = 12),
    #legend.title = element_text(size = 14),
    text = element_text(size = 14),
    legend.text = element_text(size = 12, face = "italic"),
    legend.key.size = unit(0.5, "cm"),
    plot.title = element_text(face = "bold", size = 24),
    legend.position = "left",
    legend.justification = "bottom",
    legend.direction = "vertical",
    plot.margin = unit(c(20, 100, 30, 20), "pt")
)
p_annotated

ggsave(p_annotated,
       file="output/results_final/fig2_tree.png",
       width = 10, height=8, dpi = 300, device = "png", bg = "white")

