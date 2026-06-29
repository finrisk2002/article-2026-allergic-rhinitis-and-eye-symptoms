
# Visualization of KO and pathway results

library(ggplot2)
library(dplyr)
library(gridExtra)

tse <- readRDS("input/tse_greengenes2_MGS_allergy_withFunct.rds")
variables <- c("AR", "AES", "Allergies")

sig_maaslin2_ko <- lapply(variables, function(var) 
  load_maaslin2_results_ko(var))
ar <- read.table(sig_maaslin2_ko[[1]], header = TRUE, sep = ",")
aes <- read.table(sig_maaslin2_ko[[2]], header = TRUE, sep = ",")
both <- read.table(sig_maaslin2_ko[[3]], header = TRUE, sep = ",")

signifs <- ar$feature[ar$qval < 0.05]
signifs <- c(signifs, aes$feature[aes$qval < 0.05])
signifs <- c(signifs, both$feature[both$qval < 0.05])
sf <- unique(signifs)

ar <- subset(ar, feature %in% sf)
aes <- subset(aes, feature %in% sf)
both <- subset(both, feature %in% sf)

ar$group <- "AR"
ar <- ar %>%
  mutate(sig = case_when(
    qval < 0.001 ~ "***",
    qval < 0.01  ~ "**",
    qval < 0.05  ~ "*",
    TRUE ~ ""))
aes$group <- "AES"
aes <- aes %>%
  mutate(sig = case_when(
    qval < 0.001 ~ "***",
    qval < 0.01  ~ "**",
    qval < 0.05  ~ "*",
    TRUE ~ ""))
both$group <- "CA"
both <- both %>%
  mutate(sig = case_when(
    qval < 0.001 ~ "***",
    qval < 0.01  ~ "**",
    qval < 0.05  ~ "*",
    TRUE ~ ""))

# Combine together
all_data <- bind_rows(ar, aes, both)
all_data$coef <- log2(exp(all_data$coef)) #log2FC

pb <- ggplot(all_data, aes(x = group, y = reorder(feature, coef), 
                           fill = coef)) +
  geom_tile(color = "gray10") +
  scale_fill_gradient2(
    low = "white", high = "maroon4", # no negative coefficients
    midpoint = 0) +
  #coord_flip() +
  geom_text(aes(label = sig), color = "black", size = 4) +
  scale_x_discrete(position = "top") +
  labs(fill = expression(log[2] ~ FC), x = "", y = "", 
       subtitle = "KO functions",
       title = "A") +
  theme_classic(base_size = 14) + 
  theme(plot.title = element_text(face = "bold", size = 24))


#path_aes <- read.table(sig_maaslin2_pathway[[2]], header = TRUE, sep = ",")
pathway <- cbind(as.data.frame(t(assay(altExp(tse, 
                                              "prevalent_pathway")[c("PWY-5030: L-histidine degradation III", 
                                                      "RHAMCAT-PWY: L-rhamnose degradation I"),],
                   "relabundance"))), colData(tse)$AES)
colnames(pathway) <- c("PWY-5030", "RHAMCAT-PWY", "AES")
pathway_df <- pathway |>
  pivot_longer(cols = c("PWY-5030", "RHAMCAT-PWY"), 
               names_to = "pathway", 
               values_to = "relab")
pathway_df$relab <- pathway_df$relab*100

subs <- subset(pathway_df, pathway == "RHAMCAT-PWY")
subs <- na.omit(subs)
test <- ggplot(subs, mapping = aes(x = factor(AES), y = relab,
                                   group = factor(AES),
                                   color = factor(AES))) + 
  geom_jitter(alpha = 0.3) +
  scale_color_manual(values = c("1" = "#0571b0",
                                "2" = "deeppink")) +
  stat_summary(fun = mean, na.rm = TRUE,
    geom = "text",
    aes(label = paste0("mean=", round(after_stat(y), 2))),
    vjust = 0,
    color = "black",
    size = 3.5) +
  ylim(0, 3.5) +
  ylab("relative abundance %") + xlab("") +
  labs(title = "B", subtitle = "RHAMCAT-PWY") +
  theme_bw() +
  scale_x_discrete(breaks = c(1,2), 
                    labels = c("No symptoms", "AES")) +
  theme(text = element_text(size = 14),
        plot.title = element_text(face = "bold", size = 24),
        legend.position = "none")

subs <- subset(pathway_df, pathway == "PWY-5030")
subs <- na.omit(subs)

test2 <- ggplot(subs, mapping = aes(x = factor(AES), y = relab,
                                    group = factor(AES), 
                                    color = factor(AES))) + 
  geom_jitter(alpha = 0.3) +
  scale_color_manual(values = c("1" = "#0571b0",
                                "2" = "deeppink")) +
  stat_summary(
    fun = mean, na.rm = TRUE, 
    geom = "text",
    aes(label = paste0("mean=", round(after_stat(y), 2))),
    vjust = 0,
    color = "black",
    size = 3.5) +
  ylim(0, 3.5) +
  ylab("relative abundance %") + xlab("") +
  labs(title = "", subtitle = "PWY-5030") +
  theme_bw() +
  scale_x_discrete(breaks = c(1,2), 
                   labels = c("No symptoms", "AES")) +
  theme(text = element_text(size = 14),
        plot.title = element_text(face = "bold", size = 24),
        legend.position = "none")
  
p0 <- grid::rectGrob(gp = grid::gpar(col = "white"))
pws <- grid.arrange(test, test2, ncol = 2)
pleft <- grid.arrange(pws, p0, ncol = 1, heights = c(0.2, 0.8))

fig3 <- grid.arrange(pb, pleft, widths = c(0.4, 0.6), ncol = 2)
ggsave(fig3,
       file="output/results_final/fig3_draft.png",
       width = 10, height=12, dpi = 300, device = "png", bg = "white")
