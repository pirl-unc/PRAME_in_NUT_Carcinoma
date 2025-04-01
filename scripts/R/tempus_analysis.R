library(tidyverse)
library(ggplot2)
library(cowplot)
library(gt)
library(colorRamp2)
library(ComplexHeatmap)
library(ggridges)
library(ggforce)
library(ggsignif)

# Custom Theme for plots

custom_settings <- list(
  base_family = "Helvetica",
  base_size = 12,
  plot_title_family = "Helvetica",
  plot_title_size = 16,
  plot_title_face = "bold",
  plot_title_margin = 10,
  subtitle_family = "Helvetica",
  subtitle_size = 12,
  subtitle_face = "plain",
  subtitle_margin = 14,
  strip_text_family = "Helvetica",
  strip_text_size = 12,
  strip_text_face = "plain",
  caption_family = "Helvetica",
  caption_size = 9,
  caption_face = "italic",
  caption_margin = 10,
  axis_text_size = "Helvetica",
  axis_title_family = "Helvetica",
  axis_title_size = 9,
  axis_title_face = "plain",
  axis_title_just = "cc",
  plot_margin = margin(30, 30, 30, 30),
  grid_col = "#E2E2E3",
  grid = TRUE,
  axis_col = "#E2E2E3",
  axis = FALSE,
  ticks = FALSE
)

theme_custom <- function(base_size = 12, ticks = TRUE, ylab = TRUE, xlab = TRUE,
                         legend = TRUE, legend_title = TRUE, grid = FALSE){
  half_line <- base_size/2
  ret <- do.call(theme, custom_settings) + theme(line = element_line(color = "black"),
                                                 axis.ticks = element_line(color = "black"),
                                                 axis.ticks.length = unit(4, "pt"), legend.background = element_blank(),
                                                 legend.key = element_blank(), panel.background = element_blank(),
                                                 panel.border = element_blank(), strip.background = element_blank(),
                                                 plot.background = element_blank(), axis.line = element_line(color = "black"),
                                                 legend.text = element_text(size = rel(0.8)), legend.title = element_text(size = rel(1)))
  if (!ticks) {
    ret <- ret + theme(axis.ticks = element_blank(), axis.ticks.x = element_blank(),
                       axis.ticks.y = element_blank())
  }
  if (!ylab) {
    ret <- ret + theme(axis.title.y = element_blank())
  }
  if (!xlab) {
    ret <- ret + theme(axis.title.x = element_blank())
  }
  if (!legend_title) {
    ret <- ret + theme(legend.title = element_blank())
  }
  if (!legend) {
    ret <- ret + theme(legend.position = "none")
  }
  ret
}

# Load data

df <- read_csv("work/data/tempus_data/tempus_NUT_PRAME_paper_data.csv")
genelist <- colnames(df)[-(1:4)]

# Table 1

desired_order <- c("BRD4-NUTM1", "BRD3-NUTM1", "NSD3-NUTM1", "Other")

df <- df %>%
  mutate(fusion_group = ifelse(marker_name %in% c('BRD4-NUTM1','BRD3-NUTM1','NSD3-NUTM1'), marker_name, 'Other'),
         fusion_group = factor(fusion_group, levels = desired_order))

fusion_diagnosis_nut <- df %>% 
  filter(nut_diagnosis) %>% 
  group_by(fusion_group) %>% 
  summarize(nut_diagnosis_count = n(), .groups = "drop")

fusion_total_counts <- df %>% 
  group_by(fusion_group) %>% 
  summarize(total_samples = n(), .groups = "drop")

fusion_site_summary <- df %>% 
  group_by(fusion_group, primary_site) %>% 
  summarize(primary_site_count = n(), .groups = "drop")

fusion_top_sites <- fusion_site_summary  %>% 
  left_join(fusion_total_counts, by = "fusion_group") %>%  # merge total sample counts
  left_join(fusion_diagnosis_nut, by = "fusion_group") %>%   # merge top diagnosis
  group_by(fusion_group) %>%
  mutate(site_percent = round(primary_site_count / total_samples * 100),
         nc_diagnosis_percent = round(nut_diagnosis_count / total_samples * 100)) %>%
  slice_max(primary_site_count, n = 3) %>%
  ungroup() %>%
  mutate(fusion_group = factor(fusion_group, levels = desired_order))

sites_table <- fusion_top_sites %>%
  arrange(fusion_group) %>%
  mutate(fusion_group = as.character(fusion_group)) %>%
  group_by(fusion_group) %>%
  mutate(
    total_samples = ifelse(row_number() == 1,
                           paste0(total_samples, " (", round((total_samples / nrow(df)) * 100), "%)"),
                           ""),
    nut_diagnosis_percent = ifelse(
      row_number() == 1,
      paste0(nut_diagnosis_count, " (", nc_diagnosis_percent, "%)"), ""),
    site_percent = paste0(primary_site_count, " (", round(site_percent, 1), "%)"),
    fusion_group = ifelse(row_number() == 1, fusion_group, "")) %>%
  ungroup() %>%
  select(fusion_group, total_samples, nut_diagnosis_percent, primary_site, site_percent) %>%
  gt() %>%
  tab_header(
    title = "Table: Cohort Characteristics By NUTM1 Gene Fusion Partner (n=165)",
    subtitle = ""
  ) %>%
  cols_label(
    fusion_group = "Fusion",
    total_samples = "Samples (%)",
    nut_diagnosis_percent = "NC Diagnosis (%)",
    primary_site = "Primary Site",
    site_percent = 'Samples Per Site (%)'
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything()))

#gtsave(sites_table, filename = "work/results/tempus/table1.png")

# Heatmap 3A

rna <- df %>% select(all_of(genelist)) %>% drop_na()
plt_df <- df %>% filter(!is.na(MYC))
top_genes <- names(sort(colMeans(as.matrix(rna)), decreasing = TRUE))[1:20]
hm_df <- rna %>% select(all_of(top_genes))
data_matrix <- t(as.matrix(hm_df))
dm_max <- max(data_matrix)
pal <- colorRampPalette(c("blue", "red"))
pal_cols <- pal(3)
color_scale <- colorRamp2(c(0, dm_max/2, dm_max), c(pal_cols[1], pal_cols[2], pal_cols[3]))

oncogene_list <- c("MYC", "TP53", "KRAS")
row_groups <- factor(ifelse(rownames(data_matrix) %in% oncogene_list, "oncogenes", "CTAs"), levels = c("oncogenes", "CTAs"))
column_order <- gsub("-", "::", c("BRD4-NUTM1", "BRD3-NUTM1", "NSD3-NUTM1", "Other"))

ht <- Heatmap(data_matrix,
              name = "log2(TPM+1)",
              show_row_names = TRUE,
              show_column_names = TRUE,
              col = color_scale,
              row_split = row_groups,
              column_split = factor(gsub("-", "::", plt_df$fusion_group), levels = column_order),
              cluster_column_slices = FALSE,
              cluster_row_slices = FALSE,
              column_title_rot = 30,
              row_names_gp = gpar(fontsize = 10),
              column_names_gp = gpar(fontsize = 10),
              row_title_gp = gpar(fontsize = 12),
              column_title_gp = gpar(fontsize = 12))

#png("work/results/tempus/fig3A_heatmap.png", width = 800, height = 350)
#draw(ht, newpage = TRUE)
#dev.off()

# Density scatter Fig3B

prame_threshold = 3.1

plt_df <- df %>% select(-c(1:4))

fig3B_densityscatter <- ggplot(plt_df, aes(x = NUTM1, y = PRAME, color = PRAME > prame_threshold)) +
  geom_density_2d(color = "black") +  # Ensure density lines are black
  geom_hline(yintercept = prame_threshold, linetype = "dashed", color = "#404040", size = 1) +
  geom_point(size = 2, alpha = 0.7) +
  annotate("text", x = 9.5, y = 5, label = expression(italic(PRAME)^"high"), hjust = -0.1, angle = 90, size = 4.5) +
  annotate("text", x = 9.5, y = -0.2, label = expression(italic(PRAME)^"low"), hjust = -0.1, angle = 90, size = 4.5) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  theme_custom() +
  labs(title = "",
       x = "NUTM1 Expression\nLog2(TPM+1)",
       y = "PRAME Expression\nLog2(TPM+1)") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))

#ggsave("work/results/tempus/fig3B_densityscatter.pdf", fig3B_densityscatter, units = 'in', width = 7, height = 5)


# Boxplot Fig3C

plt_df <- df

MIN_max <- min(plt_df[["PRAME"]], na.rm = TRUE)
MAX_max <- max(plt_df[["PRAME"]], na.rm = TRUE)
tenth_RANGE <- (MAX_max - MIN_max) * 0.1

ft_data <- plt_df %>%
  summarize(
    median_value = signif(median(.data[["PRAME"]], na.rm = T), 2),
    mean_value = signif(mean(.data[["PRAME"]], na.rm = T), 2),
    n_all = n(),
    n_high = sum(.data[["PRAME"]] > prame_threshold, na.rm = T),
    n_low = sum(.data[["PRAME"]] < prame_threshold, na.rm = T),
    .by = fusion_group)

fig3C_box <- plt_df %>%
  ggplot(aes(x = fusion_group, y = PRAME, color = fusion_group)) +
  geom_sina(aes(color = ifelse(PRAME > prame_threshold, "PRAMEhigh", "PRAMElow")),
            shape = 21,
            size = 1, alpha = 0.6
  ) +
  scale_color_manual(values = c("PRAMEhigh" = "red", "PRAMElow"="blue"), labels = c(expression(italic(PRAME)^"high"), expression(italic(PRAME)^"low")), name="PRAME Group")+
  geom_boxplot(
    outlier.shape = NULL,
    alpha = 0.3, size = 0.5, width = 0.3, color = "black"
  ) +
  geom_hline(yintercept = prame_threshold, color = "black", linetype = "dashed") +
  scale_x_discrete(labels = c("BRD4-NUTM1"="BRD4::NUTM1", "BRD3-NUTM1"="BRD3::NUTM1", "NSD3-NUTM1"="NSD3::NUTM1")) +
  geom_label(
    data = ft_data,
    aes(
      x = fusion_group,
      y = MIN_max - tenth_RANGE/2, label = paste0("N=", n_high)
    ),
    color = "red", size = 3,
    show.legend = FALSE) +
  geom_label(
    data = ft_data,
    aes(
      x = fusion_group,
      y = MIN_max - tenth_RANGE*1.5, label = paste0("N=", n_low)
    ),
    color = "blue", size = 3,
    show.legend = FALSE) +
  geom_label(
    data = ft_data,
    aes(
      x = fusion_group,
      y = .data[[paste0("median_value")]],
      label = .data[[paste0("median_value")]]
    ),
    color = "black", size = 3.5,
    show.legend = FALSE
  ) +
  scale_y_continuous(expand = expansion(add = c(tenth_RANGE * 1.2, tenth_RANGE * 0.5))) +
  theme_custom() +
  theme(legend.postition = "none",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  labs(x = "",
       y = "PRAME Log2(TPM+1)")

#ggsave('work/results/tempus/fig3C_boxplot.pdf', fig3C_box, units = 'in', width = 7, height = 5)

# Tissue Expression Fig3D

plt_df <- df %>%
  mutate(primary_site = replace_na(primary_site, "Unknown primary site")) %>%
  mutate(primary_site = recode(primary_site, "Skin of other and unspecified parts of face" = "Skin - Other"))

site_percentages <- plt_df %>%
  group_by(primary_site) %>%
  summarise(
    total_samples = n(),
    samples_above_3_1 = sum(PRAME > 3.1, na.rm = TRUE)
  ) %>%
  mutate(
    percent_above_3_1 = (samples_above_3_1 / total_samples) * 100
  )

plt_df$primary_site <- factor(plt_df$primary_site,
                              levels = site_percentages$primary_site[order(site_percentages$percent_above_3_1, decreasing = TRUE)])
site_percentages$primary_site <- factor(site_percentages$primary_site,
                                        levels = site_percentages$primary_site[order(site_percentages$percent_above_3_1, decreasing = TRUE)])

bar_plot <- ggplot(site_percentages, aes(x = primary_site, y = percent_above_3_1)) +
  geom_bar(stat = "identity", fill = "red2", color = "black") +
  theme_classic() +
  labs(y = "% Samples", x = NULL) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 12)
  )

dot_plot <- ggplot(plt_df) +
  geom_jitter(aes(x = primary_site, y = PRAME, color = PRAME > 3.1),
              width = 0.05, size = 2, alpha = 0.7) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  geom_hline(yintercept = 3.1, linetype = "dashed", color = "#404040", size = 0.5) +
  theme_classic() +
  labs(y = "PRAME Expression\nLog2(TPM+1)", x = NULL) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1)
  ) +
  scale_x_discrete(labels = c("(PNS) Peripheral nervous system, NOS" = "(PNS), NOS"))

grid <- plot_grid(bar_plot, dot_plot, ncol = 1, align = "v", rel_heights = c(1, 3))
#ggsave("work/results/tempus/fig3D_tissue_expression.pdf", grid, units = 'in', width = 8, height = 4)

# Fusion Counts Supplemental Figure S2A

plt_df <- df %>% group_by(marker_name) %>% summarise(count = n())

plt_df <- plt_df |>
  mutate(marker_name = gsub("-", "::", marker_name)) |>
  mutate(percent = count/165*100)

# show count of all fusion partners
count_plot <- ggplot(plt_df, aes(x = reorder(marker_name, count), y = count)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = paste0(count, "")), 
            hjust = -0.1,  # adjust horizontal position as needed
            size = 3) +  # adjust text size as needed
  coord_flip() + 
  labs(x = "Fusion", y = "Number of Samples (n=165)", title = "") +
  theme_minimal() + 
  theme(
    axis.text.y = element_text(size = 10, face = "italic"),
    axis.text.x = element_text(size = 12)
  ) 


# Supplemental Figure S2B

plt_df <- df |>
  mutate(marker_name = gsub("-", "::", marker_name)) |>
  filter(!is.na(PRAME))

MIN_max <- min(plt_df[["PRAME"]], na.rm = TRUE)
MAX_max <- max(plt_df[["PRAME"]], na.rm = TRUE)
tenth_RANGE <- (MAX_max - MIN_max) * 0.1

ft_data <- plt_df |>
  summarize(
    median_value = signif(median(.data[["PRAME"]], na.rm = T), 2),
    mean_value = signif(mean(.data[["PRAME"]], na.rm = T), 2),
    n_all = n(),
    n_high = sum(.data[["PRAME"]] > 3.1, na.rm = T),
    n_low = sum(.data[["PRAME"]] < 3.1, na.rm = T),
    .by = marker_name)
marker_list <- ft_data |> filter(n_high > 0) |> pull(marker_name)

plt_df$marker_name <- factor(plt_df$marker_name,
                             levels = ft_data$marker_name[order(ft_data$n_high, decreasing = TRUE)])

point_plot <- plt_df |>
  filter(marker_name %in% marker_list) |>
  ggplot(aes(x = marker_name, y = PRAME, color = marker_name)) +
  geom_point(aes(color = ifelse(PRAME > 3.1, "PRAMEhigh", "PRAMElow")),
             shape = 21,
             size = 2, alpha = 0.6
  ) +
  scale_color_manual(values = c("PRAMEhigh" = "red", "PRAMElow"="blue"), labels = c(expression(italic(PRAME)^"high"), expression(italic(PRAME)^"low")), name="PRAME Group")+
  geom_hline(yintercept = 3.1, color = "black", linetype = "dashed") +
  scale_x_discrete(labels = c("BRD4-NUTM1"="BRD4::NUTM1", "BRD3-NUTM1"="BRD3::NUTM1", "NSD3-NUTM1"="NSD3::NUTM1")) +
  theme_minimal() +
  theme(legend.postition = "none",
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14)) +
  labs(x = "",
       y = "PRAME Log2(TPM+1)")



