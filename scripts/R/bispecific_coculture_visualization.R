library(dplyr)
library(tidyr)
library(drc)
library(here)
library(ggplot2)


# read in data ------------------------------------------------------------

t2_data  <- read.csv('work/data/bispecific_cell_culture/PRAME_T2_peptide_pulse.txt', check.names = F, sep = '\t')
bispecific_time_data <- read.csv('work/data/bispecific_cell_culture/nut_cell_lines_bispecific_luminescence.csv', check.names = F)
conc <- read.csv('work/data/bispecific_cell_culture/PRAME_cytotoxicity_concentration_gradient.csv', check.names = F, sep=',')


# Panel B Peptide required for T Cell Killing  ---------------------------
u_avg <- t2_data %>%
  filter(well_id == 'U') %>%
  summarize(average = mean(raw_val)) %>%
  pull(average)

t_avg <- t2_data %>%
  filter(well_id == 'T') %>%
  summarize(average = mean(raw_val)) %>% 
  pull(average)

sample_vals <- t2_data %>% filter(well_id == 'sample ')

sample_vals_norm <- sample_vals %>% mutate(norm_value = (raw_val - t_avg) / (u_avg - t_avg),
                                           desc = factor(desc, levels = c('unpulsed', 'pulsed')))
sample_vals_norm$concentration <- as.numeric(sample_vals_norm$concentration)

models <- drm(
  norm_value ~ concentration, 
  desc, # grouping variable, group on pulse or no pulse
  data = sample_vals_norm, # contain concentration, measured output, and grouping variable
  fct = LL2.4() # type of fit, this a four parameter logistic curve
)

concentration_range <- seq(
  min(sample_vals_norm$concentration), 
  max(sample_vals_norm$concentration), 
  length = 100000
)


# Create a data frame for predictions
predictions <- expand.grid(
  concentration = concentration_range,
  desc = unique(sample_vals_norm$desc)
)

# Predict survival values for the fitted curves
pm <- predict(models, newdata = predictions, interval = "confidence")

predictions$p = pm[,1]
predictions$pmin = pm[,2]
predictions$pmax = pm[,3]


t2_cell_viability <- ggplot() +
  geom_point(data = sample_vals_norm,
             aes(x = concentration, y = norm_value, color = desc),
             size = 2) +
  geom_line(data = predictions,
            aes(x = concentration, y = p, color = desc, linetype = desc),
            size = 1) +
  theme_classic() +
  labs(x = 'PRAME BiTE Concentration [M]',
       y = 'Normalized T2 Cell Survival',
       title = '') + 
  scale_color_manual(values = c('#fde725', '#440154'),
                     labels = c("- Peptide Pulse", "+ Peptide Pulse")) +
  scale_linetype_manual(values = c("solid", "dashed"), guide = "none") +
  ylim(0, 1.75) +
  scale_x_log10() +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.position = "top",
        legend.margin = margin(-25, 0, 0, 0))

ggsave('work/results/bispecific/Figure5b_T2_norm_cell_viability_concentration.pdf', plot = t2_cell_viability, units = 'in', width = 4.5, height = 3.5,  dpi = 500)

 
max(sample_vals_norm$norm_value)


# Panel C Cytoxicity - Concentration Concentration ------------------------------------------

conc <- conc %>% mutate(Concentration = 10^seq(-8, -13, length.out = n()))

conc_long <- conc %>% pivot_longer(cols = -Concentration, values_to = 'normalized_cell_survival', names_to = 'cell_line')
conc_long <- conc_long %>% mutate(cell_line = gsub('\\_.*', "", conc_long$cell_line))


# Fit the 4PL model to the raw data (including replicates)
models <- drm(
  normalized_cell_survival ~ Concentration, 
  cell_line, # grouping variable, makes a model for each group provided here
  data = conc_long, # contain concentration, measured output, and grouping variable
  fct = LL.4() # type of fit, this a four parameter logistic curve
)

# Check the summary of the fit (to confirm parameters are from the 4PL model)
summary(models)

# Generate a sequence of concentrations to use for predictions to have smooth curves
concentration_range <- seq(
  min(conc_long$Concentration), 
  max(conc_long$Concentration), 
  length.out = 100000
)

# Create a data frame for predictions
predictions <- expand.grid(
  Concentration = concentration_range,
  cell_line = unique(conc_long$cell_line)
)

# Predict survival values for the fitted curves
pm <- predict(models, newdata = predictions, interval = 'confidence')
predictions$p = pm[,1]
predictions$pmin = pm[,2]
predictions$pmax = pm[,3]

# Aggregate data to calculate mean and standard deviation for each concentration and cell line
# Will use to plot the mean point as dot and then sd as lines
aggregated_data <- conc_long %>%
  group_by(cell_line, Concentration) %>%
  summarize(
    mean_survival = mean(normalized_cell_survival, na.rm = TRUE),
    sd_survival = sd(normalized_cell_survival, na.rm = TRUE),
    .groups = "drop")


NCI_cells_lines_plot <- ggplot() +
  # Add the fitted 4PL curve
  geom_line(data = predictions, 
            aes(x = Concentration, y = p, color = cell_line, linetype = cell_line), 
            size = 1) +
  # Add points for the mean
  geom_point(data = aggregated_data, 
             aes(x = Concentration, y = mean_survival, color = cell_line), 
             size = 3) +
  # Add error bars for standard deviation
  geom_errorbar(data = aggregated_data, 
                aes(x = Concentration, ymin = mean_survival - sd_survival, ymax = mean_survival + sd_survival, color = cell_line), 
                width = 0.2, size = 0.4) +
  # Log scale for x-axis
  scale_x_log10() +
  ylim(0,2) +
  # Labels and theme
  labs(
    x = "PRAME BiTE Concentration [M]", 
    y = "Normalized Cell Survival", 
    title = ""
  ) +
  theme_classic() +
  scale_color_manual(values = c('#35b779','#31688e')) +
  scale_linetype_manual(values = c("solid", "dashed"), guide = "none") +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.position = "top",
        legend.margin = margin(-25, 0, 0, 0))
  
ggsave('work/results/bispecific/Figure5c_NCI_celllines_norm_cell_viability_dose_titration.pdf', plot = NCI_cells_lines_plot,units = 'in', width = 4.5, height = 3.5, dpi = 500)

# Panel D - NUT carcinoma Cell viability after bispecific treatment  ---------------------------------------------------

# Calculate the mean and standard deviation of each "no bispecific" condition
control_summary <- bispecific_time_data %>%
  filter(bispecific_condition == "no_bispecific") %>%
  group_by(culture_condition, pulse_condition, time) %>%
  summarise(across(c(H460, `797`, `14169`, `1015`), 
                   list(mean = ~ mean(.x, na.rm = TRUE), 
                        sd = ~ sd(.x, na.rm = TRUE)), 
                   .names = "{.col}_{.fn}"),
            .groups = 'drop')

# Normalize the "bispecific" condition relative to the "no bispecific" condition
data_normalized <- bispecific_time_data %>%
  left_join(control_summary, 
            by = c("culture_condition", "pulse_condition", "time")) %>%
  mutate(
    H460_norm = (H460 / H460_mean),
    `797_norm` = (`797` / `797_mean`),
    `14169_norm` = (`14169` / `14169_mean`),
    `1015_norm` = (`1015` / `1015_mean`)
  )

# Convert to long format and calculate means/SDs for plotting
data_long <- data_normalized %>%
  filter(bispecific_condition == 'bispecific') %>% 
  dplyr::select(culture_condition, pulse_condition, time, 
         ends_with("_norm")) %>%
  pivot_longer(
    cols = ends_with("_norm"),
    names_to = "cell_line",
    values_to = "normalized_luminescence"
  ) %>%
  mutate(
    cell_line = gsub("_norm", "", cell_line), # Remove "_norm" from column names
    culture_condition = factor(culture_condition, 
                               levels = c("tumor_only", "E:T1", "E:T5"))
  )

# Aggregate to calculate mean and SD for each group
data_summary <- data_long %>%
  group_by(culture_condition, pulse_condition, time, cell_line) %>%
  summarise(
    mean_luminescence = mean(normalized_luminescence, na.rm = TRUE),
    sd_luminescence = sd(normalized_luminescence, na.rm = TRUE),
    .groups = "drop"
  )

# Create the plot
bispecific_timecourse_NC_cell_lines_plot <- ggplot(data_summary, 
       aes(x = culture_condition, 
           y = mean_luminescence, 
           color = pulse_condition)) +
  # Add points for the mean
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  # Add error bars for SD
  geom_errorbar(aes(ymin = mean_luminescence - sd_luminescence, 
                    ymax = mean_luminescence + sd_luminescence),
                width = 0.2, size = 0.8, 
                position = position_dodge(width = 0.5)) +
  # Facet grid by time and cell line
  facet_grid(time ~ cell_line, scales = "fixed", 
             labeller = labeller(
               cell_line = c(
                 "H460" = "H460",
                 "797" = "TC-797",
                 "14169" = "14169",
                 "1015" = "10-15"
               ),
               time = c(
                 "48" = "48 Hours",  
                 "96" = "96 Hours"
               )
             )) +
  # Logarithmic fold change and custom aesthetics
  theme_bw() +
  labs(
    title = "",
    y = "Luminescence Fold Change \n(BiTE/No BiTE)",
    x = "",
    color = ""
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    panel.spacing = unit(0.5, "lines")
  ) +
  scale_color_manual(
    values = c('#fde725', '#440154'),
    labels = c("- Peptide Pulse", "+ Peptide Pulse")
  ) +
  scale_x_discrete(
    labels = c(
      "tumor_only" = "Tumor Cells",
      "E:T1" = "1:1 E:T",
      "E:T5" = "5:1 E:T"
    )
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50")

ggsave('work/results/bispecific/Figure5d_PRAME_bispecific_across_cell_lines_normalized.pdf', plot = bispecific_timecourse_NC_cell_lines_plot, units = 'in', width = 6.5, height = 4.5, dpi = 500)
