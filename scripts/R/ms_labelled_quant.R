library(tidyverse)


ms <- read.csv('GitHub/PRAME_in_NUT_Carcinoma/data/shotgun_mass_spec/quant_labelled_MS_4_peptides.csv') %>% select(-c('X','X.1','X.2'))
ms$labelled <- factor(ms$labelled, levels = c('single', 'double'))
ms$peptide <- factor(ms$peptide, levels = c("YLHARLREL", "RLDQLLRHV", "SLLQHLIGL", "DVYENFRQW"))

ms$absolute_quantity[is.na(ms$absolute_quantity)] <- 0
ms$cell_line <- factor(ms$cell_line, levels = c('TC797','10-15', '14169','PER-403', 'JCM1', 'PDX'))

ms <- ms %>% mutate(alpha = ifelse(ms$absolute_quantity == 0.000, 0, 1))



ggplot(ms %>% filter(labelled == 'double'), aes(x = factor(cell_line), y = peptide)) +
  geom_point(aes(fill = absolute_quantity, alpha = as.factor(alpha)), 
             shape = 21, size = 7, stroke = 1,
             show.legend = c(alpha = FALSE)) +  
  scale_x_discrete(position = "top") +
  scale_fill_gradient(low = "blue", high = "red", 
                      name = "Absolute Quantity (fM)",
                      limits = c(min(ms$absolute_quantity, na.rm = TRUE), 
                                 max(ms$absolute_quantity, na.rm = TRUE))) + 
  theme_minimal() +
  theme(legend.position = "right", 
        panel.grid.major = element_blank(),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12, color='black'),
        axis.text.x = element_text(size = 12, color='black'),
        axis.text.y = element_text(size = 10, color='black')) +
  labs(x = "", y = "Peptide", fill = "Absolute Quantity")





lower <- 1e-6
tiny_neg <- -0.3  # how far below the baseline “undetected” bars should go

ms_plot <- ms %>%
  mutate(
    # Treat missing or 0 as undetected
    detected = ifelse(is.na(absolute_quantity) | absolute_quantity == 0, FALSE, TRUE),
    absolute_quantity = ifelse(is.na(absolute_quantity) | absolute_quantity == 0, lower, absolute_quantity),
    y_log  = log10(absolute_quantity),
    y_plot = ifelse(detected, y_log - log10(lower), tiny_neg)
  )

break_vals  <- c(1e-6, 1e-4, 1e-2, 1)
breaks_plot <- log10(break_vals) - log10(lower)

pd <- position_dodge2(width = 0.8, preserve = "single")

ggplot(ms_plot,
       aes(x = peptide,
           y = y_plot,
           fill = cell_line)) +
  geom_col(position = pd, width = 0.7) +
  scale_y_continuous(
    limits = c(tiny_neg, max(breaks_plot)),
    breaks = c(tiny_neg, breaks_plot),
    labels = c("Not detected", break_vals),
    expand = c(0, 0)
  ) +
  theme_classic() +
  facet_wrap(~labelled) +
  ylab("Absolute Quantity (fM)") + 
  coord_flip() + 
  theme(axis.text = element_text(color = 'black')) + 
  scale_fill_manual(values = c('TC797' = '#FF9933', 'PER-403' = '#6699FF', '10-15' = '#CC66CC',
                               '14169' = '#33FFCC', 'JCM1' = '#FFCC00','PDX' = '#CC6666' ),
                    labels = c('TC797' = 'TC-797', 'PER-403' = 'PER403'),
                    name = 'Cell Line')

