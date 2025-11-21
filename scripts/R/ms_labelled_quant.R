library(tidyverse)


ms <- read.csv('work/data/immunopeptidomics_mass_spec/quant_labelled_MS_4_peptides.csv')
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

