library(here)
library(tidyverse)


tma <- read.csv('~/work/data/tma/TMA_PRAME_C_French.csv')
tma$PRAME <- factor(tma$PRAME, levels = c('3', '2', '1', '0'))

prame_labels <- c("0" = "Negative", "1" = "1-25%", "2" = "26-75%", "3" = "≥ 75%")
prame_colors <- c("0" = "#0000CC", "1" = "#660099", "2" = "#990066", "3" = "#CC0000")

tma_count <- tma %>%
  group_by(PRAME) %>%
  summarise(n = n())


ggplot(tma, aes(x = "", fill = factor(PRAME))) + 
  geom_bar() + 
  scale_fill_manual(values = prame_colors,
                    labels = prame_labels,
                    name = "% PRAME\nPositive Cells") +
  ylab("Number of Samples") +
  xlab('') +
  theme_classic()
  
