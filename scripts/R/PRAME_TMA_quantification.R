library(here)
library(tidyverse)


tma <- read.csv('data/source_data/TMA_stain/TMA_PRAME_C_French.csv')
prame_labels <- c("0" = "< 1%", "1" = "1-25%", "2" = "26-75%", "3" = "≥ 75%")
prame_colors <- c("0" = "#31688EFF", "1" = "#CC99FF", "2" = "#FF9933", "3" = "#CC0000")

tma_count <- tma %>%
  group_by(PRAME) %>%
  summarise(n = n())


# Compute proportions and label positions
tma_count <- tma_count %>%
  arrange(desc(PRAME)) %>%
  mutate(prop = round(n / sum(n) * 100)) %>%
  mutate(ypos = cumsum(prop) - 0.5 * prop)


ggplot(tma_count, aes(x = "", y = prop, fill = factor(PRAME))) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = prame_colors, labels = prame_labels) +
  theme_void() +
  labs(fill = "% PRAME+ Cells") +
  geom_text(aes(y = ypos, label = paste0(round(prop), "%", "\n(N = ", n, ")")), size = 6)

