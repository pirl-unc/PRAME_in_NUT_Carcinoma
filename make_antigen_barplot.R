library(ggplot2)
library("NatParksPalettes")
library(svglite)

antigen_counts <- read.delim("antigen_counts.tsv")

p <- ggplot(antigen_counts, aes(x=reorder(Line, -pMHC_Count), y=pMHC_Count, fill=Antigen_Source)) 
p <- p + geom_bar(stat='identity')
p <- p + theme_classic() 
p <- p + ylab("pMHC Count")
p <- p + xlab("Cell Line")
p <- p + theme(plot.title = element_text(hjust = 0.5))
p <- p + labs(fill="Antigen Source")
p <- p + theme(legend.title = element_text(hjust = 0.5))
p <- p + scale_fill_manual(values=natparks.pals("SmokyMtns"))

ggsave("antigen_source_barplot.pdf", p, width=8, height=8, unit="in")
ggsave("antigen_source_barplot.svg", p, width=8, height=8, unit="in")
