library(ggplot2)

nutm1_heatmap_matrix <- read.delim("nutm1_matrix.tsv")
prame_heatmap_matrix <- read.delim("prame_matrix.tsv")

p <- ggplot(nutm1_heatmap_matrix, aes(x=pos, y=line, fill=peptide_status)) 
p <- p + geom_tile(height=0.95) 
p <- p + theme_classic() 
p <- p + ylab("Sample") 
p <- p + xlab("NUTM1 Amino Acid Position")

ggsave("NUTM1_peptide_heatmap.pdf", p, width=8, height=8, unit="in")

p2 <- ggplot(prame_heatmap_matrix, aes(x=pos, y=line, fill=peptide_status)) 
p2 <- p2 + geom_tile(height=0.95) 
p2 <- p2 + theme_classic() 
p2 <- p2 + ylab("Sample") 
p2 <- p2 + xlab("PRAME Amino Acid Position")

ggsave("PRAME_peptide_heatmap.pdf", p2, width=8, height=8, unit="in")
