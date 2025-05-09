library(ggplot2)
library(gridExtra)
library(NatParksPalettes)
library(svglite)

cols <- c("TC-797"="#0067A2", "10-15"="#DFCB91", "14169"="#CB7223", "JCM1"="#289A84", "PDX"="#7FA4C2", "PER403" = "#AF7E56")

nutm1_heatmap_matrix <- read.delim("nutm1_matrix.tsv")
prame_heatmap_matrix <- read.delim("prame_matrix.tsv")

nutm1p <- ggplot(nutm1_heatmap_matrix, aes(x=pos, y=line, alpha=as.numeric(peptide_status), fill=line)) 
nutm1p <- nutm1p + geom_tile(height=1, width=1) 
nutm1p <- nutm1p + theme_classic() 
nutm1p <- nutm1p + ylab("NUT") 
nutm1p <- nutm1p + xlab("Amino Acid Position") 
nutm1p <- nutm1p + scale_fill_discrete(guide = guide_legend()) 
nutm1p <- nutm1p + theme(legend.position="bottom") 
nutm1p <- nutm1p + guides(alpha = "none") 
nutm1p <- nutm1p + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) 
nutm1p <- nutm1p + scale_x_continuous(limits = c(1,1132), expand = c(0, 0)) 
nutm1p <- nutm1p + guides(fill="none")
nutm1p <- nutm1p + scale_fill_manual(values=cols)

pramep <- ggplot(prame_heatmap_matrix, aes(x=pos, y=line, alpha=as.numeric(peptide_status), fill=line))
pramep <- pramep + geom_tile(height=1, width=1)
pramep <- pramep + theme_classic()
pramep <- pramep + ylab("PRAME")
pramep <- pramep + scale_fill_discrete(guide = guide_legend()) 
pramep <- pramep + guides(alpha = "none")
pramep <- pramep + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
pramep <- pramep + scale_x_continuous(limits = c(1,509), expand = c(0, 0))
pramep <- pramep + xlab("")
pramep <- pramep + guides(fill=guide_legend(ncol=2, title="Cell Line"))
pramep <- pramep + theme(legend.title.align=0.5)
pramep <- pramep + scale_fill_manual(values=cols)

merged1 <- grid.arrange(pramep, nutm1p, ncol=1, nrow=2, heights=c(1,1))

ggsave("peptide_heatmap.pdf", merged1, width=8, height=2, unit="in")

SLL <- subset(prame_heatmap_matrix, pos > 421 & pos < 429)

sllp <- ggplot(SLL, aes(x=pos, y=line, alpha=as.numeric(peptide_status), fill=line))
sllp <- sllp + theme_classic()
sllp <- sllp + geom_tile()
sllp <- sllp + guides(fill="none")
sllp <- sllp + guides(alpha="none")
sllp <- sllp + scale_x_continuous(limits = c(421,429), expand = c(0, 0))
sllp <- sllp + xlab("")
sllp <- sllp + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
sllp <- sllp + ylab("")
sllp <- sllp + scale_x_continuous(breaks = 421:429)
sllp <- sllp + theme(axis.line = element_blank())
sllp <- sllp + geom_tile(color="black")
sllp <- sllp + annotate("rect", xmin=421.5, xmax=424.5, ymin=0.5, ymax=6.5, fill="grey", alpha=.75)
sllp <- sllp + annotate("rect", xmin=425.5, xmax=428.5, ymin=0.5, ymax=6.5, fill="grey", alpha=.75)
sllp <- sllp + annotate("rect", xmin=424.5, xmax=425.5, ymin=2.5, ymax=3.5, fill="white", color="black", size=.25)
sllp <- sllp + annotate("rect", xmin=424.5, xmax=425.5, ymin=3.5, ymax=4.5, fill="#289A84", color="black", size=.25)
sllp <- sllp + annotate("text", x=425, y=3.5, label="*")
sllp <- sllp + xlab("SLLQHLIGL")
sllp <- sllp + annotate("rect", xmin=424.5, xmax=425.5, ymin=0, ymax=7, fill="transparent", color="black", size=2)
sllp <- sllp + scale_fill_manual(values=cols)

ggsave("ssl_heatmap.pdf", sllp, width=8, height=2, unit="in")

merged2 <- grid.arrange(sllp, pramep, nutm1p, ncol=1, nrow=3, heights=c(1,1,1))

ggsave("all_heatmap.pdf", merged2, width=8, height=3, unit="in")
ggsave("all_heatmap.svg", merged2, width=8, height=3, unit="in")
