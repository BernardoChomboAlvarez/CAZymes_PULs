library(ggplot2)
library(tidyr)

setwd("C:/Users/chomb/Desktop")
metadata = read.csv("FULL_sortedMAGSHYPHEN.csv",header=FALSE, sep = ',')
newnames = c('MAGs','Ecosystem','Taxa','PUL_Counts','CAZy_Counts')
colnames(metadata) = newnames

png("CAZY_eco.boxplot.png", width=2500, height=3200, res = 500)
plot = ggplot(metadata,aes(x=factor(Ecosystem),y=CAZy_Counts)) +
  geom_boxplot(aes(fill=Ecosystem)) +
  labs(x='Ecosystem',y="CAZymes Counts") +
  scale_fill_manual(values = c("#191970","#00FFFF","#088F8F","#1434A4","#800020","#E97451","#C04000","#808000","#AAFF00","#008000","#32CD32","#0FFF50","#FFBF00","#0ECAB2","#960ECA","#CAC00E","#5F5E4B","#ADFF08","#FFAF00","#B78700","#044F0B","#E60F5E","#830B37","#E97451")) +
  geom_point(alpha = 0.1, position=position_jitter(height=0.01, width=0.01)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "none") +
  coord_flip() +
  theme(panel.grid.major = element_blank())

print(plot)

dev.off()