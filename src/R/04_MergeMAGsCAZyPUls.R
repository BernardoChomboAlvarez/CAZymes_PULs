## Libraries
suppressPackageStartupMessages(library(textshape))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))

## Set working directory
setwd('/home/chombo/Documents/SEGOVIA/CAZymes_PULs/data/')

## Extract the PULs counts per MAG
PULsperMAG.df <- read.table(file = "PULsperMAG.tsv",sep = "\t",header = FALSE)
colnames(PULsperMAG.df) <- c("genome_id","PULsFrequency")

## Extract the CAZys counts per MAG
CAZysperMAG.df <- read.table(file = "CAZysperMAG.tsv",sep = "\t",header = FALSE)
colnames(CAZysperMAG.df) <- c("genome_id","CAZysFrequency")

## Load MAGs metadata
MAGs.meta.df <- read.table(file = "1069_MAGs_metadata.tsv",sep = "\t",header = TRUE)

## Merge the three data frames
MAGs.meta.df <- merge(MAGs.meta.df, PULsperMAG.df,
                      by = 'genome_id', all = TRUE)
MAGs.meta.df <- merge(MAGs.meta.df, CAZysperMAG.df,
                      by = 'genome_id', all = TRUE)

## Rewrite MAGs metadata
write.table(MAGs.meta.df,file = "1069_MAGs_metadata.tsv",sep = "\t",quote = FALSE,row.names = FALSE,col.names = TRUE)

