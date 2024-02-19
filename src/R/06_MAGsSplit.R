## Libraries
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(textshape))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(vegan))
suppressPackageStartupMessages(library(pracma))

## Set working directory
setwd("/home/chombo/Documents/SEGOVIA/CAZymes_PULs/data/")

## Load the MAGs metadata and extract the desired columns
# The PCoA done with the MAGs that had more PULs and CAZys diversity were the ones selected
# Here it is important to remark that the diversity is not equal to the amount of related items
# on the other hand, it means that we took in count how many DIFFERENT items were inside each MAG
## Read the MAGs metadata and extract the columns: class, order, family, ecosystem
MAGs.meta.df <- read.table(file = "1069_MAGs_metadata.tsv",sep = "\t",header = TRUE)
MAGs.meta.df <- MAGs.meta.df %>%
    select(genome_id,PULsFrequency,CAZysFrequency,ecosystem_category,class,order,family)
colnames(MAGs.meta.df) <- c("MAGs","PULS","CAZY","Ecosystem","Class","Order","Family")
MAGs.meta.df$Ecosystem <- str_replace(MAGs.meta.df$Ecosystem," ","_")

## Import the data for a matrix
MAGS.arr <- as.vector(read.table("1069MAGS_list.txt"))[[1]]


## PULS
## Extract the PULs that are present >= 25 MAGs
all.puls <- data.frame()
for (mag in MAGS.arr) {
    ## Filter the 5 MAGs that do not have significant PULs associated with them
    if (mag == "3300005326_35" || mag == "3300005370_1" || mag == "3300017650_9" || mag == "3300027284_65" || mag == "3300027607_50") {
        next
    }
    
    df.temp <- read.table(file = paste(mag,"/",mag,".matches_counts.out",sep = ""),sep = "\t",)
    colnames(df.temp) <- c("pul","frequency")
    df.temp["mag"] <- mag
    df.temp <- df.temp %>%
        select(mag,pul,everything())
    
    ## Create a dataframe with the nodes: mag = from, pul = to, freq = value
    all.puls <- rbind(all.puls,df.temp)
}


## List of MAGs that have >= puls than their mean
mags.data <- as.data.frame(table(all.puls$mag))
# summary(mags.data$Freq)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.00   14.00   25.00   29.65   39.25  173.00
mags.data.cut <- mags.data[mags.data$Freq >= 30,]
# length(mags.data.cut$Var1)
# [1] 445
puls.filtered.mags.vec <- as.vector(mags.data.cut$Var1)

## List of PULs that are present >= their mean
puls.data <- as.data.frame(table(all.puls$pul))
# summary(puls.data$Freq)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00    2.00    6.00   23.15   18.00  695.00
puls.data.cut <- puls.data[puls.data$Freq >= 23,]
# length(puls.data.cut$Var1)
# [1] 303
filtered.puls.vec <- as.vector(puls.data.cut$Var1)


## CAZYMES
## Extract the CAZys
all.cazys <- data.frame()
for (mag in MAGS.arr) {
    ## Filter the 5 MAGs that do not have significant PULs associated with them
    if (mag == "3300007500_43" || mag == "3300011404_8" || mag == "3300012067_9" || mag == "3300012151_17" || mag == "3300012165_16" || mag == "3300007500_43") {
        next
    }
    
    df.temp <- read.table(file = paste(mag,"/output_dbcan4/",mag,".CAZyCounts.single.txt",sep = ""),sep = "\t",)
    colnames(df.temp) <- c("cazy","frequency")
    df.temp["mag"] <- mag
    df.temp <- df.temp %>%
        select(mag,cazy,everything())
    
    ## Create a dataframe with the nodes: mag = from, cazy = to, freq = value
    all.cazys <- rbind(all.cazys,df.temp)
}


## List of MAGs that have >= cazys than their mean
mags.data <- as.data.frame(table(all.cazys$mag))
# summary(mags.data$Freq)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00   17.00   26.00   32.28   43.00  108.00 
mags.data.cut <- mags.data[mags.data$Freq >= 32,]
# length(mags.data.cut$Var1)
# [1] 453
cazys.filtered.mags.vec <- as.vector(mags.data.cut$Var1)

## List of CAZys that are present >= their median
cazys.data <- as.data.frame(table(all.cazys$cazy))
# summary(cazys.data$Freq)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.0     6.0    24.0    81.0    89.5  1041.0
cazys.data.cut <- cazys.data[cazys.data$Freq >= 24,]
# length(cazys.data.cut$Var1)
# [1] 214
filtered.cazys.vec <- as.vector(cazys.data.cut$Var1)


## Filter the MAGs metadata according to the merge of the two filtered mags arrays
merge.arr <- c(cazys.filtered.mags.vec,puls.filtered.mags.vec)
merge.arr <- unique(merge.arr)
MAGs.meta.cleaned.df <- MAGs.meta.df[is.element(MAGs.meta.df$MAGs,merge.arr),]

## Write the output
write.table(MAGs.meta.cleaned.df,file = "747_MAGS_metadata.tsv",sep = "\t",quote = FALSE,col.names = TRUE,row.names = FALSE)
