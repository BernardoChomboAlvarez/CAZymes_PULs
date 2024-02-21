## Libraries
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(textshape))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(vegan))
suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(circlize))

## FUNCTIONS
## EXTRACT CAZYMES
get.cazymes <- function(mags) {
    all.cazys <- data.frame()
    for (mag in mags) {
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
    return(all.cazys)
}

## EXTRACT PULS
get.puls <- function(mags) {
    all.puls <- data.frame()
    for (mag in mags) {
        # ## Filter the 5 MAGs that do not have significant PULs associated with them
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
    return(all.puls)
}

## DATA PREPARATION
## Extracting all the metadata for the subsecuent analysis
## There are no other metadata dataframes than the extracted here
## Set working directory
setwd("/home/chombo/Documents/SEGOVIA/CAZymes_PULs/data/")

## Read the MAGs metadata and extract the columns: class, order, family, ecosystem
MAGS.meta <- read.table(file = "1069_MAGs_metadata.tsv",sep = "\t",header = TRUE,row.names = 1)
MAGS.meta <- MAGS.meta %>%
    select(PULsFrequency,CAZysFrequency,class,order,family,specie,ecosystem_category)
colnames(MAGS.meta) <- c("PULS","CAZYS","Class","Order","Family","Specie","Ecosystem")
MAGS.meta$Ecosystem <- str_replace(MAGS.meta$Ecosystem," ","_")
MAGS.meta["mag"] <- rownames(MAGS.meta)

## Read the 343 MAGs new list
MAGS.arr <- as.vector(read.table("343_MAGslist.txt"))[[1]]

## Filter metadata
MAGS.meta <- MAGS.meta[is.element(MAGS.meta$mag,MAGS.arr),]

## Extract the counts
all.cazys <- get.cazymes(mags = MAGS.arr)

## Create the matrix
mtx <- all.cazys %>%
    spread(cazy, frequency) %>%
    column_to_rownames('mag') %>%
    as.matrix() %>%
    replace(is.na(.), 0)

#mtx <- vegdist(mtx,"bray")
dend <- as.dendrogram(hclust(vegdist(mtx,"bray")))

## Colors
dend <- dend %>%
    color_branches(k = 19) %>%
    raise.dendrogram (15) %>%
    hang.dendrogram()

## Control quality plot
ggplot(dend, labels = FALSE) + scale_y_reverse(expand = c(0.2, 0)) + coord_polar(theta="x")
png(filename = "CAZy_dend_test.png",width = 4000,height = 2000,res = 100)
factoextra::fviz_dend(dend, cex = 0.5, k = 20, 
          color_labels_by_k = TRUE, rect = TRUE)
dev.off()

## Write the tree
phy <- ape::as.phylo(hclust(vegdist(mtx,"bray")))
ape::write.tree(phy = phy,file='CAZy343_ward_Class.tree') 

## Repeat the whole process with PULs
all.puls <- get.puls(mags = MAGS.arr)

## Create the matrix
mtx <- all.puls %>%
    spread(pul, frequency) %>%
    column_to_rownames('mag') %>%
    as.matrix() %>%
    replace(is.na(.), 0)

#mtx <- vegdist(mtx,"bray")
dend <- as.dendrogram(hclust(vegdist(mtx,"bray")))

## Colors
dend <- dend %>%
    color_branches(k = 19) %>%
    raise.dendrogram (15) %>%
    hang.dendrogram()

## Control quality plot
ggplot(dend, labels = FALSE) + scale_y_reverse(expand = c(0.2, 0)) + coord_polar(theta="x")
png(filename = "PULs_dend_test.png",width = 4000,height = 2000,res = 100)
factoextra::fviz_dend(dend, cex = 0.5, k = 20, 
                      color_labels_by_k = TRUE, rect = TRUE)
dev.off()

## Write the tree
phy <- ape::as.phylo(hclust(vegdist(mtx,"bray")))
ape::write.tree(phy = phy,file='PULs343_ward_Class.tree') 
