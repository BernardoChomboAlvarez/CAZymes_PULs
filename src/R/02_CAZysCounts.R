## Libraries
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(textshape))
suppressPackageStartupMessages(library(stringr))

## Set working directory
setwd("/home/chombo/Documents/SEGOVIA/CAZymes_PULs/data/")

## Import MAGs list
MAGS.arr <- as.vector(read.table("1069MAGS_list.txt"))[[1]]

## Create a df with all the CAZy counts according to each MAG
CAZYS.df <- data.frame()
for (mag in MAGS.arr) {
    ## Filter the 5 MAGs that do not have significant CAZys associated with them
    if (mag == "3300007500_43") {
        next
    }
    
    ## Read the HMMER output for those matches with an evalue <= 1e-18 and a coverage >= 0.8
    df.temp <- read.table(file = paste(mag,"/output_dbcan4/hmmer.out",sep = ""),sep = "\t",skip = 1)
    colnames(df.temp) <- c("CAZy","length","geneID","genelength","evalue","start","end","genestart","geneend","coverage")
    df.temp <- df.temp %>%
        filter(evalue <= 1e-18,coverage >= 0.8) %>%
        select(CAZy,geneID,evalue,coverage)
    
    ## Extract the gene ID
    ids <- as.vector(df.temp$geneID)
    
    ## Read the overview.txt file to compare it with the hmmer.txt file
    df.temp <- read.table(paste(mag,"/output_dbcan4/overview.txt",sep = ""),header = FALSE,skip = 1)
    colnames(df.temp) <- c("geneID","EC","HMMER","dbCAN_sub","DIAMOND","Signalp","Tools")
    df.temp["MAG"] <- mag
    df.temp <- df.temp[is.element(df.temp$geneID,ids),]
    df.temp <- df.temp %>%
        filter(Tools >= 2) %>%
        select(MAG,HMMER,Signalp)
    df.temp$HMMER <- gsub("\\(.*?\\)", "", df.temp$HMMER)
    df.temp$Signalp <- gsub("\\(.*\\)", "",df.temp$Signalp)
    df.temp <- df.temp[df.temp$HMMER != "-",]
    
    ## Save the CAZYs frequency
    freq <- as.data.frame(table(df.temp$HMMER))
    freq <- freq[order(freq$Freq,decreasing = TRUE),]
    write.table(freq,file = paste(mag,"/output_dbcan4/",mag,".CAZyCounts.txt",sep = ""),sep = "\t",quote = FALSE,col.names = FALSE,row.names = FALSE)
    
    ## Clean the CAZys and do another df
    df.temp$HMMER <- gsub("\\+.*", "", df.temp$HMMER)
    
    ## Save the CAZYs frequency
    freq <- as.data.frame(table(df.temp$HMMER))
    freq <- freq[order(freq$Freq,decreasing = TRUE),]
    write.table(freq,file = paste(mag,"/output_dbcan4/",mag,".CAZyCounts.single.txt",sep = ""),sep = "\t",quote = FALSE,col.names = FALSE,row.names = FALSE)
    
    if (length(df.temp) == 0) {
        print(mag)
    }
    
    ## Store them in a single df
    CAZYS.df <- rbind(CAZYS.df,df.temp)
}
colnames(CAZYS.df) <- c("MAG","CAZys","Signalp")

## 61654 total CAZys
## 424 different CAZys
## 1064 final MAGs

## Get the summary
CAZysperMAG.df <- as.data.frame(table(CAZYS.df$MAG))
colnames(CAZysperMAG.df) <- c("MAGs","CAZys_Freq")
excludedMAGS <- data.frame(MAGs = c("3300007500_43","3300011404_8","3300012067_9","3300012151_17","3300012165_16"),
                           CAZys_Freq = c(0,0,0,0,0))
CAZysperMAG.df <- rbind(CAZysperMAG.df,excludedMAGS)

## Search for MAGs without CAZys
# setdiff(MAGS.arr,c(CAZysperMAG.df$MAGs))
# [1] "3300011404_8"  "3300012067_9"  "3300012151_17"
# [4] "3300012165_16"

## Save the data
write.table(CAZysperMAG.df,file = "CAZysperMAG.tsv",sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)

## Get the PULs diversity
cazys.diversity <- as.data.frame(table(CAZYS.df$CAZys))
write.table(cazys.diversity,file = "CAZYS.diversity.tsv",sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)


## PART 2
## Now we are going to separate the diversity according to the taxonomic ranks (class,order,family) and for the ecosystem distribution
## Read the MAGs metadata
MAGS.meta <- read.table(file = "1069_MAGs_metadata.tsv",sep = "\t",header = TRUE)
MAGS.meta <- MAGS.meta %>%
    select(genome_id,class,order,family,ecosystem_category)
MAGS.meta$ecosystem_category <- gsub(" ","_",MAGS.meta$ecosystem_category)

## Remove the excluded MAGs
MAGS.meta <- MAGS.meta[!is.element(MAGS.meta$genome_id,c("3300007500_43","3300011404_8","3300012067_9","3300012151_17","3300012165_16")),]

## PULs per Class
class.arr <- MAGS.meta %>%
    distinct(class)
class.arr <- class.arr[order(class.arr$class,decreasing = FALSE),]

full.div.df <- data.frame()
for (class in class.arr) {
    ## Extract all the MAGs annotated in each Class
    mags.arr.temp <- as.vector(MAGS.meta$genome_id[MAGS.meta$class == class])
    
    ## For each MAG, extract all the PUL counts
    CAZYsperMAG.temp <- CAZysperMAG.df[is.element(CAZysperMAG.df$MAGs,mags.arr.temp),]    
    
    ## Extract the complete diversity
    class.div.df <- data.frame(class = c(class),
                               frequency = sum(CAZYsperMAG.temp$CAZys_Freq))
    full.div.df <- rbind(full.div.df,class.div.df)
    
    ## Extract the PULs diversity
    rawCAZys.temp <- CAZYS.df[is.element(CAZYS.df$MAG,mags.arr.temp),]
    cazyscounts.temp <- as.data.frame(table(rawCAZys.temp$CAZys))
    colnames(cazyscounts.temp) <- c("CAZys","frequency")
    cazyscounts.temp <- cazyscounts.temp[order(cazyscounts.temp$frequency,decreasing = TRUE),]
    
    ## Save the data
    write.table(CAZYsperMAG.temp,file = paste("Taxa_Class/",class,"/",class,"_CAZysperMAG.tsv",sep = ""),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
    write.table(rawCAZys.temp,file = paste("Taxa_Class/",class,"/",class,"_CAZysDiversity.tsv",sep = ""),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
    write.table(cazyscounts.temp,file = paste("Taxa_Class/",class,"/",class,"_CAZysCounts.tsv",sep = ""),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
}
write.table(full.div.df,file = paste("Taxa_Class/","CAZY.TaxaClassDiversity.tsv",sep = ""),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)

## PULs per Order
order.arr <- MAGS.meta %>%
    distinct(order)
order.arr <- order.arr[order(order.arr$order,decreasing = FALSE),]

full.div.df <- data.frame()
for (ord in order.arr) {
    ## Extract all the MAGs annotated in each Class
    mags.arr.temp <- as.vector(MAGS.meta$genome_id[MAGS.meta$order == ord])
    
    ## For each MAG, extract all the PUL counts
    CAZYsperMAG.temp <- CAZysperMAG.df[is.element(CAZysperMAG.df$MAGs,mags.arr.temp),]    
    
    ## Extract the complete diversity
    order.div.df <- data.frame(order = c(ord),
                               frequency = sum(CAZYsperMAG.temp$CAZys_Freq))
    full.div.df <- rbind(full.div.df,order.div.df)
    
    ## Extract the PULs diversity
    rawCAZys.temp <- CAZYS.df[is.element(CAZYS.df$MAG,mags.arr.temp),]
    cazyscounts.temp <- as.data.frame(table(rawCAZys.temp$CAZys))
    colnames(cazyscounts.temp) <- c("CAZys","frequency")
    cazyscounts.temp <- cazyscounts.temp[order(cazyscounts.temp$frequency,decreasing = TRUE),]
    
    ## Save the data
    write.table(CAZYsperMAG.temp,file = paste("Taxa_Order/",ord,"/",ord,"_CAZysperMAG.tsv",sep = ""),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
    write.table(rawCAZys.temp,file = paste("Taxa_Order/",ord,"/",ord,"_CAZysDiversity.tsv",sep = ""),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
    write.table(cazyscounts.temp,file = paste("Taxa_Order/",ord,"/",ord,"_CAZysCounts.tsv",sep = ""),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
}
write.table(full.div.df,file = paste("Taxa_Order/","CAZY.TaxaOrderDiversity.tsv",sep = ""),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)

## PULs per Family
fam.arr <- MAGS.meta %>%
    distinct(family)
fam.arr <- fam.arr[order(fam.arr$family,decreasing = FALSE),]

full.div.df <- data.frame()
for (fam in fam.arr) {
    ## Extract all the MAGs annotated in each Class
    mags.arr.temp <- as.vector(MAGS.meta$genome_id[MAGS.meta$family == fam])
    
    ## For each MAG, extract all the PUL counts
    CAZYsperMAG.temp <- CAZysperMAG.df[is.element(CAZysperMAG.df$MAGs,mags.arr.temp),]    
    
    ## Extract the complete diversity
    family.div.df <- data.frame(order = c(fam),
                               frequency = sum(CAZYsperMAG.temp$CAZys_Freq))
    full.div.df <- rbind(full.div.df,family.div.df)
    
    ## Extract the PULs diversity
    rawCAZys.temp <- CAZYS.df[is.element(CAZYS.df$MAG,mags.arr.temp),]
    cazyscounts.temp <- as.data.frame(table(rawCAZys.temp$CAZys))
    colnames(cazyscounts.temp) <- c("CAZys","frequency")
    cazyscounts.temp <- cazyscounts.temp[order(cazyscounts.temp$frequency,decreasing = TRUE),]
    
    ## Save the data
    write.table(CAZYsperMAG.temp,file = paste("Taxa_Family/",fam,"/",fam,"_CAZysperMAG.tsv",sep = ""),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
    write.table(rawCAZys.temp,file = paste("Taxa_Family/",fam,"/",fam,"_CAZysDiversity.tsv",sep = ""),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
    write.table(cazyscounts.temp,file = paste("Taxa_Family/",fam,"/",fam,"_CAZysCounts.tsv",sep = ""),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
}
write.table(full.div.df,file = paste("Taxa_Family/","CAZY.TaxaFamilyDiversity.tsv",sep = ""),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)

## CAZys per Ecosystem
eco.arr <- MAGS.meta %>%
    distinct(ecosystem_category)
eco.arr <- eco.arr[order(eco.arr$ecosystem_category,decreasing = FALSE),]

full.div.df <- data.frame()
for (eco in eco.arr) {
    ## Extract all the MAGs annotated in each Class
    mags.arr.temp <- as.vector(MAGS.meta$genome_id[MAGS.meta$ecosystem_category == eco])
    
    ## For each MAG, extract all the PUL counts
    CAZYsperMAG.temp <- CAZysperMAG.df[is.element(CAZysperMAG.df$MAGs,mags.arr.temp),]    
    
    ## Extract the complete diversity
    eco.div.df <- data.frame(order = c(eco),
                                frequency = sum(CAZYsperMAG.temp$CAZys_Freq))
    full.div.df <- rbind(full.div.df,eco.div.df)
    
    ## Extract the PULs diversity
    rawCAZys.temp <- CAZYS.df[is.element(CAZYS.df$MAG,mags.arr.temp),]
    cazyscounts.temp <- as.data.frame(table(rawCAZys.temp$CAZys))
    colnames(cazyscounts.temp) <- c("CAZys","frequency")
    cazyscounts.temp <- cazyscounts.temp[order(cazyscounts.temp$frequency,decreasing = TRUE),]
    
    ## Save the data
    write.table(CAZYsperMAG.temp,file = paste("Ecosystem/",eco,"/",eco,"_CAZysperMAG.tsv",sep = ""),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
    write.table(rawCAZys.temp,file = paste("Ecosystem/",eco,"/",eco,"_CAZysDiversity.tsv",sep = ""),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
    write.table(cazyscounts.temp,file = paste("Ecosystem/",eco,"/",eco,"_CAZysCounts.tsv",sep = ""),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
}
write.table(full.div.df,file = paste("Ecosystem/","CAZY.EcosystemDiversity.tsv",sep = ""),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
