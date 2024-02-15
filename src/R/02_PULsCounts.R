## Libraries
library(dplyr)

## Set working directory
setwd("/home/chombo/Documents/SEGOVIA/CAZymes_PULs/data/")

## Import MAGs list
MAGS.arr <- as.vector(read.table("1069MAGS_list.txt"))[[1]]

## Create a df with all the PUL counts according to each MAG
PULS.df <- data.frame()
for (mag in MAGS.arr) {
    ## Filter the 5 MAGs that do not have significant PULs associated with them
    if (mag == "3300005326_35" || mag == "3300005370_1" || mag == "3300017650_9" || mag == "3300027284_65" || mag == "3300027607_50") {
        next
    }
    df.temp <- read.table(paste(mag,"/",mag,".matches_evalue.out",sep = ""))
    colnames(df.temp) <- c("PULs","evalue")
    df.temp["MAG"] <- mag
    
    ## Save the PULs frequency
    freq <- as.data.frame(table(df.temp$PULs))
    freq <- freq[order(freq$Freq,decreasing = TRUE),]
    write.table(freq,file = paste(mag,"/",mag,".matches_counts.out",sep = ""),sep = "\t",quote = FALSE,col.names = FALSE,row.names = FALSE)
    
    ## Store them in a single df
    PULS.df <- rbind(PULS.df,df.temp)
}

## 36172 total PULs
## 1064 final MAGs

## Get the summary
PULSperMAG.df <- as.data.frame(table(PULS.df$MAG))
colnames(PULSperMAG.df) <- c("MAGs","PULs_Freq")
excludedMAGS <- data.frame(MAGs = c("3300005326_35","3300005370_1","3300017650_9","3300027284_65","3300027607_50"),
                           PULs_Freq = c(0,0,0,0,0))
PULSperMAG.df <- rbind(PULSperMAG.df,excludedMAGS)

## Save the data
write.table(PULSperMAG.df,file = "PULsperMAG.tsv",sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)


## PART 2
## Now we are going to separate the diversity according to the taxonomic ranks (class,order,family) and for the ecosystem distribution
## Read the MAGs metadata
MAGS.meta <- read.table(file = "1069_MAGs_metadata.tsv",sep = "\t",header = TRUE)
MAGS.meta <- MAGS.meta %>%
    select(genome_id,class,order,family,ecosystem_category)
MAGS.meta$ecosystem_category <- gsub(" ","_",MAGS.meta$ecosystem_category)

## Remove the excluded MAGs
MAGS.meta <- MAGS.meta[!is.element(MAGS.meta$genome_id,c("3300005326_35","3300005370_1","3300017650_9","3300027284_65","3300027607_50")),]

## PULs per Class
class.arr <- MAGS.meta %>%
    distinct(class)
class.arr <- class.arr[order(class.arr$class,decreasing = FALSE),]

full.div.df <- data.frame()
for (class in class.arr) {
    ## Extract all the MAGs annotated in each Class
    mags.arr.temp <- as.vector(MAGS.meta$genome_id[MAGS.meta$class == class])
    
    ## For each MAG, extract all the PUL counts
    PULSperMAG.temp <- PULSperMAG.df[is.element(PULSperMAG.df$MAGs,mags.arr.temp),]    
    
    ## Extract the complete diversity
    class.div.df <- data.frame(class = c(class),
                               frequency = sum(PULSperMAG.temp$PULs_Freq))
    full.div.df <- rbind(full.div.df,class.div.df)
    
    ## Extract the PULs diversity
    rawPULs.temp <- PULS.df[is.element(PULS.df$MAG,mags.arr.temp),]
    pulscounts.temp <- as.data.frame(table(rawPULs.temp$PULs))
    colnames(pulscounts.temp) <- c("PULs","frequency")
    pulscounts.temp <- pulscounts.temp[order(pulscounts.temp$frequency,decreasing = TRUE),]
    
    ## Save the data
    write.table(PULSperMAG.temp,file = paste("Taxa_Class/",class,"/",class,"_PULsperMAG.tsv",sep = ""),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
    write.table(rawPULs.temp,file = paste("Taxa_Class/",class,"/",class,"_PULsDiversity.tsv",sep = ""),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
    write.table(pulscounts.temp,file = paste("Taxa_Class/",class,"/",class,"_PULsCounts.tsv",sep = ""),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
}
write.table(full.div.df,file = paste("Taxa_Class/","TaxaClassDiversity.tsv",sep = ""),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)

## PULs per Order
order.arr <- MAGS.meta %>%
    distinct(order)
order.arr <- order.arr[order(order.arr$order,decreasing = FALSE),]

full.div.df <- data.frame()
for (ord in order.arr) {
    ## Extract all the MAGs annotated in each Order
    mags.arr.temp <- as.vector(MAGS.meta$genome_id[MAGS.meta$order == ord])
    
    ## For each MAG, extract all the PUL counts
    PULSperMAG.temp <- PULSperMAG.df[is.element(PULSperMAG.df$MAGs,mags.arr.temp),]    

    ## Extract the complete diversity
    order.div.df <- data.frame(order = c(ord),
                               frequency = sum(PULSperMAG.temp$PULs_Freq))
    full.div.df <- rbind(full.div.df,order.div.df)
    
    ## Extract the PULs diversity
    rawPULs.temp <- PULS.df[is.element(PULS.df$MAG,mags.arr.temp),]
    pulscounts.temp <- as.data.frame(table(rawPULs.temp$PULs))
    colnames(pulscounts.temp) <- c("PULs","frequency")
    pulscounts.temp <- pulscounts.temp[order(pulscounts.temp$frequency,decreasing = TRUE),]
    
    ## Save the data
    write.table(PULSperMAG.temp,file = paste("Taxa_Order/",ord,"/",ord,"_PULsperMAG.tsv",sep = ""),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
    write.table(rawPULs.temp,file = paste("Taxa_Order/",ord,"/",ord,"_PULsDiversity.tsv",sep = ""),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
    write.table(pulscounts.temp,file = paste("Taxa_Order/",ord,"/",ord,"_PULsCounts.tsv",sep = ""),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
}
write.table(full.div.df,file = paste("Taxa_Order/","TaxaOrderDiversity.tsv",sep = ""),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)

## PULs per Family
fam.arr <- MAGS.meta %>%
    distinct(family)
fam.arr <- fam.arr[order(fam.arr$family,decreasing = FALSE),]

full.div.df <- data.frame()
for (fam in fam.arr) {
    ## Extract all the MAGs annotated in each Family
    mags.arr.temp <- as.vector(MAGS.meta$genome_id[MAGS.meta$family == fam])
    
    ## For each MAG, extract all the PUL counts
    PULSperMAG.temp <- PULSperMAG.df[is.element(PULSperMAG.df$MAGs,mags.arr.temp),]    
    
    ## Extract the complete diversity
    fam.div.df <- data.frame(family = c(fam),
                               frequency = sum(PULSperMAG.temp$PULs_Freq))
    full.div.df <- rbind(full.div.df,fam.div.df)
    
    ## Extract the PULs diversity
    rawPULs.temp <- PULS.df[is.element(PULS.df$MAG,mags.arr.temp),]
    pulscounts.temp <- as.data.frame(table(rawPULs.temp$PULs))
    colnames(pulscounts.temp) <- c("PULs","frequency")
    pulscounts.temp <- pulscounts.temp[order(pulscounts.temp$frequency,decreasing = TRUE),]
    
    ## Save the data
    write.table(PULSperMAG.temp,file = paste("Taxa_Family/",fam,"/",fam,"_PULsperMAG.tsv",sep = ""),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
    write.table(rawPULs.temp,file = paste("Taxa_Family/",fam,"/",fam,"_PULsDiversity.tsv",sep = ""),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
    write.table(pulscounts.temp,file = paste("Taxa_Family/",fam,"/",fam,"_PULsCounts.tsv",sep = ""),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
}
write.table(full.div.df,file = paste("Taxa_Family/","TaxaFamilyDiversity.tsv",sep = ""),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)

## PULs per Ecosystem
eco.arr <- MAGS.meta %>%
    distinct(ecosystem_category)
eco.arr <- eco.arr[order(eco.arr$ecosystem_category,decreasing = FALSE),]

full.div.df <- data.frame()
for (eco in eco.arr) {
    ## Extract all the MAGs annotated in each Ecosystem
    mags.arr.temp <- as.vector(MAGS.meta$genome_id[MAGS.meta$ecosystem_category == eco])
    
    ## For each MAG, extract all the PUL counts
    PULSperMAG.temp <- PULSperMAG.df[is.element(PULSperMAG.df$MAGs,mags.arr.temp),]    
    
    ## Extract the complete diversity
    eco.div.df <- data.frame(ecosystem = c(eco),
                             frequency = sum(PULSperMAG.temp$PULs_Freq))
    full.div.df <- rbind(full.div.df,eco.div.df)
    
    ## Extract the PULs diversity
    rawPULs.temp <- PULS.df[is.element(PULS.df$MAG,mags.arr.temp),]
    pulscounts.temp <- as.data.frame(table(rawPULs.temp$PULs))
    colnames(pulscounts.temp) <- c("PULs","frequency")
    pulscounts.temp <- pulscounts.temp[order(pulscounts.temp$frequency,decreasing = TRUE),]
    
    ## Save the data
    write.table(PULSperMAG.temp,file = paste("Ecosystem/",eco,"/",eco,"_PULsperMAG.tsv",sep = ""),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
    write.table(rawPULs.temp,file = paste("Ecosystem/",eco,"/",eco,"_PULsDiversity.tsv",sep = ""),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
    write.table(pulscounts.temp,file = paste("Ecosystem/",eco,"/",eco,"_PULsCounts.tsv",sep = ""),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
}
write.table(full.div.df,file = paste("Ecosystem/","EcosystemDiversity.tsv",sep = ""),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
