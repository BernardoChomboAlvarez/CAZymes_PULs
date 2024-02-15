## Libraries
library(stringr)
library(dplyr)

## First rearrange the file so the taxonomic ranges are splitted
# paste <(cut -f1,15 $FILE | tr ";" "\t" | awk '{ if (NR != 1) print }' | sed 's/[a-z]__//g' | awk -F "\t" '{ print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$8 }') <(awk -F "\t" '{if (NR != 1) print}' $FILE | cut --complement -d$'\t' -f15) 

## Read the file
setwd(dir = '/home/chombo/Documents/SEGOVIA/CAZymes_PULs/data/')
colnames.vec <- c("division","phylum","class","order","family","specie","genome_id","metagenome_id","genome_length","num_contigs","n50","num_16s","num_5s","num_23s","num_trna","completeness","contamination","quality_score","mimag_quality","otu_id","ecosystem_category","ecosystem_type","habitat","longitude","latitude")
MAGS.df <- read.table('MAGs_test.tsv',sep = "\t",col.names = colnames.vec)

## Filter the file
## Completness >= 99%
## Contamination <= 0.05%
MAGS.filtered <- subset(MAGS.df,completeness >= 99 & contamination <= 0.05 & family != "")
MAGS.filtered <- MAGS.filtered %>%
    select(genome_id,metagenome_id,everything())
MAGS.filtered <- MAGS.filtered[order(MAGS.filtered$genome_id),]

## Rename specific cases
MAGS.filtered$class[MAGS.filtered$class == "Bacilli_A"] <- "Bacilli"
MAGS.filtered$class[MAGS.filtered$class == "Thermoplasmata_A"] <- "Thermoplasmata"
MAGS.filtered$order[MAGS.filtered$order == "Bacillales_A"] <- "Bacillales"

## Save the filtered file
## Output: 1069 MAGs
write.table(MAGS.filtered,file = "1069_MAGs_metadata.tsv",sep = "\t",row.names = FALSE)

## Filter diversity by ecosystem type, class, order and family
## Ecosystem
eco.div <- as.data.frame(table(MAGS.filtered$ecosystem_category))
colnames(eco.div) <- c("ecosystem","frequency")

## Class
class.div <- as.data.frame(table(MAGS.filtered$class))
colnames(class.div) <- c("class","frequency")

## Order
order.div <- as.data.frame(table(MAGS.filtered$order))
colnames(order.div) <- c("order","frequency")

## Family
fam.div <- as.data.frame(table(MAGS.filtered$family))
colnames(fam.div) <- c("family","frequency")

# write.table(MAGS.filtered$genome_id,file = "1069MAGS_list.txt",sep = "",row.names = FALSE,quote = FALSE,col.names = FALSE)
# write.table(eco.div$ecosystem,file = "1069Eco_list.txt",sep = "",row.names = FALSE,quote = FALSE,col.names = FALSE)
# write.table(class.div$class,file = "1069Class_list.txt",sep = "",row.names = FALSE,quote = FALSE,col.names = FALSE)
# write.table(order.div$order,file = "1069Order_list.txt",sep = "",row.names = FALSE,quote = FALSE,col.names = FALSE)
# write.table(fam.div$family,file = "1069Fam_list.txt",sep = "",row.names = FALSE,quote = FALSE,col.names = FALSE)
