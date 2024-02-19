## Libraries
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(textshape))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(vegan))
suppressPackageStartupMessages(library(pracma))


## FUNCTIONS
plot.pcoa <- function(mtx,metadata,labels,component) {
    if (length(component) >= 6) {
        short.name <- substr(component, 1, 3) 
    } else {
        short.name <- component
    }
    
    PULs <- mtx
    Bray=vegdist(PULs,"bray")
    PULs_bray=capscale(Bray~-1)
    PULs_bray_eig = eigenvals(PULs_bray)
    percentage_variance_explained <- PULs_bray_eig / sum(PULs_bray_eig)
    sum_percentage_variance_explained <- cumsum(PULs_bray_eig / sum(PULs_bray_eig))
    
    xlabel= as.numeric(format(round((percentage_variance_explained[1]*100), 2), nsmall = 2))
    xlabel= sprintf("%.2f %%", xlabel)
    xlabel= paste ("PCo1 (", xlabel, ")")
    ylabel= as.numeric(format(round((percentage_variance_explained[2]*100), 2), nsmall = 2))
    ylabel= sprintf("%.2f %%", ylabel)
    ylabel= paste ("PCo2 (", ylabel, ")")
    
    gsites<- c(metadata[component])
    Permanova_test <- adonis(formula = PULs ~ gsites[[1]])
    pval <- Permanova_test$aov.tab$`Pr(>F)`[1]
    H_CLustering=hclust(vegdist(PULs,"bray"))
    
    png(paste("bray_",short.name,"_PULs_PCoA.png",sep = ""), width=3200, height=3200, res = 900)
    plot(PULs_bray,
         family="Times",
         type="n",
         xlab="",
         ylab="",
         ylim=c(-3.0, 2.5),
         xlim=c(-2.2, 5.5),
         cex.axis=0.4,
         tck = -0.01,
         mgp = c(3, 0.2, 0),
         xaxp  = c(-4, 4, 8),
         panel.first=grid(col = "white",lty=0))
    title(ylab=ylabel,
          line=1.2,
          cex.lab=0.4)
    title(xlab=xlabel,
          line=1.2, cex.lab=0.4)
    title(main = paste("PCoA PULs, ",component,sep = ""),
          cex.main = 0.7)
    abline(h=0,
           v=0,
           col = "white",
           lty = 1,
           lwd = 1.5)
    abline(h=-10:10,
           v=-10:10,
           col = "lightgray",
           lty = "dotted",
           lwd = 0.5)
    par(lty=2)
    
    # Adding the confidence intervals at two different levels (97.5% and 95%)
    ordiellipse(PULs_bray, groups = gsites[[1]],  kind = "sd",conf= 0.95, draw ="polygon", lwd = 0.5)
    par(lty=2)
    # ordiellipse(PULs_bray, groups = gsites,  kind = "sd",conf= 0.95, draw ="polygon", col = "grey90", lwd = 0.5, lty=2)
    # par(lty=2)
    
    ordicluster(PULs_bray, H_CLustering, prune= 4, display = "sites", col = "grey", lwd = 0.5)
    points(PULs_bray, cex= 0.2, pch=21, col="black", bg= metadata$Colvec, lwd = 0.3)
    
    lgend.x=4.5
    legend.y=2.0
    
    if (length(labels) <= 25) {
        legend(lgend.x-0.5,legend.y, pt.cex = 0.25 , pt.lwd = 0.2, labels,bty = "n" ,pch = 21,col="black",pt.bg = colvec, cex = 0.25)
        legend((lgend.x-0.9),(legend.y+0.3), cex = 0.3 , component,bty = "n" ,pch = 21,col="white",pt.bg = c("white"))
        text(x= 3.0 , y= -2.35 ,expression(paste("Permanova")), cex = 0.3)
        text(x= 3.0 , y= -2.55 ,paste("p-value = ", pval, sep = ""), cex = 0.3)
    } else if (length(labels) <= 80) {
        legend(lgend.x-2.0,(legend.y+0.3), pt.cex = 0.25 , pt.lwd = 0.2, labels,bty = "n" ,pch = 21,col="black",pt.bg = colvec, cex = 0.20,ncol=2)
        legend((lgend.x-2.3),(legend.y+0.6), cex = 0.3 , component,bty = "n" ,pch = 21,col="white",pt.bg = c("white"))
        text(x= 0.5 , y= -2.35 ,expression(paste("Permanova")), cex = 0.3)
        text(x= 0.5 , y= -2.55 ,paste("p-value = ", pval, sep = ""), cex = 0.3)
    } else {
        legend(lgend.x-2.0,(legend.y+0.3), pt.cex = 0.2 , pt.lwd = 0.15, labels,bty = "n" ,pch = 21,col="black",pt.bg = colvec, cex = 0.15,ncol=2)
        legend((lgend.x-2.3),(legend.y+0.6), cex = 0.25 , component,bty = "n" ,pch = 21,col="white",pt.bg = c("white"))
        text(x= 0.5 , y= -2.35 ,expression(paste("Permanova")), cex = 0.3)
        text(x= 0.5 , y= -2.55 ,paste("p-value = ", pval, sep = ""), cex = 0.3)
    }
    dev.off()
}

## DATA PREPARATION
## Extracting all the metadata for the subsecuent analysis
## There are no other metadata dataframes than the extracted here
## Set working directory
setwd("/home/chombo/Documents/SEGOVIA/CAZymes_PULs/data/")

## Read the MAGs metadata and extract the columns: class, order, family, ecosystem
MAGS.meta <- read.table(file = "1069_MAGs_metadata.tsv",sep = "\t",header = TRUE,row.names = 1)
MAGS.meta <- MAGS.meta %>%
    select(PULsFrequency,CAZysFrequency,class,order,family,ecosystem_category)
colnames(MAGS.meta) <- c("PULS","CAZYS","Class","Order","Family","Ecosystem")
MAGS.meta$Ecosystem <- str_replace(MAGS.meta$Ecosystem," ","_")
MAGS.meta["Row.names"] <- rownames(MAGS.meta)

## Filter by Class
classes.arr <- c("Actinobacteria","Bacilli","Bacteroidia","Clostridia","Deinococci","Gammaproteobacteria","Negativicutes","Spirochaetia","Synergistia")
MAGS.meta <- MAGS.meta[is.element(MAGS.meta$Class,classes.arr),]

## Filter by Order
orders.arr <- c("Acidaminococcales","Bacteroidales","Betaproteobacteriales","Chitinophagales","Deinococcales","Enterobacterales","Erysipelotrichales","Flavobacteriales","Lachnospirales","Lactobacillales","Oscillospirales","Pseudomonadales","Veillonellales")
MAGS.meta <- MAGS.meta[is.element(MAGS.meta$Order,orders.arr),]

## Extract MAGs
MAGS.arr <- rownames(MAGS.meta)


all.puls <- data.frame()
for (mag in MAGS.arr) {
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


## List of MAGs that have >= puls than their mean
mags.data <- as.data.frame(table(all.puls$mag))
# summary(mags.data$Freq)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.00   14.00   25.00   29.65   39.25  173.00
mags.data.cut <- mags.data[mags.data$Freq >= 25,]
# length(mags.data.cut$Var1)
# [1] 453
filtered.mags.vec <- as.vector(mags.data.cut$Var1)

## List of PULs that are present >= their mean
puls.data <- as.data.frame(table(all.puls$pul))
# summary(puls.data$Freq)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00    2.00    6.00   23.15   18.00  695.00
puls.data.cut <- puls.data[puls.data$Freq >= 21,]
# length(puls.data.cut$Var1)
# [1] 278
filtered.puls.vec <- as.vector(puls.data.cut$Var1)

# ===================================================================
## At the end we will be working with a matrix 445 MAGs x 303 PULs
# ===================================================================


## UNIVERSAL MATRIX GENERATION
## Use the current data
## Filter the original data by the filtered vectors
matrix.data <- all.puls[is.element(all.puls$mag,filtered.mags.vec),]
matrix.data <- matrix.data[is.element(matrix.data$pul,filtered.puls.vec),]

## Generate the matrix
mtx <- matrix.data %>%
    spread(pul, frequency) %>%
    column_to_rownames('mag') %>%
    as.matrix() %>%
    replace(is.na(.), 0)


## ECOSYSTEM DATASET
## Use the MAGs metadata and order that dataframe by 'Ecosystem' column
## Assign the MAGS.meta content to another variable
ECO.meta <- MAGS.meta[is.element(MAGS.meta$Row.names,filtered.mags.vec),]
ECO.meta <- ECO.meta[order(ECO.meta$Ecosystem,decreasing = FALSE),]
unique.df <- as.data.frame(unique(ECO.meta$Ecosystem))
colnames(unique.df) <- c("Ecosystem")
colvec<- c("#191970","#00FFFF","#088F8F","#1434A4","#CD7F32","#800020","#E97451","#6E260E","#8B0000","#C04000","#808000","#AAFF00","#228B22","#7CFC00")#,"#008000","#32CD32","#0FFF50")#,"#FFBF00") #,"#9370DB","#8A2BE2","#9932CC","#DA70D6","#FF69B4"
unique.df["Colvec"] <- colvec

## Merge both data frames
ECO.meta <- merge(ECO.meta,unique.df,
                  by = "Ecosystem", all = TRUE)
rownames(ECO.meta) <- ECO.meta$Row.names

## Order the matrix
mags.sorted <- as.vector(ECO.meta$Row.names)
PULs <- mtx[match(mags.sorted,as.character(rownames(mtx))),]
PULs <- as.data.frame(PULs)
PULs <- tibble::rownames_to_column(PULs, "MAG")
PULs$MAG <- gsub("^X","",PULs$MAG)
PULs2 <- PULs[,-1]
rownames(PULs2) <- PULs[,1]
PULs <- PULs2
PULs <- PULs %>%
    replace(is.na(.), 0)

## Define the labels
eco.labels <- ECO.meta %>%
    distinct(Ecosystem) %>%
    unname()
eco.labels <- c(eco.labels)[[1]]
eco.labels <- sort(eco.labels)

## Plot the PCoA
plot.pcoa(mtx = PULs,metadata = ECO.meta,labels = eco.labels,component = "Ecosystem")


## CLASS DATASET
## Use the MAGs metadata and order that dataframe by 'Class' column
## Assign the MAGS.meta content to another variable
CLASS.meta <- MAGS.meta[is.element(MAGS.meta$Row.names,filtered.mags.vec),]
CLASS.meta <- CLASS.meta[order(CLASS.meta$Class,decreasing = FALSE),]
unique.df <- as.data.frame(unique(CLASS.meta$Class))
colnames(unique.df) <- c("Class")
colvec<- c("#800020", "#E97451", "#C04000", "#7CFC00", "#8A2BE2","#FFD700")#, "#4169E1", "#7B68EE", "#D2691E", "#1E90FF", "#FF1493")
unique.df["Colvec"] <- colvec

## Merge both data frames
CLASS.meta <- merge(CLASS.meta,unique.df,
                  by = "Class", all = TRUE)
rownames(CLASS.meta) <- CLASS.meta$Row.names

## Order the matrix
mags.sorted <- as.vector(CLASS.meta$Row.names)
PULs <- mtx[match(mags.sorted,as.character(rownames(mtx))),]
PULs <- as.data.frame(PULs)
PULs <- tibble::rownames_to_column(PULs, "MAG")
PULs$MAG <- gsub("^X","",PULs$MAG)
PULs2 <- PULs[,-1]
rownames(PULs2) <- PULs[,1]
PULs <- PULs2
PULs <- PULs %>%
    replace(is.na(.), 0)

## Define the labels
class.labels <- CLASS.meta %>%
    distinct(Class) %>%
    unname()
class.labels <- c(class.labels)[[1]]
class.labels <- sort(class.labels)

## Plot the PCoA
plot.pcoa(mtx = PULs,metadata = CLASS.meta,labels = class.labels,component = "Class")


## ORDER DATASET
## Use the MAGs metadata and order that dataframe by 'Order' column
## Assign the MAGS.meta content to another variable
ORDER.meta <- MAGS.meta[is.element(MAGS.meta$Row.names,filtered.mags.vec),]
ORDER.meta <- ORDER.meta[order(ORDER.meta$Order,decreasing = FALSE),]
unique.df <- as.data.frame(unique(ORDER.meta$Order))
colnames(unique.df) <- c("Order")
#colvec <- c("#191970", "#00FFFF", "#088F8F", "#1434A4", "#CD7F32", "#800020", "#E97451", "#6E260E", "#8B0000", "#C04000", "#808000", "#AAFF00", "#228B22", "#7CFC00", "#008000", "#32CD32", "#0FFF50", "#FFBF00", "#9370DB", "#8A2BE2", "#9932CC", "#DA70D6", "#FF69B4", "#FF4500", "#FFD700", "#00CED1", "#800080", "#FF00FF", "#FF6347", "#FF8C00", "#4169E1", "#7B68EE", "#D2691E", "#1E90FF", "#FF1493", "#808080", "#A9A9A9", "#DC143C", "#00FF00", "#4169E1", "#FFFF00", "#00CED1", "#FF69B4", "#BDB76B", "#FF7F50", "#008080", "#DDA0DD", "#E6E6FA", "#ADFF2F", "#CD853F", "#00FFFF", "#FFD700", "#FF6347", "#C71585", "#000080", "#FF1493", "#0000FF", "#808080", "#FF4500", "#00FF00", "#008080", "#4682B4", "#2E8B57", "#FFFACD", "#FF69B4", "#FF4500", "#00FFFF", "#0000FF", "#FF1493", "#008080", "#808080", "#FFD700", "#00CED1", "#7B68EE", "#CD853F", "#ADFF2F", "#BDB76B", "#00FF00", "#808000", "#A9A9A9")
colvec <- c("#00FFFF", "#CD7F32", "#800020","#6E260E","#AAFF00", "#228B22", "#7CFC00", "#32CD32","#FFBF00", "#9370DB", "#9932CC", "#FFD700", "#D2691E")#, "#FF1493", "#808080", "#A9A9A9", "#DC143C", "#00FF00", "#4169E1", "#FFFF00", "#00CED1", "#FF69B4", "#BDB76B", "#FF7F50", "#008080", "#DDA0DD", "#E6E6FA", "#ADFF2F", "#CD853F", "#00FFFF", "#FFD700", "#FF6347", "#C71585", "#000080", "#FF1493", "#0000FF", "#808080", "#FF4500", "#00FF00", "#008080", "#4682B4", "#2E8B57", "#FFFACD", "#FF69B4", "#FF4500")
unique.df["Colvec"] <- colvec

## Merge both data frames
ORDER.meta <- merge(ORDER.meta,unique.df,
                    by = "Order", all = TRUE)
rownames(ORDER.meta) <- ORDER.meta$Row.names

## Order the matrix
mags.sorted <- as.vector(ORDER.meta$Row.names)
PULs <- mtx[match(mags.sorted,as.character(rownames(mtx))),]
PULs <- as.data.frame(PULs)
PULs <- tibble::rownames_to_column(PULs, "MAG")
PULs$MAG <- gsub("^X","",PULs$MAG)
PULs2 <- PULs[,-1]
rownames(PULs2) <- PULs[,1]
PULs <- PULs2
PULs <- PULs %>%
    replace(is.na(.), 0)

## Define the labels
order.labels <- ORDER.meta %>%
    distinct(Order) %>%
    unname()
order.labels <- c(order.labels)[[1]]
order.labels <- sort(order.labels)

## Plot the PCoA
plot.pcoa(mtx = PULs,metadata = ORDER.meta,labels = order.labels,component = "Order")


## FAMILY DATASET
## Use the MAGs metadata and order that dataframe by 'Family' column
## Assign the MAGS.meta content to another variable
FAM.meta <- MAGS.meta[is.element(MAGS.meta$Row.names,filtered.mags.vec),]
FAM.meta <- FAM.meta[order(FAM.meta$Family,decreasing = FALSE),]
unique.df <- as.data.frame(unique(FAM.meta$Family))
colnames(unique.df) <- c("Family")
#colvec <- c("#191970", "#00FFFF", "#088F8F", "#1434A4", "#CD7F32", "#800020", "#E97451", "#6E260E", "#8B0000", "#C04000", "#808000", "#AAFF00", "#228B22", "#7CFC00", "#008000", "#32CD32", "#0FFF50", "#FFBF00", "#9370DB", "#8A2BE2", "#9932CC", "#DA70D6", "#FF69B4", "#FF4500", "#FFD700", "#00CED1", "#800080", "#FF00FF", "#FF6347", "#FF8C00", "#4169E1", "#7B68EE", "#D2691E", "#1E90FF", "#FF1493", "#808080", "#A9A9A9", "#DC143C", "#00FF00", "#4169E1", "#FFFF00", "#00CED1", "#FF69B4", "#BDB76B", "#FF7F50", "#008080", "#DDA0DD", "#E6E6FA", "#ADFF2F", "#CD853F", "#00FFFF", "#FFD700", "#FF6347", "#C71585", "#000080", "#FF1493", "#0000FF", "#808080", "#FF4500", "#00FF00", "#008080", "#4682B4", "#2E8B57", "#FFFACD", "#FF69B4", "#FF4500", "#00FFFF", "#0000FF", "#FF1493", "#008080", "#808080", "#FFD700", "#00CED1", "#7B68EE", "#CD853F", "#ADFF2F", "#BDB76B", "#00FF00", "#808000", "#A9A9A9")
colvec <- c("#191970", "#00FFFF", "#088F8F", "#1434A4", "#CD7F32", "#800020", "#E97451", "#6E260E", "#8B0000", "#C04000", "#808000", "#AAFF00", "#228B22", "#7CFC00", "#008000", "#32CD32", "#0FFF50", "#FFBF00", "#9370DB", "#8A2BE2", "#9932CC", "#DA70D6", "#FF69B4", "#FF4500", "#FFD700", "#00CED1", "#800080", "#FF00FF", "#FF6347", "#FF8C00", "#4169E1", "#7B68EE", "#D2691E", "#1E90FF", "#FF1493", "#808080")#, "#A9A9A9", "#DC143C", "#00FF00", "#4169E1", "#FFFF00", "#00CED1", "#FF69B4", "#BDB76B", "#FF7F50", "#008080", "#DDA0DD", "#E6E6FA", "#ADFF2F", "#CD853F", "#00FFFF", "#FFD700", "#FF6347", "#C71585", "#000080", "#FF1493", "#0000FF", "#808080", "#FF4500", "#00FF00", "#008080", "#4682B4", "#2E8B57", "#FFFACD", "#32CD32")#, "#0FFF50", "#FFBF00", "#9370DB", "#8A2BE2", "#9932CC", "#DA70D6", "#FF69B4", "#FF4500", "#FFD700", "#00CED1", "#800080", "#FF00FF", "#FF6347", "#FF8C00", "#4169E1", "#7B68EE", "#D2691E", "#1E90FF", "#FF1493", "#808080", "#A9A9A9", "#DC143C", "#00FF00", "#4169E1", "#FFFF00", "#00CED1", "#FF69B4", "#BDB76B", "#FF7F50", "#008080", "#DDA0DD", "#E6E6FA", "#ADFF2F", "#CD853F", "#FF6347", "#008B8B")
unique.df["Colvec"] <- colvec

## Merge both data frames
FAM.meta <- merge(FAM.meta,unique.df,
                    by = "Family", all = TRUE)
rownames(FAM.meta) <- FAM.meta$Row.names

## Order the matrix
mags.sorted <- as.vector(FAM.meta$Row.names)
PULs <- mtx[match(mags.sorted,as.character(rownames(mtx))),]
PULs <- as.data.frame(PULs)
PULs <- tibble::rownames_to_column(PULs, "MAG")
PULs$MAG <- gsub("^X","",PULs$MAG)
PULs2 <- PULs[,-1]
rownames(PULs2) <- PULs[,1]
PULs <- PULs2
PULs <- PULs %>%
    replace(is.na(.), 0)

## Define the labels
family.labels <- FAM.meta %>%
    distinct(Family) %>%
    unname()
family.labels <- c(family.labels)[[1]]
family.labels <- sort(family.labels)

## Plot the PCoA
plot.pcoa(mtx = PULs,metadata = FAM.meta,labels = family.labels,component = "Family")


## Create the heatmap
# puls_heatmap <- t(mtx)
# 
# heatmap.df <- read.table(file = "1069_MAGs_metadata.tsv",sep = "\t",header = TRUE,row.names = 1)
# heatmap.df <- heatmap.df %>%
#     select(class,ecosystem_category)
# colnames(heatmap.df) <- c("Class","Ecosystem")
# heatmap.df <- filter(heatmap.df,!is.element(rownames(heatmap.df),c("3300005326_35","3300005370_1","3300017650_9","3300027284_65","3300027607_50")))
# rownames(heatmap.df) <- colnames(puls_heatmap)
#     # heatmap.df[!is.element(heatmap.df$MAG,c("3300005326_35","3300005370_1","3300017650_9","3300027284_65","3300027607_50")),]
# 
# ## Create the heatmap
# plot <- pheatmap(
#     mtx,
#     cluster_rows = FALSE,
#     cluster_cols = TRUE,
#     show_rownames = TRUE,
#     show_colnames = FALSE,
#     #annotation_col = heatmap.df
# )
# print(plot)
