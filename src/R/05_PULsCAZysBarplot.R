## Libraries
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))

## Functions
plot.boxplot <- function(metadata,component,type,colvec) {
    ## BOXPLOT
    ## Load the corresponding dataset
    name <- component
    if (component != "Ecosystem") {
        component <- paste0("Taxa_",component)
    }
    
    counts.df <- read.table(file = paste0(component,"/",type,".",name,"Diversity.tsv"),sep = "\t",header = FALSE)
    colnames(counts.df) <- c(name,"counts")
    counts.df[[name]] <- gsub("_"," ",counts.df[[name]])
    labels.arr <- rev(as.vector(counts.df[[name]]))
    labels.counts.arr <- rev(as.vector(paste0(counts.df[[name]],paste0("(",counts.df$counts,")"))))
    
    if (length(labels.arr) >= 60) {
        png(paste0(type,"_",name,".boxplot.png"), width=2500, height=4500, res = 500)
        plot <- ggplot(metadata,aes(x=factor(.data[[name]], levels = labels.arr),y=.data[[type]])) +
            geom_boxplot(aes(fill=.data[[name]])) +
            labs(x=name,y=paste0(type," Counts")) +
            scale_fill_manual(values = colvec) +
            geom_point(alpha = 0.1, position=position_jitter(height=0.005, width=0.005),size = 0.1) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "none",axis.text.y = element_text(size = 5)) +
            coord_flip() + 
            scale_x_discrete(breaks = labels.arr,
                             labels = labels.counts.arr) +
            theme(panel.grid.major = element_blank()) +
            annotate("text",  x=labels.arr[5], y = max(metadata[[type]]), label = paste0("n = ",sum(counts.df$counts)), vjust=1, hjust=1,size = 3)
        
        print(plot)
        
    } else {
        png(paste0(type,"_",name,".boxplot.png"), width=2500, height=3200, res = 500)
        plot <- ggplot(metadata,aes(x=factor(.data[[name]], levels = labels.arr),y=.data[[type]])) +
            geom_boxplot(aes(fill=.data[[name]])) +
            labs(x=name,y=paste0(type," Counts")) +
            scale_fill_manual(values = colvec) +
            geom_point(alpha = 0.1, position=position_jitter(height=0.01, width=0.01)) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "none") +
            coord_flip() + 
            scale_x_discrete(breaks = labels.arr,
                             labels = labels.counts.arr) +
            theme(panel.grid.major = element_blank()) +
            annotate("text",  x=labels.arr[5], y = max(metadata[[type]]), label = paste0("n = ",sum(counts.df$counts)), vjust=1, hjust=1,size = 3)
        
        print(plot)
    }

    dev.off()
}

## Set working directory
setwd('/home/chombo/Documents/SEGOVIA/CAZymes_PULs/data/')

## Load the MAGs metadata and extract the desired columns
MAGs.meta.df <- read.table(file = "1069_MAGs_metadata.tsv",sep = "\t",header = TRUE)
MAGs.meta.df <- MAGs.meta.df %>%
    select(genome_id,PULsFrequency,CAZysFrequency,ecosystem_category,class,order,family)
colnames(MAGs.meta.df) <- c("MAGs","PULS","CAZY","Ecosystem","Class","Order","Family")


## BOXPLOT
## CAZymes and PULs
colvec <- c(
    "#191970", "#00FFFF", "#088F8F", "#1434A4", "#800020", "#E97451", "#C04000", "#808000", "#AAFF00", "#008000", 
    "#32CD32", "#0FFF50", "#FFBF00", "#0ECAB2", "#960ECA", "#CAC00E", "#5F5E4B", "#ADFF08", "#FFAF00", "#B78700", 
    "#044F0B", "#E60F5E", "#830B37", "#E97451", "#9370DB", "#8A2BE2", "#9932CC", "#DA70D6", "#FF69B4", "#FF4500", 
    "#FFD700", "#00CED1", "#800080", "#FF00FF", "#FF6347", "#FF8C00", "#4169E1", "#7B68EE", "#D2691E", "#1E90FF", 
    "#FF1493", "#808080", "#A9A9A9", "#DC143C", "#00FF00", "#4169E1", "#FFFF00", "#00CED1", "#FF69B4", "#BDB76B", 
    "#FF7F50", "#008080", "#DDA0DD", "#E6E6FA", "#ADFF2F", "#CD853F", "#00FFFF", "#FFD700", "#FF6347", "#C71585", 
    "#000080", "#FF1493", "#0000FF", "#808080", "#FF4500", "#00FF00", "#008080", "#4682B4", "#2E8B57", "#FFFACD", 
    "#32CD32", "#0FFF50", "#FFBF00", "#9370DB", "#8A2BE2", "#9932CC", "#DA70D6", "#FF69B4", "#FF4500", "#FFD700", 
    "#00CED1", "#800080", "#FF00FF", "#FF6347", "#FF8C00", "#4169E1", "#7B68EE", "#D2691E", "#1E90FF", "#FF1493", 
    "#808080", "#A9A9A9", "#DC143C", "#00FF00", "#4169E1", "#FFFF00", "#00CED1", "#FF69B4", "#BDB76B", "#FF7F50", 
    "#008080", "#DDA0DD", "#E6E6FA", "#ADFF2F", "#CD853F", "#FF6347", "#008B8B",
    "#4B0082", "#8B008B", "#B8860B", "#8B4513", "#2E8B57", "#800000", "#8B0000", "#483D8B", "#2F4F4F", "#2F4F4F", 
    "#00CED1", "#5F9EA0", "#D2691E", "#A52A2A", "#8B4513", "#2F4F4F", "#6495ED", "#0000CD", "#00008B", "#FF8C00", 
    "#8B0000", "#483D8B", "#2F4F4F", "#2F4F4F", "#00CED1", "#5F9EA0", "#D2691E", "#A52A2A", "#8B4513", "#2F4F4F", 
    "#6495ED", "#0000CD", "#00008B", "#FF8C00"
)

## Ecosystem
plot.boxplot(metadata = MAGs.meta.df,component = "Ecosystem",type = "CAZY",colvec = colvec)
plot.boxplot(metadata = MAGs.meta.df,component = "Ecosystem",type = "PULS",colvec = colvec)

## Class
plot.boxplot(metadata = MAGs.meta.df,component = "Class",type = "CAZY",colvec = colvec)
plot.boxplot(metadata = MAGs.meta.df,component = "Class",type = "PULS",colvec = colvec)

## Order
plot.boxplot(metadata = MAGs.meta.df,component = "Order",type = "CAZY",colvec = colvec)
plot.boxplot(metadata = MAGs.meta.df,component = "Order",type = "PULS",colvec = colvec)

## Family
plot.boxplot(metadata = MAGs.meta.df,component = "Family",type = "CAZY",colvec = colvec)
plot.boxplot(metadata = MAGs.meta.df,component = "Family",type = "PULS",colvec = colvec)
