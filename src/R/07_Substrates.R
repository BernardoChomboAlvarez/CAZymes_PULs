## Libraries
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(textshape))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(readxl))


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

## PLOTS
## CAZYMES BARPLOT
plot.cazydiv <- function(all.cazys,component) {
    ## Obtain labels
    counts.df <- all.cazys %>%
        group_by(!!sym(component)) %>%
        summarize(total_frequency = sum(frequency, na.rm = TRUE))
    colnames(counts.df) <- c(component,"counts")
    #counts.df[[component]] <- gsub("_"," ",counts.df[[component]])
    labels.arr <- rev(as.vector(counts.df[[component]]))
    labels.counts.arr <- rev(as.vector(paste0(counts.df[[component]],paste0("(",counts.df$counts,")"))))
    
    
    if (component == "Class") {
        png(filename = paste0("CAZyFamilies_",component,"Diversity.png"),width = 3000,height = 1500,res = 350)
    } else {
        png(filename = paste0("CAZyFamilies_",component,"Diversity.png"),width = 3000,height = 2000,res = 350)
    }

    plot <- ggplot(data = all.cazys,
                   aes(fill = seed, x=factor(.data[[component]],levels = labels.arr), y=frequency),
                   show.legend = FALSE) +
        geom_bar(position="fill", stat="identity",show.legend = TRUE) +
        scale_y_continuous(labels = scales::percent) +
        theme_bw() +
        scale_fill_manual(values = c("#581845", "#900C3F", "#C70039", "#FF5733", "#FFC300","#DAF7A6")) +
        coord_flip() +
        ylab('Percentage') +
        xlab(component) +
        theme(legend.key.size = unit(0.2, 'cm'), #change legend key size
              legend.key.height = unit(0.2, 'cm'), #change legend key height
              legend.key.width = unit(0.2, 'cm'), #change legend key width
              legend.title = element_text(size=9), #change legend title font size
              legend.text = element_text(size=7),
              legend.position = 'right',
              legend.direction = 'vertical') +
        scale_x_discrete(breaks = labels.arr,
                         labels = labels.counts.arr) +
        guides(fill = guide_legend(ncol = 2,title = "CAZymes Families")) +
        labs(caption = paste0("n = ",sum(counts.df$counts)))
    
    print(plot)
    
    
    dev.off()
}

## SUBSTRATES DIVERSITY PLOT
plot.substrate <- function(metadata,type) {
    ## Create an array with the unique values of the 'Family' column
    fam.arr <- c(sort(unique(metadata$Family)))
    
    for (fam in fam.arr) {
        ## Create a temp dataset per family
        df.temp <- metadata %>%
            filter(Family == fam)
        
        ## Obtain labels
        counts.df <- df.temp %>%
            group_by(Family) %>%
            summarize(total_frequency = sum(frequency, na.rm = TRUE))
        colnames(counts.df) <- c("family","counts")
        labels.arr <- rev(c(counts.df$family))
        labels.counts.arr <- rev(c(paste0(counts.df$family,paste0("(",counts.df$counts,")"))))
        
        png(filename = paste0(type,"_",fam,"_SubstrateDiversity.png"),width = 3000,height = 800,res = 200)
        plot <- ggplot(data = df.temp,
                       aes(fill = substrate, x=factor(Family,levels = labels.arr), y=frequency),
                       show.legend = FALSE) +
            geom_bar(position="fill", stat="identity",show.legend = TRUE) +
            scale_y_continuous(labels = scales::percent) +
            theme_bw() +
            scale_fill_manual(values = c("#00FFCC", "#FF9900", "#33FF00", "#9900FF", "#FF0066", "#0066FF", "#CC3300", "#FFCC00", "#00FFCC", "#CC00FF", "#FF3399", "#99FF33", "#3399FF", "#FF6600", "#00FF66", "#6600FF", "#FF9966", "#66FF99", "#9966FF", "#99FF00", "#00FF33", "#0099FF", "#CC6600", "#00CC66", "#6600CC", "#FF3300", "#00FF00", "#0033CC", "#FF6666", "#66FF66", "#CCCCFF", "#FF3300", "#00FF00", "#0033CC", "#FF6666", "#66FF66", "#CCCCFF", "#CC6600", "#00CC66", "#6600CC", "#CC3300", "#00CC00", "#003399", "#FF9966", "#99FF66", "#9966FF", "#99FF00", "#00FF33", "#0099FF", "#663300", "#006633", "#330066", "#996633", "#339966", "#663399", "#666600", "#006666", "#660066", "#996666", "#669966", "#666699", "#669900", "#009966", "#990066", "#996699", "#669933", "#339933", "#333366", "#993333", "#336699", "#663366", "#333399", "#993366", "#336633", "#663333", "#333333", "#996600", "#006699", "#660099", "#333300", "#009933", "#990033", "#330099", "#339933", "#993333", "#336699", "#663366", "#333399", "#993366", "#336633", "#663333", "#333333", "#996600", "#006699", "#660099", "#333300", "#009933", "#990033", "#330099", "#33CC66", "#6633CC", "#CC6633", "#33CC99", "#9933CC", "#CC9933", "#33CC00", "#009933", "#3300CC", "#66CC33", "#CC3300", "#66CC66", "#CC66CC", "#66CC99", "#CC99CC", "#669900", "#00CC99", "#996633", "#33CC33", "#0099CC", "#9900CC", "#336600", "#00CC66", "#6633CC", "#33CC99", "#9933CC", "#CC9933", "#33CC00", "#009933", "#3300CC", "#66CC33", "#CC3300", "#66CC66", "#CC66CC", "#66CC99", "#CC99CC", "#669900", "#00CC99", "#996633", "#33CC33", "#0099CC", "#9900CC", "#336600", "#00CC66", "#6633CC", "#33CC99", "#9933CC", "#CC9933", "#33CC00", "#009933", "#3300CC", "#66CC33", "#CC3300", "#66CC66", "#CC66CC", "#66CC99", "#CC99CC", "#669900", "#00CC99", "#996633", "#33CC33", "#0099CC", "#9900CC", "#336600", "#00CC66", "#6633CC", "#33CC99", "#9933CC", "#CC9933", "#33CC00", "#009933", "#3300CC", "#66CC33", "#CC3300", "#66CC66", "#CC66CC", "#66CC99", "#CC99CC")) +
            coord_flip() +
            ylab('Percentage') +
            xlab('Family') +
            theme(legend.key.size = unit(0.2, 'cm'), #change legend key size
                  legend.key.height = unit(0.2, 'cm'), #change legend key height
                  legend.key.width = unit(0.2, 'cm'), #change legend key width
                  legend.title = element_text(size=9), #change legend title font size
                  legend.text = element_text(size=7),
                  legend.position = 'right',
                  legend.direction = 'vertical') +
            scale_x_discrete(breaks = labels.arr,
                             labels = labels.counts.arr) +
            guides(fill = guide_legend(ncol = 2,title = paste(type,fam,"Substrate",sep = " "))) +
            labs(caption = paste0("n = ",sum(counts.df$counts))) + 
            theme(legend.text=element_text(size=4))
        
        print(plot)
        dev.off()
    }
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
MAGS.meta["mag"] <- rownames(MAGS.meta)

## Read the 343 MAGs new list
MAGS.arr <- as.vector(read.table("343_MAGslist.txt"))[[1]]

## Filter MAGs.meta
MAGS.meta <- MAGS.meta[is.element(MAGS.meta$mag,MAGS.arr),]

## Extract PULs and CAZymes
all.puls <- get.puls(mags = MAGS.arr)
all.cazys <- get.cazymes(mags = MAGS.arr)

## Merge the metadata with the counts
all.cazys <- merge(all.cazys,MAGS.meta, by = "mag") %>%
    select(!c("PULS","CAZYS"))
all.puls <- merge(all.puls,MAGS.meta, by = "mag") %>%
    select(!c("PULS","CAZYS"))

## Set the seeds for CAZymes
all.cazys <- all.cazys %>%
    mutate(seed = gsub("\\d+","",cazy))
all.cazys$seed <- gsub("_","",all.cazys$seed)

# ================================================================================
# 18194 total CAZymes
# 17314 total PULs
# 343 total MAGs
# ================================================================================

## Plot CAZymes family diversity per category chosen
plot.cazydiv(all.cazys = all.cazys,component = "Ecosystem")
plot.cazydiv(all.cazys = all.cazys,component = "Class")
plot.cazydiv(all.cazys = all.cazys,component = "Order")
plot.cazydiv(all.cazys = all.cazys,component = "Family")

## SUBSTRATES
## PULS
## Load the PULs substrate dataset
puls.substrate <- read_excel("PULSUB_05102023.xlsx")
puls.substrate <- puls.substrate %>%
    select(PULID,substrate_final,cazymes_predicted_dbcan)
colnames(puls.substrate) <- c("pul","substrate","predictedCAZY")

## CLEAN MAGS BY THEIR SUBSTRATE
all.puls <- merge(all.puls,puls.substrate, by = "pul",all.x = TRUE) %>%
    select(!c("predictedCAZY"))
all.puls$substrate <- replace_na(all.puls$substrate,"Unknown")
all.puls <- all.puls %>%
    filter(substrate != "Unknown")

# ================================================================================
## Unknown PULs = 7812
## Left PULs = 7938
## MAGs = 343
## Substrates = 143 w/Unknown
# ================================================================================

substrates.unique <- all.puls %>%
    filter(substrate != "Unknown") %>%
    group_by(substrate) %>%
    summarize(counts = sum(frequency))
## summary(substrates.unique$counts)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00    7.00   24.50   66.92   74.75  651.00

substrates.unique <- substrates.unique %>%
    filter(counts >= 24)
substrates.arr <- c(substrates.unique$substrate)
all.puls <- all.puls[is.element(all.puls$substrate,substrates.arr),]

## Plot PULs substrate diversity
plot.substrate(metadata = all.puls, type = "PULs")


## CAZYMES
## Working with CAZymes is a little bit different because they have their own substrate binding module
## i.e., they are associated with an specific substrate. That is why the approach will be to group in
## a list all the possible substrates and substrate binding modules by each of the 175 available CAZymes
## in the 343 MAGs
## Load the CAZymes substrate dataset
cazymes.substrate <- read.table("CAZYSUB_08252022.tsv", sep = "\t", header = TRUE)
cazymes.substrate <- cazymes.substrate %>%
    select(Family,Substrate_high_level,Name)
colnames(cazymes.substrate) <- c("cazy","substrate","enzyme")

## CLEAN MAGS BY THEIR SUBSTRATE
## Create an array with te unique CAZymes of the 343 MAGs
all.cazys.copy <- all.cazys
all.cazys.copy$cazy <- gsub("_\\d+","",all.cazys.copy$cazy)

## Remove the GTs
all.cazys.copy <- all.cazys.copy %>%
    filter(!grepl("GT",cazy))
unique.cazys <- c(sort(unique(all.cazys.copy$cazy)))

## Create a list with the metadata of all the unique CAZymes
cazys.sub.list <- list()
cazy.sub.unknown <- c()
for (i in unique.cazys) {
    df.temp <- subset(cazymes.substrate,is.element(cazymes.substrate$cazy,i))
    cazys.sub.list[[i]] <- df.temp
    
    ## Store the CAZymes with unknown substrates for further analysis
    if (length(row.names(df.temp)) == 0) {
        cazy.sub.unknown <- c(cazy.sub.unknown,i)
    }
}
## Some CAZymes were not found in the substrate dataset 
## They were searched manually
cazys.sub.list[["PL27"]] <- data.frame(cazy = c("PL27","PL27"),
                                   substrate = c("L-rhamnose","D-gluconarate"),
                                   enzyme = c("L-rhamnose-α-1","4-D-glucuronate lyase"))
cazys.sub.list[["PL35"]] <- data.frame(cazy = c("PL35"),
                                   substrate = c("chondroitin"),
                                   enzyme = c("Chondroitin lyase"))
cazys.sub.list[["PL40"]] <- data.frame(cazy = c("PL40"),
                                   substrate = c("ulvan"),
                                   enzyme = c("Ulvan lyase"))

## Merge the CAZymes with the substrates
## Extract the substrates from each cazy data frame and merge them so they can be incorporated to the all.cazys data frame
trans.cazys.sub <- data.frame()
for (cazy in unique.cazys) {
    sub.vec <- sort(unique(cazys.sub.list[[cazy]]$substrate))
    sub.str <- paste(sub.vec,collapse = ",")
    trans.cazys.sub <- rbind(trans.cazys.sub,data.frame(cazy = cazy,substrate = sub.str))
}

## Merge the substrates with the all.cazys data frame
all.cazys.copy <- merge(all.cazys.copy,trans.cazys.sub,
                        by = "cazy",all.x = TRUE)

# ================================================================================
## GT CAZymes = 6647
## Left CAZymes = 11547
## MAGs = 343
## Substrates = 87
# ================================================================================

## Filter the substrates by their median
substrates.unique <- all.cazys.copy %>%
    group_by(substrate) %>%
    summarize(counts = sum(frequency))
# summary(substrates.unique$counts)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.0     9.5    54.0   132.7   137.0  1036.0

substrates.unique <- substrates.unique %>%
    filter(counts >= 54)
substrates.arr <- c(substrates.unique$substrate)
all.cazys.copy <- all.cazys.copy[is.element(all.cazys.copy$substrate,substrates.arr),]

## Plot the CAZymes substrate diversity
plot.substrate(metadata = all.cazys.copy, type = "CAZy")