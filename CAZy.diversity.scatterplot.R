library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(gtools)


filter_dfs <- function(df,family){
  df.temp <- df %>% 
    filter(grepl(family,CAZymes)) # warning! counts more than expected. Missmatches
  
  sums.temp <- sum(df.temp$Counts,NA,na.rm = TRUE)
  print(paste0("family ",family,": ",sums.temp))
  
  diversity.vector <- unique(df.temp$CAZymes)
  countsvec <- c()
  for (cazy in diversity.vector){
    temp <- df %>%
      filter(grepl(cazy,CAZymes))
    sums <- sum(temp$Counts,NA,na.rm = TRUE)
    countsvec <- c(countsvec,sums)
  }
  
  diversity.vector <- sapply(X = diversity.vector, FUN = function(t) gsub(pattern = "_", replacement = ".", x = t, fixed = TRUE))
  df.merge <- data.frame(diversity.vector,countsvec)
  df.merge <- df.merge[mixedorder(as.character(df.merge$diversity.vector)),]
  colnames(df.merge) <- c('CAZy','Sums')
  
  #print(df.merge)
  return(df.merge)
}

full_sums <- function(df){
  df$CAZymes <- gsub("_",".",df$CAZymes)
  diversity.vector <- unique(df$CAZymes)
  diversity.vector <- diversity.vector[diversity.vector != 'SLH']
  diversity.vector <- diversity.vector[diversity.vector != 'cohesin']
  diversity.vector <- sapply(X = diversity.vector, FUN = function(t) gsub(pattern = "_", replacement = ".", x = t, fixed = TRUE))
  
  seedvec <- c()
  idvec <- c()
  sumsvec <- c()
  meansvec <- c()
  countsvec <- c()
  for (cazy in diversity.vector){
    temp <- filter(df,CAZymes == cazy)
    counts <- nrow(temp)
    sums <- sum(temp$Counts,NA,na.rm = TRUE)
    means <- sums/counts
    seed <- str_split(cazy, '[1-9]+')[[1]][1]
    id <- str_split(cazy, '[A-Z]+')[[1]][2]
    
    seedvec <- c(seedvec,seed)
    idvec <- c(idvec,id)
    sumsvec <- c(sumsvec,sums)
    meansvec <- c(meansvec,means)
    countsvec <- c(countsvec,counts)
    
  }
  
  df.merge <- data.frame(diversity.vector,seedvec,idvec,countsvec,sumsvec,meansvec)
  df.merge <- df.merge[mixedorder(as.character(df.merge$diversity.vector)),][-1,]
  colnames(df.merge) <- c('CAZymes','Seed','ID','Counts','Sums','Means')
  
  #print(df.merge)
  return(df.merge)
}

setwd("C:/Users/chomb/Desktop/")
df <- read.table("CAZy.diversitylist.txt",sep = "\t",header = FALSE,row.names = NULL)
colnames(df) <- c('CAZymes','Counts','Taxa','Ecosystem')
famvec <- c('AA','GH','GT','CE','CBM','PL')
sum.full <- sum(df$Counts,NA,na.rm = TRUE)

# Test 1
#df.uniq.counts <- data.frame()

# for (item in famvec){
#   df.merge <- filter_dfs(df,item)
#   df.uniq.counts <- rbind(df.uniq.counts,df.merge)
# }
#sum.cleaned <- sum(df.uniq.counts$Sums,NA,na.rm = TRUE)

# Test 2 THIS FUNCTIONS !!!!
df.full <- full_sums(df)

df.full <- df.full %>% 
  replace(is.na(.), 0)
sums.full <- sum(df.full$Sums,NA,na.rm = TRUE)

png(filename = "CAZy.sums.png",width = 2000, height = 2700,res = 400)
plot <- ggplot(data = df.full,
               aes(x=factor(df.full$Seed),
                   y=as.numeric(df.full$ID)),show.legend = FALSE) +
  geom_point(aes(size = Sums,
                 colour = factor(Seed)),na.rm = TRUE) +
  scale_y_continuous(name = "Subtypes",
                     limits = c(0,170),
                     breaks = seq(0, 170, by = 10)) +
  xlab('CAZymes families') +
  theme_bw() +
  labs(size="Frequency") +
  guides(color = FALSE) +
  theme(legend.position = c(0.9, 0.91),
        legend.background = element_rect(fill = "white", color = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"))

print(plot)
dev.off() #cohesin and SHL not included (freq. of 10 between them)


