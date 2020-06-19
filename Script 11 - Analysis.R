####################################################
#####    PART 1: loading required libraries    #####
####################################################

library(ggplot2)
library(tidyverse)
library(plotly)
library(readr)
library(stringr)
library(mudata)



#####################################################
#####           PART 2: loading data             ####
#####################################################

# Collecting phylogeny data
setwd("C:/Users/elisc/AppData/Local/Packages/CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc/LocalState/rootfs/home/eli/BEP/Genomes_of_interest/Anvio")
main.df1 <- read.csv('Additional_info.txt',sep='\t')

# Collecting BUSCO data
setwd("C:/Users/elisc/AppData/Local/Packages/CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc/LocalState/rootfs/home/eli/BEP/Rhizobiales/BUSCO_output")
Rhi.BUSCO <- read.delim('BUSCO_summary.txt')
setwd("C:/Users/elisc/AppData/Local/Packages/CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc/LocalState/rootfs/home/eli/BEP/MEGARHO/BUSCO_output")
Rho.BUSCO <- read.delim('BUSCO_summary.txt')
BUSCO.all <- rbind(Rho.BUSCO,Rhi.BUSCO)

# Collecting GTA data
setwd("C:/Users/elisc/AppData/Local/Packages/CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc/LocalState/rootfs/home/eli/BEP/Rhizobiales/GTA_results")
Rhi.GTA_results <- read.csv('GTA_clustering_results.csv')
setwd("C:/Users/elisc/AppData/Local/Packages/CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc/LocalState/rootfs/home/eli/BEP/MEGARHO/GTA_results")
Rho.GTA_results <- read.csv('GTA_clustering_results.csv')
GTA.all <- rbind(Rho.GTA_results,Rhi.GTA_results)
GTA.all$no.of.homologs[GTA.all$no.of.homologs > 11] <- 11
GTA.all <- select(GTA.all, 1:3)

temp1 <- main.df1[main.df1$GTA == 'False',]
temp2 <- data.frame('Assembly.Accession' = temp1$genome_id, 'no.of.homologs' = rep(0,length(temp1$genome_id)), 'GTA.Cluster' = rep('False',length(temp1$genome_id)))

GTA.all <- rbind(GTA.all,temp2)

# Collecting CRISPR & Cas data
setwd("C:/Users/elisc/AppData/Local/Packages/CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc/LocalState/rootfs/home/eli/BEP/Genomes_of_interest/CRISPRCasFinder_results")
CRISPRCAS <- read_delim("C:/Users/elisc/AppData/Local/Packages/CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc/LocalState/rootfs/home/eli/BEP/Genomes_of_interest/CRISPRCasFinder_results/TOTAL.CRISPR.summary.tsv", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)

# Merging data frames
colnames(main.df1)[1] <- 'Assembly.Accession'
main.df1 <- merge(select(main.df1, 1:5), BUSCO.all, by.x = 'Assembly.Accession')
main.df1 <- merge(main.df1, select(GTA.all, 1:3), by.x = 'Assembly.Accession', all.x = TRUE)
main.df1 <- merge(main.df1, CRISPRCAS, all.x = TRUE)

# Making sure NAs are converted
main.df1$no.of.homologs[is.na(main.df1$no.of.homologs)] <- 0
main.df1$GTA.Cluster[is.na(main.df1$GTA.Cluster)] <- 'False'




####################################################
#####          PART 4: making subsets          #####
####################################################

# subset with complete CRISPR 
CRISPR.subset <- subset(main.df2, subset = Complete.CRISPR == 'True')

# subset with only True or False GTA clusters (disregarding completeness in this plot yields a better interpretable result)
GTA.subset <- subset(main.df2, subset = Complete.CRISPR == 'TRUE')
GTA.subset <- subset(GTA.subset, subset = GTA.Cluster == 'True' | GTA.Cluster == 'False')



####################################################
#####       PART 5: Analysing the results      #####
####################################################

# BUSCO score distribution
hist.BUSCO <- ggplot(BUSCO.all, aes(x=BUSCO.score))
hist.BUSCO <- hist.BUSCO +
  geom_histogram(stat='count', binwidth = 1) +
  ggtitle('BUSCO score distribution') +
  theme_minimal() +
  geom_vline(xintercept=80, linetype='dashed')
print(hist.BUSCO)

# distribution of GTA genes as histogram
hist.GTA1 <- ggplot(GTA.all, aes(x=no.of.homologs)) +
  geom_histogram(binwidth = 1) +
  scale_x_continuous(breaks=0:11) +
  ylab('Genome count') + xlab('Number of GTA homologs in a genome') +
  ggtitle('Distribution of GTA homologs') +
  theme_minimal()
print(hist.GTA1) 

# distribution of CRISPR systems as double histogram
hist.GTA2 <- ggplot(GTA.subset, aes(x=Nb.of.unique.spacers, fill=GTA.Cluster)) +
  geom_histogram(binwidth = 1, alpha=0.6, position = 'identity') +
  coord_cartesian(xlim = c(0,50)) +
  ggtitle('Spacer count depending on GTA presence') + 
  theme_minimal() +
  theme(legend.position="bottom")
print(hist.GTA2)

# histogram showing the distribution of CRISPR systems in Rhizobiales
hist.CRISPR.dist <- ggplot(Rho.CRISPR, aes(x=CRISPR.Evidence_level, y = (..count..)/sum(..count..))) +
  geom_histogram(binwidth = 1, alpha=0.6, position = 'identity') +
  ggtitle('CRISPR evidence level distribution in Rhodobacterales') + 
  scale_y_continuous(labels = scales::percent) +
  theme_minimal()
print(hist.CRISPR.dist)

# histogram showing the occurence of different CAS systems in Rhodobacterales
hist.CAS <- ggplot(Rho.CAS, aes(x = CAS.type)) +
  geom_bar() +
  coord_flip() +
  ylab('CRISPR system') +
  theme_minimal()
print(hist.CAS)


#TEMP#
setwd("C:/Users/elisc/AppData/Local/Packages/CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc/LocalState/rootfs/home/eli/BEP/MEGARHO/CRISPRCasFinder_results")
MEGARho.CAS <- read.csv('TOTAL.CAS.summary.csv')
hist.CAS2 <- ggplot(MEGARho.CAS, aes(x = CAS.type)) +
  geom_bar() +
  coord_flip() +
  theme_minimal()
print(hist.CAS2)




Subset1 <- main.df2[which(main.df2$GTA.Cluster == 'True'),]
Subset2 <- main.df2[which(main.df2$GTA.Cluster == 'False'),]
final_vector_1 <- c()
final_vector_2 <- c()
for (i in test1$System) { 
  final_vector_1 <- c( final_vector_1, str_match(strsplit(i, ",")[[1]], "[A-Z].*[A-Z]"))
}
table(final_vector_1)
for (i in test2$System) { 
  final_vector_2 <- c( final_vector_2, str_match(strsplit(i, ",")[[1]], "[A-Z].*[A-Z]"))
}
table(final_vector_2)
barplot(table(final_vector_1))
barplot(table(final_vector_2))







test1 <- main.df2[which(main.df2$Complete.CRISPR == TRUE & main.df2$GTA.Cluster == 'True'),]
test2 <- main.df2[which(main.df2$Complete.CRISPR == TRUE & main.df2$GTA.Cluster == 'False'),]
final_vector_1 <- c()
final_vector_2 <- c()
for (i in test1$System) { 
  final_vector_1 <- c( final_vector_1, str_match(strsplit(i, ",")[[1]], "[A-Z].*[A-Z]"))
}
rename.values(final_vector_1, 'CAS-TypeIA'="CRISPR TypeIA", 'CAS-TypeIC'="CRISPR TypeIC", "CAS-TypeIE"="CRISPR TypeIE", "CAS-TypeIIU"="CRISPR TypeIIU")
table(final_vector_1)
for (i in test2$System) { 
  final_vector_2 <- c( final_vector_2, str_match(strsplit(i, ",")[[1]], "[A-Z].*[A-Z]"))
}
table(final_vector_2)
tempdf3 <- merge(tempdf1,tempdf2,all.y = T)
barplot(table(final_vector_1))
barplot(table(final_vector_2))
