### R script to plot statistics generated from a VCF file
### code is from https://speciationgenomics.github.io/filtering_vcfs/
### thanks and citations go to Mark Ravinet and Joana Meier, as well as the excellent course that they ran!
rm(list=ls())
setwd("yourdirectory")
#load packages
library(tidyverse)
library(ggplot2)
#set theme for graphs
my_theme <- theme_light(base_size = 22) + theme(
  plot.title = element_text(color="black", size=22, face="plain"),
  axis.title.x = element_text(color="black", size=22, face="plain"),
  axis.title.y = element_text(color="black", size=22, face="plain")
) + theme(plot.title = element_text(hjust = 0.5))

#plot variant quality, recommend phred >30
var_qual <- read_delim("youroutfile.lqual", delim = "\t",
                       col_names = c("chr", "pos", "qual"), skip = 1)
a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  my_theme + ggtitle("Phred quality score") +  xlab("Quality") + ylab("Density") + ylim(0,0.008)
a 
#plot read depth
var_depth <- read_delim("./youroutfile.ldepth.mean", delim = "\t",
                        col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  my_theme + ggtitle("Read depth") +  xlab("Mean depth") + ylab("Density") + xlim(0,15)
a
summary(var_depth$mean_depth)
#plot variant missingness
#note that vcftools inverts missingness
var_miss <- read_delim("./youroutfile.lmiss", delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  my_theme + ggtitle("Missingness") +  xlab("Missing") + ylab("Density")
a 
summary(var_miss$fmiss)
#plot allele frequency - to help identify appropriate MAF
var_freq <- read_delim("./youroutfile.frq", delim = "\t",
                      col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))
a <- ggplot(var_freq, aes(maf)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  my_theme + ggtitle("Distribution of MAF") +  xlab("MAF") + ylab("Density")
a 
summary(var_freq$maf)
#plot individual depth
ind_depth <- read_delim("./youroutfile.idepth", delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)
a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  my_theme + ggtitle("Mean depth per individual") +  xlab("depth") + ylab("count") 
a 
ind_depth
#plot proportion of missing data
ind_miss  <- read_delim("./youroutfile.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  my_theme + ggtitle("Proportion missing data") +  xlab("proportion missing") + ylab("count") 
a 
ind_miss
#plot heterozygosity and inbreeding coefficient
ind_het <- read_delim("./youroutfile.het", delim = "\t",
                      col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)
a <- ggplot(ind_het, aes(f)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  my_theme + ggtitle("Inbreeding coefficient") +  xlab("inbreeding coefficient") + ylab("count") 
a
ind_het
