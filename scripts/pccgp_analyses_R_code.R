# Genotypic characterization of the U.S. peanut core collection				
# Authors: Paul I. Otyama??????1, Roshan Kulkarni?????????1, Kelly Chamberlain§, Peggy Ozias-Akins**, Ye Chu**,
# Lori Lincoln??????, Gregory E. MacDonald??????, Noelle Anglin§§, Sudhansu Dash***, Michelle Graham?????????, 
# Steven B. Cannon?????????, Ethalinda K.S. Cannon??????	

# R code that generated Main Figures 3,4,5 and File S10.

#.............
# Figure 3: Genetic structure of 518 samples selected as representatives at >= 98% sequence identity.
# Accessions are grouped into five clusters represented by distinct colors. The X-axis represents 
# accessions ordered according to their positions in the phylogenetic tree analysis. 
# The Y-axis represents proportions of cluster assignment based on Q values from fastStructure analysis.

# Step 1
# Run fastSTRUCTURE analysis to generate meanQ files. Refer to MS for details                                                               ##  
# install pophelper package from github                                                                          ##
# http://www.royfrancis.com/pophelper/articles/index.html                                                        ##
# BiocManager::install("devtools")                                                                               ##
# devtools::install_github('royfrancis/pophelper') 

# load libraries
library(pophelper)
library(ggplot2)
library(dplyr)
library(VennDiagram)

# Step 2
# transfer all meanQ files into a single directory 
setwd("/Users/piotyama/OneDrive/peanut/03_str/cen98/")
path="/Users/piotyama/OneDrive/peanut/03_str/cen98/"

# load my Q files to qlist
afiles <- dir(path=path, pattern="*.meanQ",recursive=FALSE)
alist <- readQ(afiles, filetype="basic")
head(alist[[2]])

# Add custom individual labels
inds2 <- read.delim("/Users/piotyama/OneDrive/peanut/05_metaData/structureindlabels_cen98.txt",header=T, stringsAsFactors=FALSE)

# change rownames to suitable custom labels
rownames(alist[[8]]) <- inds2$labels # K5
rownames(alist[[7]]) <- inds2$labels # K4

# step 3: Plot K4 & K5 Q results
# a) Sort by cluster
plotQMultiline(alist[8], useindlab=T, spl=130, lpp=4, showyaxis=T,showticks=T,indlabsize=6,
               sortind= "all",height=15,barsize=0.8, imgtype="pdf",
               showlegend=T,legendkeysize=8,legendtextsize=10,
               dpi=600,outputfilename="cen98-K5-sorted-by-color") # K5 sortind=alll

plotQMultiline(alist[7], useindlab=T, spl=130, lpp=4, showyaxis=T,showticks=T,indlabsize=6,
               sortind= "all",height=15,barsize=0.8, imgtype="pdf",
               showlegend=T,legendkeysize=8,legendtextsize=10,
               dpi=600,outputfilename="cen98-K4-sorted-by-color") # K4 sortind=alll

# b) Sort by custom (tree order ie same order as input file) 
plotQMultiline(alist[8], useindlab=F, spl=130, lpp=4, showyaxis=T,showticks=T,indlabsize=6,
               sortind= "label",height=15,barsize=0.8, imgtype="pdf",
               dpi=600,outputfilename="cen98-K5-tree-order") # K5 useindlab=F

plotQMultiline(alist[7], useindlab=F, spl=130, lpp=4, showyaxis=T,showticks=T,indlabsize=6,
               sortind= "label",height=15,barsize=0.8, imgtype="pdf",
               dpi=600,outputfilename="cen98-K4-tree-order") # K4 useindlab=F

# c) Add Groups if dsired. Using phylogenetic clade groups as an example and Sort by groups
clades <- inds2[,3,drop=F]
clades$clades <- as.character(clades$clades)

plotQMultiline(alist[8], useindlab=F,spl=130, lpp=4, showyaxis=T,showticks=T,indlabsize=5,
               height=15,barsize=0.8, 
               grplab = clades, ordergrp=TRUE, grplabsize = 1.0, grplabbgcol = "lightcyan2",grplabcol = "black",
               imgtype="pdf", dpi=600,outputfilename="cen98_K5_tree_order_w_clades_final") # K 5 useindlab=T

plotQMultiline(alist[7], useindlab=T,spl=130, lpp=4, showyaxis=T,showticks=T,indlabsize=6,
               sortind= "label",height=10,barsize=0.9, 
               grplab = clades, ordergrp=TRUE, grplabsize = 4.0, grplabbgcol = "lightcyan2",grplabcol = "black",
               imgtype="png", dpi=600,outputfilename="cen98_K4_tree_order_w_clades_final") # K 4 useindlab=T

## END
#.............

# File S10 (figure) [SF10_K5_membership.pdf] shows the proportion of accessions 
# assigned to clusters 1-5 in a Structure analysis (Figure 3), for K=5 clusters.

# Data
count.data <- data.frame(
  subpops = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"),
  n = c(166, 164, 19, 112, 57),
  prop = c(32.04, 31.66, 3.67, 21.62, 11.00)
)
count.data

# Add label position
count.data <- count.data %>%
  arrange(desc(subpops)) %>%
  mutate(lab.ypos = cumsum(prop) - 0.5*prop)
count.data

mycols <- c("chartreuse3", "#c6c6ff", "blue", "red", "yellow")

mem <- ggplot(count.data, aes(x = "", y = prop, fill = subpops)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  geom_text(aes(y = lab.ypos, label = n), color = "black", size = 6, angle = 10)+
  scale_fill_manual(values = mycols) +
  theme_void()

# legend title and labels
plot <- mem + theme(legend.title = element_text(size=14, face="bold")) + 
  theme(legend.text = element_text(size=14)) + 
  theme(legend.position=c(1.1, 0.6))

ggsave(plot=plot,"K5_membership.pdf", device = "pdf",dpi = "print", height = 4, width = 8, unit = "in", scale = 1.5)

## END
#.............

# Figure 4 Principal Component Analysis of 1120 samples based on 2063 unlinked SNP markers. 
# The X-axis represents PC 1 and the Y-axis represents PC 2. Samples are colored 
# and grouped according to: A. clade membership as defined in the phylogenetic and network analyses, 
# B. botanical varieties, C. market type, D. growth Habit, E. pod shape, and F. collection type.

# STEP 1: Load the R packages
library(gdsfmt)
library(SNPRelate)
library(SeqArray)
library(dplyr)
library(MASS)
library(tidyr)
library(ggplot2)
library(reshape2)
library(cowplot)
library(factoextra)
library(FactoMineR)

setwd("C:/Users/Polight/OneDrive/peanut/09_Genotypic_analyses")

# STEP 2: Data preparation
# The SNPRelate package provides a function snpgdsVCF2GDS() to reformat a VCF file
vcf.fn <- "aradu_araip_synth.div.HK28.main_01212020.vcf"

# Reformat
snpgdsVCF2GDS(vcf.fn, "test.gds", method="biallelic.only")

# Summary
snpgdsSummary("test.gds")
# source vcf: Double check that source is updated to the current version of vcf
# File has all 1120 samples with 13410 SNPs

# Principle Component Analysis
# Open the GDS file
genofile <- snpgdsOpen("test.gds")

# LD-based SNP pruning
# It is suggested to use a pruned set of SNPs which are in approximate linkage equilibrium with 
# each other to avoid the strong influence of SNP clusters in principal component analysis 
# and relatedness analysis.

set.seed(1000)

# Try different LD thresholds for sensitivity analysis
snpset <- snpgdsLDpruning(genofile, autosome.only = F, remove.monosnp = T, maf = 0.05,
                          slide.max.bp = 1000000, method = "corr", ld.threshold=0.2, verbose = T) 

# NOTE: Recursively removes biallelic SNPs within a sliding window of 1mbp in LD at the specified threshold. 
# Filter also for mono morphic SNPs and MAF 5%
# RESULT: 2,050 markers are selected in total

# Extract sample id and variant id
head(samp.id <- seqGetData(genofile, "sample.id"))
head(variant.id <- seqGetData(genofile, "variant.id"))

# select only snps that were selected from LD pruning using snpgdsLDpruning
seqSetFilter(genofile, variant.sel = snpset.id)

# Export this as a vcf file for other downstream analyses requiring unlinked markers eg fst calculation
vcf.fn2 <- "pruned_snpdata.HK28.main_01212020.vcf"
seqGDS2VCF(genofile, vcf.fn2, seqSetFilter(genofile, variant.sel = snpset.id))

# STEP 3: Run PCA
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=4, autosome.only = F)

# calculate the percent of variation accounted for by the top principal components
# variance proportion (%)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

# STEP 4: Get population annotation and order ID by sample ID
pop_code <- read.delim("population.groups.annotations_sorted.txt", na.strings = "#N/A")
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id")) # Get sample id

# make a merged table of eigenvectors and population annotions
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],
                  EV4 = pca$eigenvect[,4],
                  EV5 = pca$eigenvect[,5],
                  EV6 = pca$eigenvect[,6],
                  EV7 = pca$eigenvect[,7],
                  EV8 = pca$eigenvect[,8],
                  EV9 = pca$eigenvect[,9],
                  EV10 = pca$eigenvect[,10],
                  stringsAsFactors = FALSE)

write.csv(tab3,"pca_masterfile_results_annot.csv") # write results to file plotting as desired
## END
#...................

# Figure 5 Plots of FST (fixation index) values among genetic groupings, 
# to determine stratification in the core collection. Cluster identities are as shown in the 
# phylogenetic and PCA analyses. The pairwise population differentiation (FST index) 
# was calculated using Hierfstat for a set of unlinked markers and plotted as heatmaps. 
# Accessions were classified into groups of: A. clade membership as defined in the phylogenetic 
# and network analyses, B. botanical varieties, C. market type, D. growth Habit, E. pod shape, 
# and F. continent of seed origin.

# STEP 1: Load required libraries
library(hierfstat)
library(pegas)
library(adegenet)
library(SeqArray)
library(ggplot2) 
library(ggcorrplot) 
library(RColorBrewer) 
library(ggpubr) 
library(egg) # key for ggarange with nice alignments

setwd("/Users/piotyama/OneDrive/peanut/09_Genotypic_analyses")

# STEP 2: Import Unlinked SNPsets for Fst analysis
d <- read.vcf("pruned_snpdata.HK28.main_01212020.vcf", which.loci = 1:1e5) # specify a number greater loci present such that all are read
ind.desc <- read.csv("sample.id.csv", sep = ",", header = T)
ind.desc <- read.csv("sample.id.vcf.csv", sep = ",", header = T)
str(d)
str(ind.desc)

# STEP 3: Transform data into hierfstat df
dat <- genind2hierfstat(loci2genind(d), pop=ind.desc$Continents)
dat[1:25,1:10] 
sum(is.na(dat$pop)) # pop lable cannot be NA! 
dat <- dat[!is.na(dat$pop),] # Remove the two cases where pop is defined as NA
dat <- dat[!dat$pop == "Unknown",] # Remove all cases where pop is defined as unknown
wcfst_continents <- pairwise.WCfst(dat,diploid=TRUE) # pairwise WC Fst for continental divisions

# New identifiers to pop
dat_clade <- genind2hierfstat(loci2genind(d), pop=ind.desc$clades)
dat_clade[1:15,1:10] 
sum(is.na(dat_clade$pop))
dat_clade <- dat_clade[!is.na(dat_clade$pop),] # Remove the two cases where pop is defined as NA
wcfst_clades <- pairwise.WCfst(dat_clade, diploid=TRUE) # pairwise WC Fst for clade membership
#-----------------

dat_mstf <- genind2hierfstat(loci2genind(d), pop=ind.desc$mainstem_flower)
dat_mstf[1:15,1:10]
dat_mstf <- dat_mstf[!is.na(dat_mstf$pop),] # Remove the two cases where pop is defined as NA
dat_mstf <- dat_mstf[!dat_mstf$pop == "Unknown",] # Remove all cases where pop is defined as unknown
wcfst_mstf <- pairwise.WCfst(dat_mstf, diploid=TRUE) # pairwise WC Fst for main stem flower groups
#-----------------

dat_mkt <- genind2hierfstat(loci2genind(d), pop=ind.desc$Market_type)
dat_mkt[1:15,1:10]
dat_mkt <- dat_mkt[!is.na(dat_mkt$pop),] # Remove the two cases where pop is defined as NA
dat_mkt <- dat_mkt[!dat_mkt$pop == "Unknown",] # Remove all cases where pop is defined as unknown
wcfst_mkt <- pairwise.WCfst(dat_mkt, diploid=TRUE) # pairwise WC Fst for main stem flower groups
#-----------------

dat_ghh <- genind2hierfstat(loci2genind(d), pop=ind.desc$Growth_habit_Holbrrok)
dat_ghh[1:15,1:10]
dat_ghh <- dat_ghh[!is.na(dat_ghh$pop),] # Remove the two cases where pop is defined as NA
dat_ghh <- dat_ghh[!dat_ghh$pop == "Unknown",] # Remove all cases where pop is defined as unknown
wcfst_ghh <- pairwise.WCfst(dat_ghh, diploid=TRUE) # pairwise WC Fst for main stem flower groups
#-----------------

dat_bvg <- genind2hierfstat(loci2genind(d), pop=ind.desc$Botanical_variety_GRIN)
dat_bvg[1:15,1:10]
dat_bvg <- dat_bvg[!is.na(dat_bvg$pop),] # Remove the two cases where pop is defined as NA
dat_bvg <- dat_bvg[!dat_bvg$pop == "Unknown",] # Remove all cases where pop is defined as unknown
wcfst_bvg <- pairwise.WCfst(dat_bvg, diploid=TRUE) # pairwise WC Fst for main stem flower groups
#-----------------

dat_psh <- genind2hierfstat(loci2genind(d), pop=ind.desc$Pod_shape_Holbrook)
dat_psh[1:15,1:10]
dat_psh <- dat_psh[!is.na(dat_psh$pop),] # Remove the two cases where pop is defined as NA
dat_psh <- dat_psh[!dat_psh$pop == "Unknown",] # Remove all cases where pop is defined as unknown
wcfst_psh <- pairwise.WCfst(dat_psh, diploid=TRUE) # pairwise WC Fst for main stem flower groups
#-----------------

# save pairwise fst results for all 7 groups for future reference and plot as desired
# STEP 4: Reporting results as heatmaps

wcfst_clades[is.na(wcfst_clades)] <- 0
A <- ggcorrplot(wcfst_clades, type = "lower",  ggtheme = theme_half_open(), legend.title = "Fst", 
                colors = c("blue",  "honeydew",  "red"), hc.order = TRUE, lab = TRUE, digits = 2, show.diag = TRUE, show.legend = FALSE)

wcfst_bvg[is.na(wcfst_bvg)] <- 0
B <- ggcorrplot(wcfst_bvg, type = "lower",  ggtheme = theme_half_open(), legend.title = "Fst", 
                colors = c("blue",  "honeydew",  "red"), hc.order = TRUE, lab = TRUE, digits = 2, show.diag = TRUE, show.legend = FALSE)

wcfst_mkt[is.na(wcfst_mkt)] <- 0
C <- ggcorrplot(wcfst_mkt, type = "lower",  ggtheme = theme_half_open(), legend.title = "Fst", 
                colors = c("blue",  "honeydew",  "red"),  hc.order = TRUE, lab = TRUE, digits = 2, show.diag = TRUE, show.legend = FALSE)

wcfst_ghh[is.na(wcfst_ghh)] <- 0
D <- ggcorrplot(wcfst_ghh, type = "lower",  ggtheme = theme_half_open(), legend.title = "Fst", 
                colors = c("blue",  "honeydew",  "red"),  hc.order = TRUE, lab = TRUE, digits = 2, show.diag = TRUE, show.legend = FALSE)

wcfst_psh[is.na(wcfst_psh)] <- 0
E <- ggcorrplot(wcfst_psh, type = "lower",  ggtheme = theme_half_open(), legend.title = "Fst", 
                colors = c("blue",  "honeydew",  "red"),  hc.order = TRUE, lab = TRUE, digits = 2, show.diag = TRUE, show.legend = FALSE)

wcfst_continents[is.na(wcfst_continents)] <- 0
F <- ggcorrplot(wcfst_continents, type = "lower",  ggtheme = theme_half_open(), legend.title = "Fst", 
                colors = c("blue",  "honeydew",  "red"), hc.order = TRUE, lab = TRUE, digits = 2, show.diag = TRUE, show.legend = TRUE) # Show legend

# Use ggarrange from egg for nice alignments and save using ggsave as pdf
fsts_a_f <- ggarrange(A,B,C,D,E,F, labels = c("A", "B", "C", "D", "E", "F"), nrow = 2) # no main stem flower in analysis
ggsave(plot=fsts_a_f,"fst_heatmaps_final.pdf", device = "pdf",dpi = "print", height = 4, width = 8, unit = "in", scale = 2.0)
ggsave(plot=fsts_a_f,"Fig7_fst.pdf", device = "pdf",dpi = "print", height = 4, width = 8, unit = "in", scale = 2.0) # paper

## END


