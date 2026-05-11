# Crotalus genomic landscape, RAxML & SVDquartets
Author: Keaka Farleigh, Ph.D.
Date: March 31st, 2025
Email: keakafarleigh@virginia.edu

## Purpose
This script will use [RAxML](https://bioconda.github.io/recipes/raxml/README.html) to build a maximum-likelihood phylogeny.

## Location 

This analysis was performed in the Schield lab on Xenomorph at the University of Virginia. Xenomorph is a system76 workstation running Ubuntu v22.04.3.


## Input data
This script uses a phylip file generated using Dr. Edgardo Ortiz's `[vc2phylip](https://github.com/edgardomortiz/vcf2phylip)` script.

### Set up environment
First, we set up the directory structure and our environment.

```
conda activate raxml

mkdir /media/queen/extradrive1/crotalus_genomic_landscape/analysis/raxml

cd /media/queen/extradrive1/crotalus_genomic_landscape/analysis/raxml

```

### Filter vcf 
We will filter the vcf to contain no missing data. We will not thin SNPs based on the recommendation of Dr. Laura Kubatko, the senior author of SVDquartets.

```
# For everything

nohup bcftools view --threads 24 -m2 -M2 -U -v snps -i 'F_MISSING=0.0' -O z -o /media/queen/extradrive1/crotalus_genomic_landscape/vcf/crotalus_genus.noCV13.auto+chrZ-tmp.snps.focal.vcf.gz /media/queen/extradrive1/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.final.noCV13.auto+chrZ.snps.vcf.gz > ./log/raxml_allfilter.log &
nohup vcftools --gzvcf /media/queen/extradrive1/crotalus_genomic_landscape/vcf/crotalus_genus.noCV13.auto+chrZ-tmp.snps.focal.vcf.gz  --recode --stdout | bgzip -c > /media/queen/extradrive1/crotalus_genomic_landscape/vcf/crotalus_genus.noCV13.auto+chrZ.snps.complete.vcf.gz > ./log/raxml_all.log &

# Autosome only

nohup bcftools view --threads 12 -m2 -M2 -U -v snps -i 'F_MISSING=0.0' -O z -o /media/queen/extradrive1/crotalus_genomic_landscape/vcf/crotalus_genus.noCV13.auto-tmp.snps.focal.vcf.gz /media/queen/extradrive1/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.final.noCV13.auto.snps.vcf.gz > ./log/raxml_autofilter.log &
nohup vcftools --gzvcf /media/queen/extradrive1/crotalus_genomic_landscape/vcf/crotalus_genus.noCV13.auto-tmp.snps.focal.vcf.gz --recode --stdout | bgzip -c > /media/queen/extradrive1/crotalus_genomic_landscape/vcf/crotalus_genus.noCV13.auto.snps.complete.vcf.gz > ./log/raxml_auto.log &

# Z only

nohup bcftools view --threads 12 -m2 -M2 -U -v snps -i 'F_MISSING=0.0' -O z -o /media/queen/extradrive1/crotalus_genomic_landscape/vcf/crotalus_genus.noCV13.chrZ-tmp.snps.focal.vcf.gz /media/queen/extradrive1/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.final.noCV13.chrZ.snps.vcf.gz > ./log/raxml_zfilter.log &
nohup vcftools --gzvcf /media/queen/extradrive1/crotalus_genomic_landscape/vcf/crotalus_genus.noCV13.chrZ-tmp.snps.focal.vcf.gz  --recode --stdout | bgzip -c > /media/queen/extradrive1/crotalus_genomic_landscape/vcf/crotalus_genus.noCV13.chrZ.snps.complete.vcf.gz > ./log/raxml_z.log &

```

### Convert vcf to phylip
Now, we will convert the vcf to a phylip format (by default) and a nexus (-n).

```
cd /media/queen/extradrive1/crotalus_genomic_landscape/analysis/raxml

curl -LO https://raw.githubusercontent.com/edgardomortiz/vcf2phylip/master/vcf2phylip.py
curl -LO https://raw.githubusercontent.com/btmartin721/raxml_ascbias/master/ascbias.py

chmod +x ascbias.py

nohup ./vcf2phylip.py -n -i /media/queen/extradrive1/crotalus_genomic_landscape/vcf/crotalus_genus.noCV13.auto+chrZ.snps.complete.vcf.gz -o CA0346 > ./vcf2phylip_alldata.log & 

nohup ./ascbias.py -p crotalus_genus.noCV13.auto+chrZ.snps.complete.min4.phy -o crotalus_genus.noCV13.auto+chrZ.snps.complete.min4.filtered.phy > ./ascbias_auto+z.log &

mkdir autosome_only

cd autosome_only
nohup ../vcf2phylip.py -n -i /media/queen/extradrive1/crotalus_genomic_landscape/vcf/crotalus_genus.noCV13.auto.snps.complete.vcf.gz -o CA0346 > ./vcf2phylip_auto.log & 

nohup ../ascbias.py -p crotalus_genus.noCV13.auto.snps.complete.min4.phy -o crotalus_genus.noCV13.auto.snps.complete.min4.filtered.phy > ./ascbias_auto.log &

mkdir ../z_only 

cd ../z_only 

nohup ../vcf2phylip.py -n -i /media/queen/extradrive1/crotalus_genomic_landscape/vcf/crotalus_genus.noCV13.chrZ.snps.complete.vcf.gz -o CA0346 > ./vcf2phylip_z.log & 

nohup ../ascbias.py -p crotalus_genus.noCV13.chrZ.snps.complete.min4.phy -o crotalus_genus.noCV13.chrZ.snps.complete.min4.filtered.phy > ./ascbias_z.log & 

```


### Run RAxML
We will generate a concatenated tree with RAxML. 

```
conda activate raxml

nohup raxml-ng --threads 24 --all --msa /media/queen/extradrive1/crotalus_genomic_landscape/analysis/raxml/crotalus_genus.noCV13.auto+chrZ.snps.complete.min4.filtered.phy --model GTGTR4+G+ASC_LEWIS --tree pars{10} --bs-trees 100 > ./raxml_auto+z_log.log &

nohup raxml-ng --threads 12 --all --msa /media/queen/extradrive1/crotalus_genomic_landscape/analysis/raxml/autosome_only/crotalus_genus.noCV13.auto.snps.complete.min4.filtered.phy  --model GTGTR4+G+ASC_LEWIS --tree pars{10} --bs-trees 100 > ./raxml_auto_log.log &

nohup raxml-ng --threads 12 --all --msa /media/queen/extradrive1/crotalus_genomic_landscape/analysis/raxml/z_only/crotalus_genus.noCV13.chrZ.snps.complete.min4.filtered.phy  --model GTGTR4+G+ASC_LEWIS --tree pars{10} --bs-trees 100 > ./raxml_z_log.log &


```


### Run Svdquartets

Now we will run svdquartets. We need to download the binary first. We run this on our local machine since we can't get the right libraries on Xenomorph.

```
cd /Users/kfarleigh/Desktop/Projects/cgli/

/Users/kfarleigh/Desktop/Software/paup4a168_osx


exe ./crotalus_genus.noCV13.auto+chrZ.snps.complete.min4.nexus; 
outgroup 1 2; 
set outroot=mono; 
svdq; 
svdq showScores=no seed=19856136 bootstrap nreps=100 treeFile=crotalus_genus.bootstrap.auto+chrZ.tre; 
savetrees file=crotalus_genus.svdq.auto+chrZ.tre savebootp=nodelabels;

```

Do it for autosome only data.

```
cd /Users/kfarleigh/Desktop/Projects/cgli/autosome

cd autosome

/Users/kfarleigh/Desktop/Software/paup4a168_osx
 
exe ./crotalus_genus.noCV13.auto.snps.complete.min4.nexus; 
outgroup 1 2; 
set outroot=mono; 
svdq; 
svdq showScores=no seed=19856136 bootstrap nreps=100 treeFile=crotalus_genus.bootstrap.auto.tre; 
savetrees file=crotalus_genus.svdq.auto.tre savebootp=nodelabels;

```
Do it for Z only. 

```
cd /Users/kfarleigh/Desktop/Github/Crotalus_Genomic_Landscape/analysis/svdquartets

cd z

/Users/kfarleigh/Desktop/Software/paup4a168_osx

exe ./crotalus_genus.noCV13.chrZ.snps.complete.min4.nexus; 
outgroup 1 2; 
set outroot=mono; 
svdq; 
svdq showScores=no seed=19856136 bootstrap nreps=100 treeFile=crotalus_genus.bootstrap.chrZ.tre; 
savetrees file=crotalus_genus.svdq.chrZ.tre savebootp=nodelabels;

```

#### Calculating branch lengths with RAxML

Since we present the svdquartets autosome tree as the tree, we will calculate branch lengths following a tutorial from svdquartets creator (https://phylosolutions.com/tutorials/svdq-qage/svdq-qage-tutorial.html).

```
cd /Users/kfarleigh/Desktop/Projects/cgli/autosome

/Users/kfarleigh/Desktop/Software/paup4a168_osx

exe ./crotalus_genus.noCV13.auto.snps.complete.min4.nexus;

# Read in the tree, this will tell paup which individual is assigned to which species and provide paup our species tree
exe ./crotalus_genus.svdq.auto.wpartitions.tre

# Look at the tree
showtree
```

We can calculate branch lengths with RAxML by specifying the topology we got from SVDquartets.

```
conda activate raxml

cd /media/queen/extradrive1/crotalus_genomic_landscape/analysis/raxml

mkdir branch_lengths

cd ./branch_lengths

nano cgli_svdq_auto.nwk

((CA0346,CR0001),(((((CT10189615,TR0001),TR0002),(TR0003,TR0004)),(((((((SR0002,SR0004),SR0006),SR0005),SR0003),(((SR0007,SR0019),(SR0009,SR0020)),SR0008)),((((SR0010,SR0011),SR0012),(SR0013,SR0015)),SR0014)),((((SR0016,(SR0017,SR0018)),(((((((SR0021,SR0024),SR0078),(SR0022,SR0023)),(SR0030,SR0031)),((SR0025,SR0028),((SR0026,SR0029),SR0027))),SR0086),SR0088)),(((((SR0041,(SR0042,(SR0051,SR0052))),SR0045),(SR0043,(SR0050,SR0056))),(SR0046,(SR0047,SR0048))),((SR0049,((SR0079,SR0081),SR0080)),SR0082))),(((SR0044,SR0054),SR0083),(SR0053,SR0055))))),((((((((((CV0004,CV0011),CV0009),(CV0006,CV0010)),CV0008),((CV0634,CV0636),(((((((CV0853,CV0854),(CV0864,(CV0865,CV0867))),(CV0860,CV0863)),(CV0857,(CV0858,((CV0862,CV0869),CV0868)))),(CV0856,CV0866)),CV0870),CV0859))),((CV0696,CV0945),CV0697)),CV0946),(((((CV0947,CV0948),CV0950),CV1038),CV0951),CV1037)),CV0949),((((((CV0047,CV0053),CV0054),((((CV0085,((((CV0093,CV0111),CV0136),(((CV0094,(CV0095,CV0096)),CV0098),CV0105)),CV0137)),(CV0159,CV0160)),(((((CV0086,CV0155),(CV0151,CV0152)),(CV0087,CV0148)),CV0145),(CV0157,CV0158))),(((CV0135,CV0162),CV0161),CV0153))),(CV0083,(CV0101,(((((((CV0764,CV0770),((CV0766,CV0800),(CV0772,CV0787))),(CV0781,CV0793)),(CV0775,(CV0784,CV0798))),(CV0783,(CV0790,CV0796))),CV0786),CV0780)))),(((CV0089,CV0150),CV0092),((CV0141,CV0147),CV0144))),((((((CV0808,CV0809),CV0810),(((CV0816,CV0819),CV0818),CV0817)),(((CV0811,CV0813),(CV0814,CV0815)),CV0812)),(((CV0979,(CV0984,CV0988)),(CV0981,CV0982)),(CV0980,CV0983))),(((((CV0985,CV0989),CV0986),CV1149),(CV0990,CV1148)),(CV0991,(((CV1150,CV1151),((CV1153,CV1227),(CV1154,CV1229))),CV1152))))))));

# Following slide 71 of https://cme.h-its.org/exelixis/pubs/ISCBacademy2020_RAxML-NG.pdf
nohup raxml-ng --evaluate --threads 24 --tree cgli_svdq_auto.nwk --msa /media/queen/extradrive1/crotalus_genomic_landscape/analysis/raxml/autosome_only/autosome_complete/crotalus_genus.noCV13.auto.snps.complete.min4.filtered.phy --outgroup CA0346,CR0001 --model GTGTR4+G+ASC_LEWIS --prefix auto_branchlens --redo > raxml_auto_bl_log.log & 
```

### Run treesplits4

```
# First, we will filter to match the raxml input
cd /media/queen/extradrive1/crotalus_genomic_landscape/analysis/fD

conda activate pixy

nohup python3 ./genomics_general/filterGenotypes.py -i crotalus_genus_auto+chrZ.geno.gz -o crotalus_genus_auto+chrZ.10kbthin.geno.gz --thinDist 10000 > splitstree_filt.log &

nohup python3 ./genomics_general/distMat.py -g crotalus_genus_auto+chrZ.10kbthin.geno.gz -T 1 -f phased --windType cat -o ./distmat/croatlus_genus_output.dist > distmat.log &


```

### Process and visualize results

We use the `cgli_raxml+svdquartets_plotting.R` script to plot the trees in R. The script is also available on [GitHub](https://github.com/kfarleigh/Crotalus_Genomic_Landscape)


### Script

`cgli_raxml+svdquartets_plotting.R`

```
# Purpose: Plot RAxML and SVDquartets results for the crotalus genomic landscape of introgression project
# Author: Keaka Farleigh
# Email: keakafarleigh@virginia.edu
# Date: April 4th, 2025

### Load you packages, set directory
library(phytools)
library(ape)

setwd("/Users/kfarleigh/Desktop/Github/Crotalus_Genomic_Landscape/analysis/raxml")

### Read in metadata
mdat <- read.csv("../../Metadata/Crotalus_genomic_landscape_final_metadata.csv", header = T)

mdat <- mdat[,1:2]

colnames(mdat) <- c("Sample", "Species")

mdat$Species <- gsub("Crotalus ", "", mdat$Species)

mdat$color <- NA

mdat$color[which(mdat$Species == "atrox" | mdat$Species == "ruber")] <- "black"
mdat$color[which(mdat$Species == "tigris")] <- "#e4aa24"
mdat$color[which(mdat$Species == "stephensi")] <- "#f0c874"
mdat$color[which(mdat$Species == "pyrrhus")] <- "#72eec4"
mdat$color[which(mdat$Species == "mitchellii")] <- "#e03330"
mdat$color[which(mdat$Species == "angelensis")] <- "#FFB6C1"
mdat$color[which(mdat$Species == "concolor")] <- "#c17889"
mdat$color[which(mdat$Species == "helleri")] <- "#5973b9"
mdat$color[which(mdat$Species == "lutosus")] <- "#8ecddd"
mdat$color[which(mdat$Species == "oreganus")] <- "#777e98"
mdat$color[which(mdat$Species == "viridis")] <- "#78c184"

# Set outgroup samples
outgroup <- c('CA0346', "CR0001")

### Read in tree with bootstrap supports

support_tree <- read.tree('./crotalus_genus.auto+chrZ.snps.10kbthin.filtered.phy.raxml.support')

plot(support_tree)

support_tree_rooted <- root(support_tree, outgroup = outgroup, resolve.root = TRUE)

support_tree_rooted <- ladderize(support_tree_rooted)

plot(support_tree_rooted, cex = 0.7, no.margin = T,edge.width=1)

nodelabels(support_tree_rooted$node.label,frame= "none", cex = 0.5, adj = -0.25)

# Get order of individuals
plot_ord <- support_tree_rooted$tip.label


########################################
##### Plot with better node labels #####
########################################
# Order to match tree order
mdat_ord <- mdat[match(plot_ord, mdat$Sample),]

# Seperate Bootstrap and alrt values 
Nod_vals <- as.data.frame(support_tree_rooted$node.label)
Nod_vals$Bootstrap <- as.numeric(gsub('.*/', '', Nod_vals$`support_tree_rooted$node.label`))


# Set background colors based on bootstrapping thresholds 
Nod_vals$Color[Nod_vals$Bootstrap > 95] <- 'Black'
Nod_vals$Color[Nod_vals$Bootstrap <= 95] <- 'Grey'
Nod_vals$Color[Nod_vals$Bootstrap <= 90] <- 'White'
Nod_vals$Color[Nod_vals$Bootstrap <= 85] <- 'NA'


# Set the outline color so it doesn't plot blank circles 
Nod_vals$outline <- Nod_vals$Color
Nod_vals$outline[Nod_vals$outline == 'NA'] <- 'NA'
Nod_vals$outline[Nod_vals$outline != 'NA'] <- 'Black'


# Plot it
plot(support_tree_rooted, no.margin=TRUE, edge.width=2, cex=0.5, tip.color = mdat_ord$color)
nodelabels(node = 1:support_tree_rooted$Nnode+Ntip(support_tree_rooted),pch = 21, cex=1.25, col= Nod_vals$outline, bg = Nod_vals$Color)

# Get list of species for legend
specs <- unique(mdat_ord[,2:3])

# Add legends
legend("bottomleft",paste("Crotalus",specs$Species, sep = " "), pch=22,pt.bg= specs$color,cex=0.9,
pt.cex=2,title="Species",bty="n")

legend("left", c("> 95", "<= 95", "<= 90 & > 85"), pch = 21, pt.bg = c("black", "grey", "white"), 
       cex = 0.9, pt.cex = 2, title = "Bootstrap support", bty="n")


#########################
##### Autosome only #####
#########################

### Read in tree with bootstrap supports

support_tree <- read.tree('./autosome/crotalus_genus.auto.snps.10kbthin.filtered.phy.raxml.support')

plot(support_tree)

support_tree_rooted <- root(support_tree, outgroup = outgroup, resolve.root = TRUE)

support_tree_rooted <- ladderize(support_tree_rooted)

plot(support_tree_rooted, cex = 0.7, no.margin = T,edge.width=1)

nodelabels(support_tree_rooted$node.label,frame= "none", cex = 0.5, adj = -0.25)

# Get order of individuals
plot_ord <- support_tree_rooted$tip.label


########################################
##### Plot with better node labels #####
########################################
# Order to match tree order
mdat_ord <- mdat[match(plot_ord, mdat$Sample),]

# Seperate Bootstrap and alrt values 
Nod_vals <- as.data.frame(support_tree_rooted$node.label)
Nod_vals$Bootstrap <- as.numeric(gsub('.*/', '', Nod_vals$`support_tree_rooted$node.label`))


# Set background colors based on bootstrapping thresholds 
Nod_vals$Color[Nod_vals$Bootstrap > 95] <- 'Black'
Nod_vals$Color[Nod_vals$Bootstrap <= 95] <- 'Grey'
Nod_vals$Color[Nod_vals$Bootstrap <= 90] <- 'White'
Nod_vals$Color[Nod_vals$Bootstrap <= 85] <- 'NA'


# Set the outline color so it doesn't plot blank circles 
Nod_vals$outline <- Nod_vals$Color
Nod_vals$outline[Nod_vals$outline == 'NA'] <- 'NA'
Nod_vals$outline[Nod_vals$outline != 'NA'] <- 'Black'


# Plot it
plot(support_tree_rooted, no.margin=TRUE, edge.width=2, cex=0.5, tip.color = mdat_ord$color)
nodelabels(node = 1:support_tree_rooted$Nnode+Ntip(support_tree_rooted),pch = 21, cex=1.25, col= Nod_vals$outline, bg = Nod_vals$Color)

# Get list of species for legend
specs <- unique(mdat_ord[,2:3])

# Add legends
legend("bottomleft",paste("Crotalus",specs$Species, sep = " "), pch=22,pt.bg= specs$color,cex=0.9,
       pt.cex=2,title="Species",bty="n")

legend("left", c("> 95", "<= 95", "<= 90 & > 85"), pch = 21, pt.bg = c("black", "grey", "white"), 
       cex = 0.9, pt.cex = 2, title = "Bootstrap support", bty="n")

##################
##### Z only #####
##################

### Read in tree with bootstrap supports

support_tree <- read.tree('./z/crotalus_genus.chrZ.snps.10kbthin.filtered.phy.raxml.support')

plot(support_tree)

support_tree_rooted <- root(support_tree, outgroup = outgroup, resolve.root = TRUE)

support_tree_rooted <- ladderize(support_tree_rooted)

plot(support_tree_rooted, cex = 0.7, no.margin = T,edge.width=1)

nodelabels(support_tree_rooted$node.label,frame= "none", cex = 0.5, adj = -0.25)

# Get order of individuals
plot_ord <- support_tree_rooted$tip.label


########################################
##### Plot with better node labels #####
########################################
# Order to match tree order
mdat_ord <- mdat[match(plot_ord, mdat$Sample),]

# Seperate Bootstrap and alrt values 
Nod_vals <- as.data.frame(support_tree_rooted$node.label)
Nod_vals$Bootstrap <- as.numeric(gsub('.*/', '', Nod_vals$`support_tree_rooted$node.label`))


# Set background colors based on bootstrapping thresholds 
Nod_vals$Color[Nod_vals$Bootstrap > 95] <- 'Black'
Nod_vals$Color[Nod_vals$Bootstrap <= 95] <- 'Grey'
Nod_vals$Color[Nod_vals$Bootstrap <= 90] <- 'White'
Nod_vals$Color[Nod_vals$Bootstrap <= 85] <- 'NA'


# Set the outline color so it doesn't plot blank circles 
Nod_vals$outline <- Nod_vals$Color
Nod_vals$outline[Nod_vals$outline == 'NA'] <- 'NA'
Nod_vals$outline[Nod_vals$outline != 'NA'] <- 'Black'


# Plot it
plot(support_tree_rooted, no.margin=TRUE, edge.width=2, cex=0.5, tip.color = mdat_ord$color)
nodelabels(node = 1:support_tree_rooted$Nnode+Ntip(support_tree_rooted),pch = 21, cex=1.25, col= Nod_vals$outline, bg = Nod_vals$Color)

# Get list of species for legend
specs <- unique(mdat_ord[,2:3])

# Add legends
legend("bottomleft",paste("Crotalus",specs$Species, sep = " "), pch=22,pt.bg= specs$color,cex=0.9,
       pt.cex=2,title="Species",bty="n")

legend("left", c("> 95", "<= 95", "<= 90 & > 85"), pch = 21, pt.bg = c("black", "grey", "white"), 
       cex = 0.9, pt.cex = 2, title = "Bootstrap support", bty="n")


# !!!!!!!!!!!!!!!!!!! #
####################### 
##### SVDquartets #####
#######################
# !!!!!!!!!!!!!!!!!!! #

########################
##### Full dataset #####
########################

### Read in tree with bootstrap supports

support_tree <- read.nexus('/Users/kfarleigh/Desktop/Github/Crotalus_Genomic_Landscape/analysis/svdquartets/crotalus_genus.svdq.auto+chrZ.tre')

plot(support_tree)

support_tree_rooted <- root(support_tree, outgroup = outgroup, resolve.root = TRUE)

support_tree_rooted <- ladderize(support_tree_rooted)

plot(support_tree_rooted, cex = 0.7, no.margin = T,edge.width=1)

nodelabels(support_tree_rooted$node.label,frame= "none", cex = 0.5, adj = -0.25)

# Get order of individuals
plot_ord <- support_tree_rooted$tip.label


########################################
##### Plot with better node labels #####
########################################
# Order to match tree order
mdat_ord <- mdat[match(plot_ord, mdat$Sample),]

# Seperate Bootstrap and alrt values 
Nod_vals <- as.data.frame(support_tree_rooted$node.label)
Nod_vals$Bootstrap <- as.numeric(gsub('.*/', '', Nod_vals$`support_tree_rooted$node.label`))


# Set background colors based on bootstrapping thresholds 
Nod_vals$Color[Nod_vals$Bootstrap > 95] <- 'Black'
Nod_vals$Color[Nod_vals$Bootstrap <= 95] <- 'Grey'
Nod_vals$Color[Nod_vals$Bootstrap <= 90] <- 'White'
Nod_vals$Color[Nod_vals$Bootstrap <= 85] <- 'NA'


# Set the outline color so it doesn't plot blank circles 
Nod_vals$outline <- Nod_vals$Color
Nod_vals$outline[Nod_vals$outline == 'NA'] <- 'NA'
Nod_vals$outline[Nod_vals$outline != 'NA'] <- 'Black'


# Plot it
plot(support_tree_rooted, no.margin=TRUE, edge.width=2, cex=0.5, tip.color = mdat_ord$color)
nodelabels(node = 1:support_tree_rooted$Nnode+Ntip(support_tree_rooted),pch = 21, cex=1.25, col= Nod_vals$outline, bg = Nod_vals$Color)

# Get list of species for legend
specs <- unique(mdat_ord[,2:3])

# Add legends
legend("bottomleft",paste("Crotalus",specs$Species, sep = " "), pch=22,pt.bg= specs$color,cex=0.9,
       pt.cex=2,title="Species",bty="n")

legend("left", c("> 95", "<= 95 & > 90", "<= 90 & > 85"), pch = 21, pt.bg = c("black", "grey", "white"), 
       cex = 0.9, pt.cex = 2, title = "Bootstrap support", bty="n")

title("SVDquartets All Data")
########################
##### Autosome only ####
########################

### Read in tree with bootstrap supports

support_tree <- read.nexus('/Users/kfarleigh/Desktop/Github/Crotalus_Genomic_Landscape/analysis/svdquartets/autosome/crotalus_genus.svdq.auto.tre')

plot(support_tree)

support_tree_rooted <- root(support_tree, outgroup = outgroup, resolve.root = TRUE)

support_tree_rooted <- ladderize(support_tree_rooted)

plot(support_tree_rooted, cex = 0.7, no.margin = T,edge.width=1)

nodelabels(support_tree_rooted$node.label,frame= "none", cex = 0.5, adj = -0.25)

# Get order of individuals
plot_ord <- support_tree_rooted$tip.label


########################################
##### Plot with better node labels #####
########################################
# Order to match tree order
mdat_ord <- mdat[match(plot_ord, mdat$Sample),]

# Seperate Bootstrap and alrt values 
Nod_vals <- as.data.frame(support_tree_rooted$node.label)
Nod_vals$Bootstrap <- as.numeric(gsub('.*/', '', Nod_vals$`support_tree_rooted$node.label`))


# Set background colors based on bootstrapping thresholds 
Nod_vals$Color[Nod_vals$Bootstrap > 95] <- 'Black'
Nod_vals$Color[Nod_vals$Bootstrap <= 95] <- 'Grey'
Nod_vals$Color[Nod_vals$Bootstrap <= 90] <- 'White'
Nod_vals$Color[Nod_vals$Bootstrap <= 85] <- 'NA'


# Set the outline color so it doesn't plot blank circles 
Nod_vals$outline <- Nod_vals$Color
Nod_vals$outline[Nod_vals$outline == 'NA'] <- 'NA'
Nod_vals$outline[Nod_vals$outline != 'NA'] <- 'Black'


# Plot it
plot(support_tree_rooted, no.margin=TRUE, edge.width=2, cex=0.5, tip.color = mdat_ord$color)
nodelabels(node = 1:support_tree_rooted$Nnode+Ntip(support_tree_rooted),pch = 21, cex=1.25, col= Nod_vals$outline, bg = Nod_vals$Color)

# Get list of species for legend
specs <- unique(mdat_ord[,2:3])

# Add legends
legend("bottomleft",paste("Crotalus",specs$Species, sep = " "), pch=22,pt.bg= specs$color,cex=0.9,
       pt.cex=2,title="Species",bty="n")

legend("left", c("> 95", "<= 95 & > 90", "<= 90 & > 85"), pch = 21, pt.bg = c("black", "grey", "white"), 
       cex = 0.9, pt.cex = 2, title = "Bootstrap support", bty="n")

title("SVDquartets Autosome")

################
##### Zonly ####
################

### Read in tree with bootstrap supports

support_tree <- read.nexus('../svdquartets/z/crotalus_genus.svdq.chrZ.tre')

plot(support_tree)

support_tree_rooted <- root(support_tree, outgroup = outgroup, resolve.root = TRUE)

support_tree_rooted <- ladderize(support_tree_rooted)

plot(support_tree_rooted, cex = 0.7, no.margin = T,edge.width=1)

nodelabels(support_tree_rooted$node.label,frame= "none", cex = 0.5, adj = -0.25)

plot(support_tree_rooted, cex = 0.7, no.margin = T,edge.width=1)

nodelabels(text=1:support_tree_rooted$Nnode,node=1:support_tree_rooted$Nnode+Ntip(support_tree_rooted)) 

support_tree_rooted <- rotateNodes(support_tree_rooted, nodes = 142) # Still testing

# Get order of individuals
plot_ord <- support_tree_rooted$tip.label


########################################
##### Plot with better node labels #####
########################################
# Order to match tree order
mdat_ord <- mdat[match(plot_ord, mdat$Sample),]

# Seperate Bootstrap and alrt values 
Nod_vals <- as.data.frame(support_tree_rooted$node.label)
Nod_vals$Bootstrap <- as.numeric(gsub('.*/', '', Nod_vals$`support_tree_rooted$node.label`))


# Set background colors based on bootstrapping thresholds 
Nod_vals$Color[Nod_vals$Bootstrap > 95] <- 'Black'
Nod_vals$Color[Nod_vals$Bootstrap <= 95] <- 'Grey'
Nod_vals$Color[Nod_vals$Bootstrap <= 90] <- 'White'
Nod_vals$Color[Nod_vals$Bootstrap <= 85] <- 'NA'


# Set the outline color so it doesn't plot blank circles 
Nod_vals$outline <- Nod_vals$Color
Nod_vals$outline[Nod_vals$outline == 'NA'] <- 'NA'
Nod_vals$outline[Nod_vals$outline != 'NA'] <- 'Black'


# Plot it
plot(support_tree_rooted, no.margin=TRUE, edge.width=2, cex=0.5, tip.color = mdat_ord$color)
nodelabels(node = 1:support_tree_rooted$Nnode+Ntip(support_tree_rooted),pch = 21, cex=1.25, col= Nod_vals$outline, bg = Nod_vals$Color)

# Get list of species for legend
specs <- unique(mdat_ord[,2:3])

# Add legends
legend("bottomleft",paste("Crotalus",specs$Species, sep = " "), pch=22,pt.bg= specs$color,cex=0.9,
       pt.cex=2,title="Species",bty="n")

legend("left", c("> 95", "<= 95 & > 90", "<= 90 & > 85"), pch = 21, pt.bg = c("black", "grey", "white"), 
       cex = 0.9, pt.cex = 2, title = "Bootstrap support", bty="n")

title("SVDquartets Z")

```

# Crotalus genomic landscape, ADMIXTURE & PCA
Author: Keaka Farleigh, Ph.D.
Date: April 18th, 2025
Email: keakafarleigh@virginia.edu

## Purpose
This script will filter vcfs and estimate ancestry using admixture and principal component analysis. 

## Location 

This analysis was performed in the Schield lab on Xenomorph at the University of Virginia. Xenomorph is a system76 workstation running Ubuntu v22.04.3.


### Filter vcfs
Following [Schield et al. (2024)](https://www.science.org/doi/full/10.1126/science.adj8766), we will only include snps with a minor allele frequency greater than 0.1 (maf > 0.1).

```

cd /media/queen/extradrive1/crotalus_genomic_landscape/analysis

mkdir admixture
mkdir pca

cd admixture 

mkdir auto_only
mkdir z_only

cd ../pca

mkdir auto_only
mkdir z_only

# All data
nohup bcftools view --threads 8 -O z -s ^CA0346,CR0001 -i 'MAF > 0.1' -o /media/queen/extradrive1/crotalus_genomic_landscape/vcf/crotalus_genus.noCV13.auto+chrZ.snps.maf01.ingroup.vcf.gz /media/queen/extradrive1/crotalus_genomic_landscape/vcf/crotalus_genus.noCV13.auto+chrZ-tmp.snps.focal.vcf.gz > ../log/admixture_filt_alldata.log &

# Autosome only
nohup bcftools view --threads 8 -O z -s ^CA0346,CR0001 -i 'MAF > 0.1' -o /media/queen/extradrive1/crotalus_genomic_landscape/vcf/crotalus_genus.noCV13.auto.snps.maf01.ingroup.vcf.gz /media/queen/extradrive1/crotalus_genomic_landscape/vcf/crotalus_genus.noCV13.auto-tmp.snps.focal.vcf.gz > ../log/admixture_filt_auto.log &

# Z only
nohup bcftools view --threads 8 -O z -s ^CA0346,CR0001 -i 'MAF > 0.1' -o /media/queen/extradrive1/crotalus_genomic_landscape/vcf/crotalus_genus.noCV13.chrZ.snps.maf01.ingroup.vcf.gz /media/queen/extradrive1/crotalus_genomic_landscape/vcf/crotalus_genus.noCV13.chrZ-tmp.snps.focal.vcf.gz > ../log/admixture_filt_Z.log &

```

### Install plink

```
conda activate raxml

mamba install plink

```

### Convert vcf to plink format for admixture

```
cd /media/queen/extradrive1/crotalus_genomic_landscape/analysis/admixture

plink --vcf /media/queen/extradrive1/crotalus_genomic_landscape/vcf/crotalus_genus.noCV13.auto+chrZ.snps.maf01.ingroup.vcf.gz --make-bed --out crotalus_genus.noCV13.auto+chrZ.snps.maf01.ingroup --allow-extra-chr

# We need to modify the chromosome names so admixture will run

awk '{$1="0";print $0}' crotalus_genus.noCV13.auto+chrZ.snps.maf01.ingroup.bim > crotalus_genus.noCV13.auto+chrZ.snps.maf01.ingroup.bim.tmp

mv crotalus_genus.noCV13.auto+chrZ.snps.maf01.ingroup.bim.tmp crotalus_genus.noCV13.auto+chrZ.snps.maf01.ingroup.bim

cd auto_only

plink --vcf /media/queen/extradrive1/crotalus_genomic_landscape/vcf/crotalus_genus.noCV13.auto.snps.maf01.ingroup.vcf.gz --make-bed --out crotalus_genus.noCV13.auto.snps.maf01.ingroup --allow-extra-chr

awk '{$1="0";print $0}' crotalus_genus.noCV13.auto.snps.maf01.ingroup.bim > crotalus_genus.noCV13.auto.snps.maf01.ingroup.bim.tmp

mv crotalus_genus.noCV13.auto.snps.maf01.ingroup.bim.tmp  crotalus_genus.noCV13.auto.snps.maf01.ingroup.bim

cd ../z_only

plink --vcf /media/queen/extradrive1/crotalus_genomic_landscape/vcf/crotalus_genus.noCV13.chrZ.snps.maf01.ingroup.vcf.gz --make-bed --out crotalus_genus.noCV13.chrZ.snps.maf01.ingroup --allow-extra-chr

awk '{$1="0";print $0}' crotalus_genus.noCV13.chrZ.snps.maf01.ingroup.bim > crotalus_genus.noCV13.chrZ.snps.maf01.ingroup.bim.tmp

mv crotalus_genus.noCV13.chrZ.snps.maf01.ingroup.bim.tmp  crotalus_genus.noCV13.chrZ.snps.maf01.ingroup.bim

cd ..
```

## Admixture

First, we will find our K, then we will run admixture for that K.

```
# All data 
cd /media/queen/extradrive1/crotalus_genomic_landscape/analysis/admixture

nohup bash -c  'for i in {1..15}; do admixture -j4 --cv crotalus_genus.noCV13.auto+chrZ.snps.maf01.ingroup.bed $i > log${i}.out; done' > ./alldata_cvtest.log &

awk '/CV/ {print $3,$4}' *out | cut -c 4,5,7-20 | sed 's/)//g' | sed 's/://g' > alldata.cv.error

# Autosome only

cd ./auto_only

nohup bash -c  'for i in {1..15}; do admixture -j4 --cv crotalus_genus.noCV13.auto.snps.maf01.ingroup.bed $i > log${i}.out; done' > ./auto_cvtest.log &

awk '/CV/ {print $3,$4}' *out | cut -c 4,5,7-20 | sed 's/)//g' | sed 's/://g' > auto.cv.error

cd ../z_only

nohup bash -c  'for i in {1..15}; do admixture -j4 --cv crotalus_genus.noCV13.chrZ.snps.maf01.ingroup.bed $i > log${i}.out; done' > ./chrZ_cvtest.log &

awk '/CV/ {print $3,$4}' *out | cut -c 4,5,7-20 | sed 's/)//g' | sed 's/://g'> chrZ.cv.error
```

Let's plot the cross-validation results to choose our optimal K value. We will use the `admixture_plotting.R` script. *Note that this script also plots the admixture ancestry plots*.

The cross-validation error suggests a K of 10, but 11-14 are pretty close, so we will plot all of them. We will use the `.Q` files from ADMIXTURE to plot the ancestry proportions and the `.nosex` file to get the order of the individuals in those `.Q` files.


## Principal component analysis

We will use the R pacakge SNPrelate to conduct PCA following the pca.R script. 


``` 

cd /media/queen/extradrive1/crotalus_genomic_landscape/analysis/pca

```


### Scripts

The following R scripts are listed below.

	- `admixture_plotting.R`
	- `pca.R`

`admixture_plotting.R`

```

# Purpose: Plot cross-validation error values from ADMIXTURE.
# Author: Keaka Farleigh
# Date: April 21st, 2025
# Email: keakafarleigh@virginia.edu

### Load your packages, set working directory

library(tidyverse)
library(cowplot)
library(PopGenHelpR)
library(phytools)
library(ape)

setwd("/Users/kfarleigh/Desktop/Github/Crotalus_Genomic_Landscape/analysis/admixture")


##################################
##### Cross-validation error #####
##################################

### Load your data
# We will make three plots since we ran ADMIXTURE on the whole dataset, autosomes only, and Z only. 

alldata <- read.table('alldata.cv.error')
auto <- read.table('./auto_only/auto.cv.error')
chrZ <- read.table('./z_only/chrZ.cv.error')


### Format data
# We need to make the K values (column 1) numeric so they plot in the right order.

alldata$V1 <- as.numeric(alldata$V1)
auto$V1 <- as.numeric(auto$V1)
chrZ$V1 <- as.numeric(chrZ$V1)


### Plot the data 
alldata_plot <- ggplot(data = alldata, aes(x = V1, y = V2)) + geom_point(color = "navy") + 
  scale_x_continuous(name ="K", breaks=c(min(alldata$V1):max(alldata$V1))) + 
  ylab("Cross-validation error") + theme_classic() + ggtitle("All data") + theme(plot.title = element_text(hjust = 0.5))

auto_plot <- ggplot(data = auto, aes(x = V1, y = V2)) + geom_point(color = "navy") + 
  scale_x_continuous(name ="K", breaks=c(min(auto$V1):max(auto$V1))) + 
  ylab("Cross-validation error") + theme_classic() + ggtitle("Autosomes") + theme(plot.title = element_text(hjust = 0.5))

z_plot <- ggplot(data = chrZ, aes(x = V1, y = V2)) + geom_point(color = "navy") + 
  scale_x_continuous(name ="K", breaks=c(min(chrZ$V1):max(chrZ$V1))) + 
  ylab("Cross-validation error") + theme_classic() + ggtitle("Z") + theme(plot.title = element_text(hjust = 0.5))


plot_grid(alldata_plot, auto_plot, z_plot, nrow = 1)



##########################
##### Ancestry plots #####
##########################
# We will create stacked plots for each of the datasets (all, auto, Z)

# Read in .nosex file to get the order of individuals
ind_ord <- read.table('./crotalus_genus.noCV13.auto+chrZ.snps.maf01.ingroup.nosex')

### All data
all_K10 <- read.table("./crotalus_genus.noCV13.auto+chrZ.snps.maf01.ingroup.10.Q")
all_K11 <- read.table("./crotalus_genus.noCV13.auto+chrZ.snps.maf01.ingroup.11.Q")
all_K12 <- read.table("./crotalus_genus.noCV13.auto+chrZ.snps.maf01.ingroup.12.Q")
all_K13 <- read.table("./crotalus_genus.noCV13.auto+chrZ.snps.maf01.ingroup.13.Q")
all_K14 <- read.table("./crotalus_genus.noCV13.auto+chrZ.snps.maf01.ingroup.14.Q")

# Get function to create ancestry plots, this is incorporated as the Ancestry_barchart function in PopGenHelpR v1.5
source("/Users/kfarleigh/Desktop/Projects/PopGenHelpR_offline/Ancestry_barchart_v2.R")

# If on PC
source("/Users/keaka/OneDrive/Desktop/Projects/PopGenHelpR_offline/Ancestry_barchart_v2.R")
source("/Users/keaka/OneDrive/Desktop/Projects/PopGenHelpR_offline/Piechart_map_v2.R")
# Modify the ind_ord object for plotting

colnames(ind_ord)[1:2] <- c("Sample", "Population")

# Replace the "Population" column with species information from our lab metadata sheet
mdat <- read.csv("/Users/kfarleigh/Desktop/Github/Crotalus_Genomic_Landscape/Metadata/Crotalus_genomic_landscape_final_metadata.csv")

# If on PC
mdat <- read.csv("../../Metadata/Crotalus_genomic_landscape_final_metadata.csv")

colnames(mdat)[1] <- "Sample"

ind_ord_join <- left_join(ind_ord, mdat)

ind_ord_join$Population <- ind_ord_join$species

ind_ord <- ind_ord_join[,c(1,2,4,5)]

### Bind individual names to Q matrix

# K of 10
all_K10_wind <- cbind(ind_ord$Sample, all_K10)

colnames(all_K10_wind)[1] <- "Ind"

colnames(ind_ord)[3:4] <- c("Long", "Lat")

# K of 11
all_K11_wind <- cbind(ind_ord$Sample, all_K11)

colnames(all_K11_wind)[1] <- "Ind"

# K of 12
all_K12_wind <- cbind(ind_ord$Sample, all_K12)

colnames(all_K12_wind)[1] <- "Ind"

# K of 13
all_K13_wind <- cbind(ind_ord$Sample, all_K13)

colnames(all_K13_wind)[1] <- "Ind"

# K of 14
all_K14_wind <- cbind(ind_ord$Sample, all_K14)

colnames(all_K14_wind)[1] <- "Ind"

# Create a plot order, we will make one so that it places the specks next to each other and the westerns with each other

specks <- ind_ord_join[which(ind_ord_join$Population == "Crotalus tigris" | ind_ord_join$Population == "Crotalus pyrrhus" | 
                               ind_ord_join$Population == "Crotalus mitchellii" | ind_ord_join$Population == "Crotalus stephensi" |
                               ind_ord_join$Population == "Crotalus angelensis"),]

west <- ind_ord_join[which(ind_ord_join$Population == "Crotalus concolor" | ind_ord_join$Population == "Crotalus helleri" | 
                             ind_ord_join$Population == "Crotalus lutosus" | ind_ord_join$Population == "Crotalus oreganus" |
                             ind_ord_join$Population == "Crotalus viridis"),]


specks_ord <- specks[order(specks$Population),1]
west_ord <- west[order(west$Population),1]

plot_ord <- c(specks_ord, west_ord)

# List colors for reference
spec_col <- c('#c17889','#5973b9','#8ecddd','#e03330','#777e98','#72eec4','#f0c874','#e4aa24','#78c184', "black")

# Viridis was split into two pops, angelensis was with mitchellii/pyrrhus 
col_k10 <- c('#78c184', '#777e98','#72eec4', '#e03330', '#e4aa24', '#48774f', '#8ecddd', '#f0c874', '#5973b9', '#c17889')

# Helleri is now split into two as well, angelensis with mitchelli/pyrrhus
col_k11 <- c('#5973b9','#48774f', '#72eec4', '#8ecddd','#3e5081', '#777e98','#78c184','#e03330','#c17889', '#e4aa24','#f0c874')

# Splits pyrrhus in two
col_k12 <- c('#e03330','#3e5081', '#72eec4','#c17889','#4fa689', '#777e98','#f0c874',  '#e4aa24','#48774f',  '#8ecddd','#78c184','#5973b9')

# Splits viridis into a 3rd pop
col_k13 <- c('#72eec4','#f0c874','#c17889','#48774f','#777e98','#2b472f','#5973b9','#78c184', '#8ecddd','#3e5081','#e03330', '#e4aa24','#4fa689')


# Splits helleri into a 3rd pop, stephensi into 2, concolor into 2
col_k14 <- c('#2c395c','#48774f','#e03330','#777e98', '#8ecddd', '#5973b9','#c17889','#78c184','#72eec4','#e4aa24','#c0a05c','#87545f','#f0c874','#2b472f')

### Plot them
all_K10_plot <- Ancestry_barchart2(all_K10_wind, pops = ind_ord, K = 10, col = col_k10, ind.order = plot_ord) 

K10_multipanel <- all_K10_plot$`Individual Ancestry Plot` + ggtitle("K = 10") +
  theme(axis.text.x  = element_blank(), legend.position = "none", plot.title = element_text(hjust = 0.5))

all_K11_plot <- Ancestry_barchart2(all_K11_wind, pops = ind_ord, K = 11, col = col_k11, ind.order = plot_ord) 

K11_multipanel <- all_K11_plot$`Individual Ancestry Plot` + ggtitle("K = 11") +
  theme(axis.text.x  = element_blank(), legend.position = "none", plot.title = element_text(hjust = 0.5))


all_K12_plot <- Ancestry_barchart2(all_K12_wind, pops = ind_ord, K = 12, col = col_k12, ind.order = plot_ord) 

K12_multipanel <- all_K12_plot$`Individual Ancestry Plot` + ggtitle("K = 12") +
  theme(axis.text.x  = element_blank(), legend.position = "none", plot.title = element_text(hjust = 0.5))

all_K13_plot <- Ancestry_barchart2(all_K13_wind, pops = ind_ord, K = 13, col = col_k13, ind.order = plot_ord) 

K13_multipanel <- all_K13_plot$`Individual Ancestry Plot` + ggtitle("K = 13") +
  theme(axis.text.x  = element_blank(), legend.position = "none", plot.title = element_text(hjust = 0.5))

all_K14_plot <- Ancestry_barchart2(all_K14_wind, pops = ind_ord, K = 14, col = col_k14, ind.order = plot_ord) 

K14_multipanel <- all_K14_plot$`Individual Ancestry Plot` + ggtitle("K = 14") +
  theme(axis.text.x  = element_blank(), legend.position = "none", plot.title = element_text(hjust = 0.5))

# Plot everything together
plot_grid(K10_multipanel, K11_multipanel, 
          K12_multipanel, K13_multipanel, K14_multipanel,
          ncol = 1)


#####################
##### Autosomes #####
#####################


#############
##### Z #####
#############

# Read in .nosex file to get the order of individuals
ind_ord_z <- read.table('./z_only/crotalus_genus.noCV13.chrZ.snps.maf01.ingroup.nosex')

# Check and see if the ind_ord_z and ind_ord are in the same order, should all be TRUE
table(ind_ord_z$V1 == ind_ord$Sample)

### All data
z_K11 <- read.table("./z_only/crotalus_genus.noCV13.chrZ.snps.maf01.ingroup.11.Q")
z_K14 <- read.table("./z_only/crotalus_genus.noCV13.chrZ.snps.maf01.ingroup.14.Q")

# K of 11
z_K11_wind <- cbind(ind_ord$Sample, z_K11)

colnames(z_K11_wind)[1] <- "Ind"

# K of 14
z_K14_wind <- cbind(ind_ord$Sample, z_K14)

colnames(z_K14_wind)[1] <- "Ind"

# Helleri is now split into two as well, angelensis with mitchelli/pyrrhus
z_col_k11 <- c('#87545f','#48774f', '#72eec4', '#8ecddd','#5973b9', '#777e98','#78c184','#e03330','#c17889', '#e4aa24','#f0c874')

# Splits helleri into a 3rd pop, stephensi into 2, concolor into 2
z_col_k14 <- c('#2c395c','#78c184','#72eec4','#474b5b', '#8ecddd', '#5973b9','#c17889','#48774f','#e03330','#e4aa24','#4fa689','#777e98','#f0c874','#2b472f')

### Plot them
z_K11_plot <- Ancestry_barchart2(z_K11_wind, pops = ind_ord, K = 11, col = z_col_k11, ind.order = plot_ord) 

z_K11_multipanel <- z_K11_plot$`Individual Ancestry Plot` + ggtitle("K = 11, Z") +
  theme(axis.text.x  = element_blank(), legend.position = "none", plot.title = element_text(hjust = 0.5))

z_K14_plot <- Ancestry_barchart2(z_K14_wind, pops = ind_ord, K = 14, col = z_col_k14, ind.order = plot_ord) 

z_K14_multipanel <- z_K14_plot$`Individual Ancestry Plot` + ggtitle("K = 14, Z") +
  theme(axis.text.x  = element_blank(), legend.position = "none", plot.title = element_text(hjust = 0.5))

plot_grid(z_K11_multipanel, z_K14_multipanel, ncol = 1)

### Make piechart maps, can't seem to go above 10

all_K10_piemap <- Piechart_map(all_K10_wind, pops = ind_ord, K = 10, plot.type = "all", Lat_buffer = 3, Long_buffer = 3, 
                               col = col_k10, Latitude_col = 3, Longitude_col = 4)

k = 13

all_K14_piemap <- Piechart_map(all_K14_wind, pops = ind_ord, K = 14, plot.type = "all", Lat_buffer = 3, Long_buffer = 3, 
                               col = col_k14, Longitude_col = 3, Latitude_col = 4)



### Let's set the plotting order so that the individuals are the same order as the tree

support_tree <- read.nexus('/Users/kfarleigh/Desktop/crotalus_genus.svdq.auto.tre')

# Set outgroup samples
outgroup <- c('CA0346', "CR0001")

support_tree_rooted <- root(support_tree, outgroup = outgroup, resolve.root = TRUE)

# Ladderize for plotting
support_tree_rooted <- ladderize(support_tree_rooted)

# Get list of tip labels
phylo_ord <- support_tree_rooted$tip.label

# Determine if the edge is associated with a tip 
is_tip <- support_tree_rooted$edge[,2] <= length(support_tree_rooted$tip.label)

# Get an order of the plotting (bottom to top)
ordered_tips <- support_tree_rooted$edge[is_tip, 2]

# Get the order for our structure plots
plot_ord_tree <- phylo_ord[ordered_tips]

# Remove outgroup samples
plot_ord_tree <- plot_ord_tree[1:180]

# Generate our plots again using new tree order 
all_K10_treeplot <- Ancestry_barchart2(all_K10_wind, pops = ind_ord, K = 10, col = col_k10, ind.order = plot_ord_tree) 

K10_multipanel_tree <- all_K10_treeplot$`Individual Ancestry Plot` + ggtitle("K = 10") +
  theme(axis.text.y  = element_blank(), legend.position = "none", plot.title = element_text(hjust = 0.5)) + coord_flip()

all_K11_treeplot<- Ancestry_barchart2(all_K11_wind, pops = ind_ord, K = 11, col = col_k11, ind.order = plot_ord_tree) 

K11_multipanel_tree <- all_K11_treeplot$`Individual Ancestry Plot` + ggtitle("K = 11") +
  theme(axis.text.y  = element_blank(), legend.position = "none", plot.title = element_text(hjust = 0.5)) + coord_flip()


all_K12_treeplot <- Ancestry_barchart2(all_K12_wind, pops = ind_ord, K = 12, col = col_k12, ind.order = plot_ord_tree) 

K12_multipanel_tree <- all_K12_treeplot$`Individual Ancestry Plot` + ggtitle("K = 12") +
  theme(axis.text.y  = element_blank(), legend.position = "none", plot.title = element_text(hjust = 0.5)) + coord_flip()

all_K13_treeplot <- Ancestry_barchart2(all_K13_wind, pops = ind_ord, K = 13, col = col_k13, ind.order = plot_ord_tree) 

K13_multipanel_tree <- all_K13_treeplot$`Individual Ancestry Plot` + ggtitle("K = 13") +
  theme(axis.text.y  = element_blank(), legend.position = "none", plot.title = element_text(hjust = 0.5)) + coord_flip()

all_K14_treeplot <- Ancestry_barchart2(all_K14_wind, pops = ind_ord, K = 14, col = col_k14, ind.order = plot_ord_tree) 

K14_multipanel_tree <- all_K14_treeplot$`Individual Ancestry Plot` + ggtitle("K = 14") +
  theme(axis.text.y  = element_blank(), legend.position = "none", plot.title = element_text(hjust = 0.5)) + coord_flip()


plot_grid(K10_multipanel_tree, K11_multipanel_tree,
          K12_multipanel_tree, K13_multipanel_tree,
          K14_multipanel_tree, ncol = 5)

# Z 

z_K11_treeplot <- Ancestry_barchart2(z_K11_wind, pops = ind_ord, K = 11, col = z_col_k11, ind.order = plot_ord_tree) 

z_K11_multipanel_tree <- z_K11_treeplot$`Individual Ancestry Plot` + ggtitle("K = 11, Z") +
  theme(axis.text.y  = element_blank(), legend.position = "none", plot.title = element_text(hjust = 0.5)) + coord_flip()

z_K14_treeplot <- Ancestry_barchart2(z_K14_wind, pops = ind_ord, K = 14, col = z_col_k14, ind.order = plot_ord_tree) 

z_K14_multipanel_tree <- z_K14_treeplot$`Individual Ancestry Plot` + ggtitle("K = 14, Z") +
  theme(axis.text.y  = element_blank(), legend.position = "none", plot.title = element_text(hjust = 0.5)) + coord_flip()

plot_grid(z_K11_multipanel_tree, z_K14_multipanel_tree, ncol = 2)


```

`pca.R`

```
# Purpose: Perform principal component analysis (PCA) on a snp dataset
# Author: Keaka Farleigh
# Date: April 22nd, 2025
# Email: keakafarleigh@virginia.edu

### Load your package

library(SNPRelate)
library(tidyverse)
library(cowplot)


### Read in vcfs 

all_dat <- snpgdsVCF2GDS("../../vcf/crotalus_genus.noCV13.auto+chrZ.snps.maf01.ingroup.vcf.gz", "cgli_alldata.gds")

auto_dat <- snpgdsVCF2GDS("../../vcf/crotalus_genus.noCV13.auto.snps.maf01.ingroup.vcf.gz", "cgli_autodata.gds")

z_dat <- snpgdsVCF2GDS("../../vcf/crotalus_genus.noCV13.chrZ.snps.maf01.ingroup.vcf.gz", "cgli_zdata.gds")

### Run PCA
# We open the geno file (.gds) file.

# All 
all_dat_geno <- snpgdsOpen("cgli_alldata.gds")

all_dat_pca <- snpgdsPCA(all_dat_geno, num.thread = 8, autosome.only = F)

# Auto 
auto_dat_geno <- snpgdsOpen("cgli_autodata.gds")

auto_dat_pca <- snpgdsPCA(auto_dat_geno, num.thread = 8, autosome.only = F)

# Z
z_dat_geno <- snpgdsOpen("cgli_zdata.gds")

z_dat_pca <- snpgdsPCA(z_dat_geno, num.thread = 8, autosome.only = F)

# Save workspace for plotting on local machine

save.image("Cgli_pca.Rdata")

### Plot on local machine 

setwd("/Users/kfarleigh/Desktop/Github/Crotalus_Genomic_Landscape/analysis/pca/")

load("Cgli_pca.Rdata")

### Read in metadata
mdat <- read.delim("../../Metadata/Crotalus_genomic_landscape_pixypopmap.txt", header = F)

colnames(mdat) <- c("Sample", "Species")

mdat$color <- NA

mdat$color[which(mdat$Species == "atrox" | mdat$Species == "ruber")] <- "black"
mdat$color[which(mdat$Species == "tigris")] <- "#e4aa24"
mdat$color[which(mdat$Species == "stephensi")] <- "#f0c874"
mdat$color[which(mdat$Species == "pyrrhus")] <- "#72eec4"
mdat$color[which(mdat$Species == "mitchellii")] <- "#e03330"
mdat$color[which(mdat$Species == "angelensis")] <- "#FFB6C1"
mdat$color[which(mdat$Species == "concolor")] <- "#c17889"
mdat$color[which(mdat$Species == "helleri")] <- "#5973b9"
mdat$color[which(mdat$Species == "lutosus")] <- "#8ecddd"
mdat$color[which(mdat$Species == "oreganus")] <- "#777e98"
mdat$color[which(mdat$Species == "viridis")] <- "#78c184"

## Create data frames of individuals and their loadings on the first two PCs
# All
all_dat_df <- data.frame(Sample = all_dat_pca$sample.id, PC1 = all_dat_pca$eigenvect[,1], PC2 = all_dat_pca$eigenvect[,2])

all_dat_df_comb <- left_join(all_dat_df, mdat)

all_dat_pcaplot <- ggplot(data = all_dat_df_comb, aes(x = PC1, y = PC2)) + geom_point(aes(color = Species)) + theme_classic() +
  xlab(paste("Principal component 1: %", round((all_dat_pca$varprop[1]*100),2), " variance explained", sep = '')) +
  ylab(paste("Principal component 2: %", round((all_dat_pca$varprop[2]*100),2), " variance explained", sep = '')) +
  scale_color_manual(values = c("#FFB6C1","#c17889","#5973b9","#8ecddd","#e03330","#777e98","#72eec4","#f0c874","#e4aa24", "#78c184")) +
  ggtitle("All") + theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

# Auto
auto_dat_df <- data.frame(Sample = auto_dat_pca$sample.id, PC1 = auto_dat_pca$eigenvect[,1], PC2 = auto_dat_pca$eigenvect[,2])
  
auto_dat_df_comb <- left_join(auto_dat_df, mdat)

auto_dat_pcaplot <- ggplot(data = auto_dat_df_comb, aes(x = PC1, y = PC2)) + geom_point(aes(color = Species)) + theme_classic() +
  xlab(paste("Principal component 1: %", round((auto_dat_pca$varprop[1]*100),2), " variance explained", sep = '')) +
  ylab(paste("Principal component 2: %", round((auto_dat_pca$varprop[2]*100),2), " variance explained", sep = '')) +
  scale_color_manual(values = c("#FFB6C1","#c17889","#5973b9","#8ecddd","#e03330","#777e98","#72eec4","#f0c874","#e4aa24", "#78c184")) +
  ggtitle("Autosomes") + theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

# Z
z_dat_df <- data.frame(Sample = z_dat_pca$sample.id, PC1 = z_dat_pca$eigenvect[,1], PC2 = z_dat_pca$eigenvect[,2])

z_dat_df_comb <- left_join(z_dat_df, mdat)

z_dat_pcaplot <- ggplot(data = z_dat_df_comb, aes(x = PC1, y = PC2)) + geom_point(aes(color = Species)) + theme_classic() +
  xlab(paste("Principal component 1: %", round((z_dat_pca$varprop[1]*100),2), " variance explained", sep = '')) +
  ylab(paste("Principal component 2: %", round((z_dat_pca$varprop[2]*100),2), " variance explained", sep = '')) +
  scale_color_manual(values = c("#FFB6C1","#c17889","#5973b9","#8ecddd","#e03330","#777e98","#72eec4","#f0c874","#e4aa24", "#78c184")) +
  ggtitle("Z") + theme(plot.title = element_text(hjust = 0.5), legend.position = "none")


plot_grid(all_dat_pcaplot, auto_dat_pcaplot, z_dat_pcaplot, nrow = 1)

```

