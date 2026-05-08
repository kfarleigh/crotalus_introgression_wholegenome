# Pipeline for assigning *Crotalus pyrrhus* scaffolds to chromosomes

Date: September 23rd, 2024

Author: Keaka Farleigh 

Email: keakafarleigh@virginia.edu

## Overview

This README file contains details about the pipeline used to anchor *Crotalus pyrrhus* reference
genome scaffolds to chromosomes using the reference genomes of *C. viridis* (stored in house) and *C. adamanteus* (downloaded from NCBI: [GCA_039797435.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_039797435.1/)).

The [Schield lab](https://schieldlab.org/) recently generated genomic scaffolds for *C. pyrrhus*, which need to be assigned to chromosomes for downstream analysis.  

The results of this analysis will be useful for many future projects in the Schield lab, including genomic landscape projects, genome scans, and sex chromosome evolution. 


## Environment and resources

*Computer and software*

- Xenomorph
- [Mashmap](https://github.com/marbl/MashMap)
	- Version 3.1.3
	- Located at /usr/local/bin/mashmap
- [Python](https://www.python.org/downloads/)
	- Version 2.7.16
	
## Data formatting

First, we need to set up our directory structure. **This directory structure should mirror the structure on the Schield lab Google drive**.

``` 
# Navigate to extradrive1
cd /media/queen/extradrive1

# Make genome folder
mkdir ./crotalus_pyrrhus_genome

# Make the chromosome assignment folder
mkdir ./crotalus_pyrrhus_genome/chromosome_assignment

# Navigate into genome folder 
cd ./crotalus_pyrrhus_genome


# For haplotype 2
cd /media/queen/extradrive1/crotalus_pyrrhus_genome

mkdir chromosome_assignment_hap2
```

Now we create a folder that holds the *C. viridis* and *C. adamanteus* reference genomes in the `chromosome_assignment` folder. 

```
mkdir ./chromosome_assignment/reference_genomes

# For haplotype 2
mkdir chromosome_assignment_hap2/reference_genomes

```


Then we download the data into the `reference_genomes` folder. The ncbi data will be a compressed (.zip) directory (unless you use sratoolkit).We need to unzip the *C. adamanteus* data from NCBI and rename the file so that it makes sense to us. 

```
cd ./chromosome_assignment/reference_genomes

unzip ncbi_dataset.zip 

# Navigate to adamanteus data
cd ./ncbi_dataset/data/GCA_039797435.1

# Rename file
mv GCA_039797435.1_Cadamanteus_3dDNAHiC_1.2_genomic.fna C_adamanteus_genome_GCA_039797435.1.fasta


# Make fasta executable
chmod +x  C_adamanteus_genome_GCA_039797435.1.fasta

# Move it up to the reference genomes directory
mv C_adamanteus_genome_GCA_039797435.1.fasta /media/queen/extradrive1/crotalus_pyrrhus_genome/chromosome_assignment/reference_genomes


# Navigate to the reference genomes, remove extra nbci data

cd /media/queen/extradrive1/crotalus_pyrrhus_genome/chromosome_assignment/reference_genomes

rm ncbi* -r
rm md5sum.txt
rm README.md

# For haplotype 2, copy over old files 

cp /media/queen/extradrive1/crotalus_pyrrhus_genome/chromosome_assignment_hap1/reference_genomes/*.fasta ./reference_genomes/

```
	 

## Run mashmap
For sanity's sake, we will run a few comparisons that are listed below.

- pyrrhus to adamanteus; using hap1
- pyrrhus to viridis; using hap1

- pyrrhus to adamanteus; using hap2
- pyrrhus to viridis; using hap2

- hap1 to hap2



First, we make a mashmap results directory.

```
mkdir mashmap_results 

cd mashmap_results

# For haplotype 2
mkdir mashmap_results 

cd mashmap_results

```

Run mashmap; we will use three different percent identity thresholds (85, 90, 95). We compared across thresholds and found them to be conordant, so we used the 95% threshold.

Start with 85
```
# pyrrhus to adamanteus; hap1
mashmap -t 24 -r ../reference_genomes/C_adamanteus_genome_GCA_039797435.1.fasta -q ../../Hap1/Crotalus_pyrrhus__07-08-2024__final_assembly.fasta -f one-to-one -s 50000 --pi 85 -o pyrrhus_2_adamanteus_hap1_pi85.txt

# pyrrhus to viridis; hap1
mashmap -t 24 -r ../reference_genomes/CroVir_genome_L77pg_16Aug2017.final_rename.fasta -q ../../Hap1/Crotalus_pyrrhus__07-08-2024__final_assembly.fasta -f one-to-one -s 50000 --pi 85 -o pyrrhus_2_viridis_hap1_pi85.txt

# pyrrhus to adamanteus; hap2
mashmap -t 24 -r ../reference_genomes/C_adamanteus_genome_GCA_039797435.1.fasta -q ../../Hap2/Crotalus_pyrrhus__07-08-2024__final_assembly.fasta -f one-to-one -s 50000 --pi 85 -o pyrrhus_2_adamanteus_hap2_pi85.txt

# pyrrhus to viridis; hap2
mashmap -t 24 -r ../reference_genomes/CroVir_genome_L77pg_16Aug2017.final_rename.fasta -q ../../Hap2/Crotalus_pyrrhus__07-08-2024__final_assembly.fasta -f one-to-one -s 50000 --pi 85 -o pyrrhus_2_viridis_hap2_pi85.txt

# hap1 vs hap2
mashmap -t 24 -r ../../Hap2/Crotalus_pyrrhus__07-08-2024__final_assembly.fasta -q ../../Hap1/Crotalus_pyrrhus__07-08-2024__final_assembly.fasta -f one-to-one -s 50000 --pi 85 -o pyrrhus_hap1_2_hap2_pi85.txt
```

90 percent identity
```
# pyrrhus to adamanteus; hap1
mashmap -t 24 -r ../reference_genomes/C_adamanteus_genome_GCA_039797435.1.fasta -q ../../Hap1/Crotalus_pyrrhus__07-08-2024__final_assembly.fasta -f one-to-one -s 50000 --pi 90 -o pyrrhus_2_adamanteus_hap1_pi90.txt

# pyrrhus to viridis; hap1
mashmap -t 24 -r ../reference_genomes/CroVir_genome_L77pg_16Aug2017.final_rename.fasta -q ../../Hap1/Crotalus_pyrrhus__07-08-2024__final_assembly.fasta -f one-to-one -s 50000 --pi 90 -o pyrrhus_2_viridis_hap1_pi90.txt

# pyrrhus to adamanteus; hap2
mashmap -t 24 -r ../reference_genomes/C_adamanteus_genome_GCA_039797435.1.fasta -q ../../Hap2/Crotalus_pyrrhus__07-08-2024__final_assembly.fasta -f one-to-one -s 50000 --pi 90 -o pyrrhus_2_adamanteus_hap2_pi90.txt

# pyrrhus to viridis; hap2
mashmap -t 24 -r ../reference_genomes/CroVir_genome_L77pg_16Aug2017.final_rename.fasta -q ../../Hap2/Crotalus_pyrrhus__07-08-2024__final_assembly.fasta -f one-to-one -s 50000 --pi 90 -o pyrrhus_2_viridis_hap2_pi90.txt

# hap1 vs hap2
mashmap -t 24 -r ../../Hap2/Crotalus_pyrrhus__07-08-2024__final_assembly.fasta -q ../../Hap1/Crotalus_pyrrhus__07-08-2024__final_assembly.fasta -f one-to-one -s 50000 --pi 90 -o pyrrhus_hap1_2_hap2_pi90.txt
```

95 percent identity
```
# pyrrhus to adamanteus; hap1
mashmap -t 24 -r ../reference_genomes/C_adamanteus_genome_GCA_039797435.1.fasta -q ../../Hap1/Crotalus_pyrrhus__07-08-2024__final_assembly.fasta -f one-to-one -s 50000 --pi 95 -o pyrrhus_2_adamanteus_hap1_pi95.txt

# pyrrhus to viridis; hap1
mashmap -t 24 -r ../reference_genomes/CroVir_genome_L77pg_16Aug2017.final_rename.fasta -q ../../Hap1/Crotalus_pyrrhus__07-08-2024__final_assembly.fasta -f one-to-one -s 50000 --pi 95 -o pyrrhus_2_viridis_hap1_pi95.txt

# pyrrhus to adamanteus; hap2
mashmap -t 24 -r ../reference_genomes/C_adamanteus_genome_GCA_039797435.1.fasta -q ../../Hap2/Crotalus_pyrrhus__07-08-2024__final_assembly.fasta -f one-to-one -s 50000 --pi 95 -o pyrrhus_2_adamanteus_hap2_pi95.txt

# pyrrhus to viridis; hap2
mashmap -t 24 -r ../reference_genomes/CroVir_genome_L77pg_16Aug2017.final_rename.fasta -q ../../Hap2/Crotalus_pyrrhus__07-08-2024__final_assembly.fasta -f one-to-one -s 50000 --pi 95 -o pyrrhus_2_viridis_hap2_pi95.txt

# hap1 vs hap2
mashmap -t 24 -r ../../Hap2/Crotalus_pyrrhus__07-08-2024__final_assembly.fasta -q ../../Hap1/Crotalus_pyrrhus__07-08-2024__final_assembly.fasta -f one-to-one -s 50000 --pi 95 -o pyrrhus_hap1_2_hap2_pi95.txt
```

We will use a smaller window size (10kb) to explore the microchromosomes, but we will only use the 90 percent identity thresholds. First we make another directory to hold our results.

```
mkdir 10kb_window_run

# pyrrhus to adamanteus; hap1
mashmap -t 24 -r ../reference_genomes/C_adamanteus_genome_GCA_039797435.1.fasta -q ../../Hap1/Crotalus_pyrrhus__07-08-2024__final_assembly.fasta -f one-to-one -s 10000 --pi 90 -o ./10kb_window_run/pyrrhus_2_adamanteus_hap1_pi90_10kb.txt

# pyrrhus to viridis; hap1
mashmap -t 24 -r ../reference_genomes/CroVir_genome_L77pg_16Aug2017.final_rename.fasta -q ../../Hap1/Crotalus_pyrrhus__07-08-2024__final_assembly.fasta -f one-to-one -s 10000 --pi 90 -o ./10kb_window_run/pyrrhus_2_viridis_hap1_pi90_10kb.txt

# pyrrhus to adamanteus; hap2
mashmap -t 24 -r ../reference_genomes/C_adamanteus_genome_GCA_039797435.1.fasta -q ../../Hap2/Crotalus_pyrrhus__07-08-2024__final_assembly.fasta -f one-to-one -s 10000 --pi 90 -o ./10kb_window_run/pyrrhus_2_adamanteus_hap2_pi90_10kb.txt

# pyrrhus to viridis; hap2
mashmap -t 24 -r ../reference_genomes/CroVir_genome_L77pg_16Aug2017.final_rename.fasta -q ../../Hap2/Crotalus_pyrrhus__07-08-2024__final_assembly.fasta -f one-to-one -s 10000 --pi 90 -o ./10kb_window_run/pyrrhus_2_viridis_hap2_pi90_10kb.txt

# hap1 vs hap2
mashmap -t 24 -r ../../Hap2/Crotalus_pyrrhus__07-08-2024__final_assembly.fasta -q ../../Hap1/Crotalus_pyrrhus__07-08-2024__final_assembly.fasta -f one-to-one -s 10000 --pi 90 -o ./10kb_window_run/pyrrhus_hap1_2_hap2_pi90_10kb.txt
```
