# crotalus_introgression_wholegenome

*Last updated May 8th, 2026*

**WARNING, LINKS MAY NOT WORK UNTIL AFTER THE MANUSCRIPT HAS BEEN PUBLISHED.**

I am in the process of uploading scripts to this repository to ensure that they are properly annotated, but please reach out if you have any questions.

This repository contains the computational workflow and scripts for [Farleigh et al., (*accepted*)](). Please email Keaka Farleigh (keakafarleigh@virginia.edu; keakafarleigh@gmail.com) if you have any questions. 

## Citation

The citation will be here.

## Genomic data

The genomic data for this project are deposited as parts of BioProjects [PRJNA1454467](), [PRJNA593834](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA593834), and [PRJNA1150930](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1150930). The genome for this project is deposited at [PRJNA1460337]().

## Structure of this repository.

This repository contains the computational workflow for [Farleigh et al., (*accepted*)](). This ReadMe will contain the entire computational workflow for the project and serve as a sort of one-stop shop. It will also tell you which scripts are associated with each analysis so that you can download pipelines for specific analyses if you wish. I have separated this ReadMe into each of the method sections listed in the [supplemental information](). Each section begins with a general description of the methods and software that are relevant before moving into the workflow. Each script referenced herein is located in the `scripts/` directory.

You can jump to any of the relevant sections by clicking on the links below:
- [Genome assembly and annotation](https://github.com/kfarleigh/crotalus_introgression_wholegenome/tree/main#genome-assembly-and-annotation)
- [Whole genome sequencing and variant calling](https://github.com/kfarleigh/crotalus_introgression_wholegenome/tree/main#whole-genome-sequencing-and-variant-calling)
- [Phylogeny, population structure, and demography](https://github.com/kfarleigh/crotalus_introgression_wholegenome/tree/main#phylogeny-population-structure-and-demography)
- [Recombination rate and recombination hotspot identification](https://github.com/kfarleigh/crotalus_introgression_wholegenome/tree/main#recombination-rate-and-recombination-hotspot-identification)
- [Calculation of introgression and divergence statistics](https://github.com/kfarleigh/crotalus_introgression_wholegenome/tree/main#calculation-of-introgression-and-divergence-statistics)
- [Distinguishing introgression from incomplete lineage sorting](https://github.com/kfarleigh/crotalus_introgression_wholegenome/tree/main#distinguishing-introgression-from-incomplete-lineage-sorting)
- [Evolutionary and ecological divergence](https://github.com/kfarleigh/crotalus_introgression_wholegenome/tree/main#evolutionary-and-ecological-divergence)
- [Relationships between introgression, recombination, and divergent selection](https://github.com/kfarleigh/crotalus_introgression_wholegenome/tree/main#relationships-between-introgression-recombination-and-divergent-selection)
- [Identification of species barriers and fine-scale variation in introgression](https://github.com/kfarleigh/crotalus_introgression_wholegenome/tree/main#identification-of-species-barriers-and-fine-scale-variation-in-introgression)

### Genome assembly and annotation

We assigned chromosome names using a synteny-based approach as implemented in [mashmap](https://github.com/marbl/MashMap). The script to perform this is called `chromosome_assignment.md`.

### Whole genome sequencing and variant calling

We sequenced newly generated libraries on Illumina NovaSeq 6000 lanes to generate 150 bp paired-end reads. We then trimmed our sequence data with [Trimmomatic](https://github.com/usadellab/trimmomatic). We then followed [GATK's best practices workflow](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows) to call and filter SNPs. The pipeline to perform this set of analyses is called `variant_calling.md`.

*Software for Whole genome sequencing and variant calling*:
- [Trimmomatic](https://github.com/usadellab/trimmomatic)
- [GATK](https://gatk.broadinstitute.org/hc/en-us)
- [bwa mem](https://github.com/lh3/BWA)
- [samtools](https://github.com/samtools/samtools)
- [bcftools](https://github.com/samtools/bcftools)
- [parallel](https://www.gnu.org/software/parallel/)
- [tabix](https://www.htslib.org/doc/tabix.html)
- [mosdepth](https://github.com/brentp/mosdepth)
- [bedtools](https://bedtools.readthedocs.io/en/latest/)
- [picard](https://github.com/broadinstitute/picard)

#### Workflow for analysis of the genomic landscape of introgression within Crotalus - Part 1 - Processing, mapping, and variant calling

Author: Keaka Farleigh, edited from Drew Schield's workflow
Date: October 15th, 2024
Email: keakafarleigh@virginia.edu

##### Overview

In this part of the workflow, we'll quality trim raw reads, perform mapping to the reference genome, and call variants.

------------------------------------------------------------------------------------------
##### Set up environment

We'll process the data using the directories on the `extradrive1` drive:

1. First we need to create the directories on the `extradrive1` drive in the `crotalus_genomic_landscape` folder. *These paths will be updated after final filtering*


```
mkdir /media/keaka/extradrive5/crotalus_genomic_landscape/fastq
mkdir /media/keaka/extradrive5/crotalus_genomic_landscape/fastq_filtered
mkdir /media/keaka/extradrive5/crotalus_genomic_landscape/bam
mkdir /media/keaka/extradrive5/crotalus_genomic_landscape/gvcf
mkdir /media/keaka/extradrive5/crotalus_genomic_landscape/vcf
mkdir /media/keaka/extradrive5/crotalus_genomic_landscape/analysis
```

See below for a list of what each directory will contain:
- `/media/keaka/extradrive5/crotalus_genomic_landscape/fastq`
	- Raw fastq's used as input

- `/media/keaka/extradrive5/crotalus_genomic_landscape/fastq_filtered`
	- Fastq files that have been filtered with [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)

- `/media/keaka/extradrive5/crotalus_genomic_landscape/bam`
	- Bam files that were generated after mapping the filtered fastq files to the pyrrhus reference genome

- `/media/keaka/extradrive5/crotalus_genomic_landscape/gvcf`
	- gVCF files that were generated after using [GATK's](https://gatk.broadinstitute.org/hc/en-us)  `HaplotypeCaller` and `GenotypeGVCFs`
	
- `/media/keaka/extradrive5/crotalus_genomic_landscape/vcf`
	- VCF file that was filtered using [bcftools](https://samtools.github.io/bcftools/bcftools.html)

Get into the working directory:
```
cd /media/keaka/extradrive5/crotalus_genomic_landscape/
mkdir log
```

##### File formatting
Some of the files we received from the Castoe group had a .fq.gz extension instead of .fastq.gz. This is not a big deal, but we will rename those files to match everything else.

```
cd /media/keaka/extradrive5/crotalus_genomic_landscape/fastq

rename "s/.fq.gz/.fastq.gz/" *.fq.gz
 
```

There are also all kind of ways that individuals are named here, so we will format the file names so that they use a common system. Sample name will precede the underscore, followed by whether it is the forward (R1) or reverse (R2) read.

```
# Preview what we want to do. 
ls *R1*.fastq.gz | cut -d_ -f1

# Do it using a for loop; for the R1s
for i in `ls *R1*.fastq.gz`; do mv $i `echo $i | cut -d_ -f1`_R1_.fastq.gz; done

# Now we do it for the R2s 
for i in `ls *R2*.fastq.gz`; do mv $i `echo $i | cut -d_ -f1`_R2_.fastq.gz; done


```

Now, we can move on to quality trimming.

------------------------------------------------------------------------------------------
##### 1. Quality trimming using Trimmomatic

###### 1. Format sample list

We will change directories into the fastq directory. Then, we will create the `listTrimmomatic.txt` which contains the sample names without the file extensions. 

```
cd /media/keaka/extradrive5/crotalus_genomic_landscape/fastq

ls *R1* | cut -d_ -f1 >> listTrimmomatic.txt

# Run trimmomatic on lane 1 of 2 from Oregon
cd /media/queen/extradrive1/crotalus_genomic_landscape/fastq

# Create a list of the samples for trimmomatic, I remove the undetermined from this list
ls *R1* | cut -d_ -f1,2,3 >> ../listTrimmomatic_UO.txt 
```

`listTrimmomatic.txt`

Note: In preliminary runs, samples RS_1 and RS_4 stopped prematurely due to corruption in the compression of the read 2 file (tested using `gunzip -t $file`). We will run these as single end analyses.

###### 2. Format script to run with `parallel`, calling the sample list

`runTrimmomatic_parallel.sh`:

```
touch runTrimmomatic_parallel.sh

# For Oregon, newly sequenced samples 
cd /media/queen/extradrive1/crotalus_genomic_landscape

touch runTrimmomatic_parallel_UO.sh
```

Copy and paste code below into `runTrimmomatic_parallel.sh`

```
indv=$1
trimmomatic PE -phred33 -threads 6 /media/keaka/extradrive5/crotalus_genomic_landscape/fastq/${indv}_R1_.fastq.gz /media/keaka/extradrive5/crotalus_genomic_landscape/fastq/${indv}_R2_.fastq.gz /media/keaka/extradrive5/crotalus_genomic_landscape/fastq_filtered/${indv}_1_P.trim.fq.gz /media/keaka/extradrive5/crotalus_genomic_landscape/fastq_filtered/${indv}_1_U.trim.fq.gz /media/keaka/extradrive5/crotalus_genomic_landscape/fastq_filtered/${indv}_2_P.trim.fq.gz /media/keaka/extradrive5/crotalus_genomic_landscape/fastq_filtered/${indv}_2_U.trim.fq.gz LEADING:20 TRAILING:20 MINLEN:32 AVGQUAL:30 > ./log/runTrimmomatic.$indv.log

# If it is for `runTrimmomatic_parallel_UO_lane1.sh`

indv=$1
trimmomatic PE -phred33 -threads 4 /media/queen/extradrive1/crotalus_genomic_landscape/fastq/${indv}_R1_001.fastq.gz /media/queen/extradrive1/crotalus_genomic_landscape/fastq/${indv}_R2_001.fastq.gz /media/queen/extradrive1/crotalus_genomic_landscape/fastq_filtered/${indv}_1_P.trim.fq.gz /media/queen/extradrive1/crotalus_genomic_landscape/fastq_filtered/${indv}_1_U.trim.fq.gz /media/queen/extradrive1/crotalus_genomic_landscape/fastq_filtered/${indv}_2_P.trim.fq.gz /media/queen/extradrive1/crotalus_genomic_landscape/fastq_filtered/${indv}_2_U.trim.fq.gz LEADING:20 TRAILING:20 MINLEN:32 AVGQUAL:30 > ./log/runTrimmomatic.$indv.log
```

###### 3. Run script with `parallel`

```
cd /media/keaka/extradrive5/crotalus_genomic_landscape
parallel --progress --joblog ./log/logfile.trimmomatic -j 8 --workdir . ./runTrimmomatic_parallel.sh :::: listTrimmomatic.txt

# The power went out and Nostromo died, so we have to restart trimmomatic where it left off using listTrimmomatic_v2.txt
parallel --progress --joblog ./log/logfile.trimmomatic -j 8 --workdir . ./runTrimmomatic_parallel.sh :::: listTrimmomatic_v2.txt


# Run for the Oregon lanes and the outgroup samples (atrox, ruber)
/media/queen/extradrive1/crotalus_genomic_landscape

nohup parallel --progress --joblog ./log/logfile_UO.trimmomatic -j 7 --workdir . ./runTrimmomatic_parallel_UO.sh :::: listTrimmomatic_UO.txt &
```

###### 4. Remove unpaired reads from output directory to save disk space

```
# For sequences from previously published data and the Castoe lab
rm /media/keaka/extradrive5/crotalus_genomic_landscape/fastq_filtered/*U.trim.fq.gz

# For the samples sequenced in the Scheild lab
rm /media/queen/extradrive1/crotalus_genomic_landscape/fastq_filtered/*U.trim.fq.gz
rm /media/queen/extradrive1/crotalus_genomic_landscape/fastq_filtered/Undetermined*
```

###### 5. Rename trimmed reads from the University of Oregon

```
for i in `ls *_1_P.trim.fq.gz `; do mv $i `echo $i | cut -d_ -f1`_1_P.trim.fq.gz; done
for i in `ls *_2_P.trim.fq.gz `; do mv $i `echo $i | cut -d_ -f1`_2_P.trim.fq.gz; done
```

------------------------------------------------------------------------------------------
##### 2. Mapping to reference genome

We'll map the trimmed reads to the chromosome-assigned reference genome `/media/queen/extradrive1/crotalus_pyrrhus_genome/bwa_index/Crotalus_pyrrhus.hap1.final.fasta`.

###### 1. Format script to run `bwa` to map and index bam files

`runBWA_UO_parallel.sh`

```
indv=$1
bwa mem -t 4 -R "@RG\tID:$indv\tLB:Crotalus\tPL:illumina\tPU:NovaSeq6000\tSM:$indv" /media/queen/extradrive1/crotalus_pyrrhus_genome/bwa_index/Crotalus_pyrrhus.hap1.final.fasta /media/queen/extradrive1/crotalus_genomic_landscape/fastq_filtered/${indv}_1_P.trim.fq.gz /media/queen/extradrive1/crotalus_genomic_landscape/fastq_filtered/${indv}_2_P.trim.fq.gz | samtools sort -@ 4 -O bam -T $indv.temp -o /media/queen/TombBucket/crotalus_genomic_landscape/bam/$indv.bam -
samtools index -@ 4 /media/queen/TombBucket/crotalus_genomic_landscape/bam/$indv.bam
```

###### 2. Run script with `parallel`

```
# Create the list for BWA
cd /media/queen/extradrive1/crotalus_genomic_landscape
cat listTrimmomatic_UO.txt | cut -d_ -f1 > listBWA_UO.txt
 
parallel --progress --joblog ./log/logfile.bwa -j 10 --workdir . ./runBWA_UO_parallel.sh :::: listBWA_UO.txt > runBWA_UO_lane1.log
```

------------------------------------------------------------------------------------------
##### 3. Variant calling

We'll call genomic variants using `gatk`, starting with individual variant calls using `HaplotypeCaller`, followed by cohort variant calls using `GenotypeGVCFs`.

###### Format sample list

`listGATK.txt`

###### 1. Individual variant calls using HaplotypeCaller

####### 1. Create listGATK.txt

```
cd /media/queen/TombBucket/crotalus_genomic_landscape/bam

ls *.bam | cut -d .  -f1 >> ../listGATK.txt
 
cd ../
```

####### 2. Format script to run GATK HaplotypeCaller

`runGATKHaplotypeCaller_parallel.sh`

```
# Install gatk4, we will create a gatk environment 
conda create -n gatk4
conda activate gatk4
conda install conda-forge::mamba
mamba install gatk4

# Create the haplotypecaller shell script 
nano runGATKHaplotypeCaller_parallel.sh

# Copy and paste the code below into the runGATKHaplotypeCaller_parallel.sh file

indv=$1
gatk --java-options "-Xmx12g -XX:ConcGCThreads=1" HaplotypeCaller -R /media/queen/extradrive1/crotalus_pyrrhus_genome/Crotalus_pyrrhus.hap1.final.fasta --ERC GVCF -I /media/queen/TombBucket/crotalus_genomic_landscape/bam/$indv.bam -O /media/queen/TombBucket/crotalus_genomic_landscape/gvcf/$indv.raw.snps.indels.chr1.g.vcf -L chr1_scaffold_3_3contigs > ./log/runGATKHaplotypeCaller.$indv.log
gatk --java-options "-Xmx12g -XX:ConcGCThreads=1" HaplotypeCaller -R /media/queen/extradrive1/crotalus_pyrrhus_genome/Crotalus_pyrrhus.hap1.final.fasta --ERC GVCF -I /media/queen/TombBucket/crotalus_genomic_landscape/bam/$indv.bam -O /media/queen/TombBucket/crotalus_genomic_landscape/gvcf/$indv.raw.snps.indels.chr2.g.vcf -L chr2_scaffold_1_1contigs > ./log/runGATKHaplotypeCaller.$indv.log
gatk --java-options "-Xmx12g -XX:ConcGCThreads=1" HaplotypeCaller -R /media/queen/extradrive1/crotalus_pyrrhus_genome/Crotalus_pyrrhus.hap1.final.fasta --ERC GVCF -I /media/queen/TombBucket/crotalus_genomic_landscape/bam/$indv.bam -O /media/queen/TombBucket/crotalus_genomic_landscape/gvcf/$indv.raw.snps.indels.chr3.g.vcf -L chr3_scaffold_2_1contigs > ./log/runGATKHaplotypeCaller.$indv.log
gatk --java-options "-Xmx12g -XX:ConcGCThreads=1" HaplotypeCaller -R /media/queen/extradrive1/crotalus_pyrrhus_genome/Crotalus_pyrrhus.hap1.final.fasta --ERC GVCF -I /media/queen/TombBucket/crotalus_genomic_landscape/bam/$indv.bam -O /media/queen/TombBucket/crotalus_genomic_landscape/gvcf/$indv.raw.snps.indels.chr4.g.vcf -L chr4_scaffold_5_1contigs > ./log/runGATKHaplotypeCaller.$indv.log
gatk --java-options "-Xmx12g -XX:ConcGCThreads=1" HaplotypeCaller -R /media/queen/extradrive1/crotalus_pyrrhus_genome/Crotalus_pyrrhus.hap1.final.fasta --ERC GVCF -I /media/queen/TombBucket/crotalus_genomic_landscape/bam/$indv.bam -O /media/queen/TombBucket/crotalus_genomic_landscape/gvcf/$indv.raw.snps.indels.chr5.g.vcf -L chr5_scaffold_6_1contigs > ./log/runGATKHaplotypeCaller.$indv.log
gatk --java-options "-Xmx12g -XX:ConcGCThreads=1" HaplotypeCaller -R /media/queen/extradrive1/crotalus_pyrrhus_genome/Crotalus_pyrrhus.hap1.final.fasta --ERC GVCF -I /media/queen/TombBucket/crotalus_genomic_landscape/bam/$indv.bam -O /media/queen/TombBucket/crotalus_genomic_landscape/gvcf/$indv.raw.snps.indels.chr6.g.vcf -L chr6_scaffold_7_1contigs > ./log/runGATKHaplotypeCaller.$indv.log
gatk --java-options "-Xmx12g -XX:ConcGCThreads=1" HaplotypeCaller -R /media/queen/extradrive1/crotalus_pyrrhus_genome/Crotalus_pyrrhus.hap1.final.fasta --ERC GVCF -I /media/queen/TombBucket/crotalus_genomic_landscape/bam/$indv.bam -O /media/queen/TombBucket/crotalus_genomic_landscape/gvcf/$indv.raw.snps.indels.chr7.g.vcf -L chr7_scaffold_8_1contigs > ./log/runGATKHaplotypeCaller.$indv.log
gatk --java-options "-Xmx12g -XX:ConcGCThreads=1" HaplotypeCaller -R /media/queen/extradrive1/crotalus_pyrrhus_genome/Crotalus_pyrrhus.hap1.final.fasta --ERC GVCF -I /media/queen/TombBucket/crotalus_genomic_landscape/bam/$indv.bam -O /media/queen/TombBucket/crotalus_genomic_landscape/gvcf/$indv.raw.snps.indels.chr8.g.vcf -L chr8_scaffold_10_1contigs > ./log/runGATKHaplotypeCaller.$indv.log
gatk --java-options "-Xmx12g -XX:ConcGCThreads=1" HaplotypeCaller -R /media/queen/extradrive1/crotalus_pyrrhus_genome/Crotalus_pyrrhus.hap1.final.fasta --ERC GVCF -I /media/queen/TombBucket/crotalus_genomic_landscape/bam/$indv.bam -O /media/queen/TombBucket/crotalus_genomic_landscape/gvcf/$indv.raw.snps.indels.chr9.g.vcf -L chr9_scaffold_9_1contigs > ./log/runGATKHaplotypeCaller.$indv.log
gatk --java-options "-Xmx12g -XX:ConcGCThreads=1" HaplotypeCaller -R /media/queen/extradrive1/crotalus_pyrrhus_genome/Crotalus_pyrrhus.hap1.final.fasta --ERC GVCF -I /media/queen/TombBucket/crotalus_genomic_landscape/bam/$indv.bam -O /media/queen/TombBucket/crotalus_genomic_landscape/gvcf/$indv.raw.snps.indels.chr10.g.vcf -L chr10_scaffold_13_1contigs > ./log/runGATKHaplotypeCaller.$indv.log
gatk --java-options "-Xmx12g -XX:ConcGCThreads=1" HaplotypeCaller -R /media/queen/extradrive1/crotalus_pyrrhus_genome/Crotalus_pyrrhus.hap1.final.fasta --ERC GVCF -I /media/queen/TombBucket/crotalus_genomic_landscape/bam/$indv.bam -O /media/queen/TombBucket/crotalus_genomic_landscape/gvcf/$indv.raw.snps.indels.chr11.g.vcf -L chr11_scaffold_11_1contigs > ./log/runGATKHaplotypeCaller.$indv.log
gatk --java-options "-Xmx12g -XX:ConcGCThreads=1" HaplotypeCaller -R /media/queen/extradrive1/crotalus_pyrrhus_genome/Crotalus_pyrrhus.hap1.final.fasta --ERC GVCF -I /media/queen/TombBucket/crotalus_genomic_landscape/bam/$indv.bam -O /media/queen/TombBucket/crotalus_genomic_landscape/gvcf/$indv.raw.snps.indels.chr12.g.vcf -L chr12_scaffold_14_1contigs > ./log/runGATKHaplotypeCaller.$indv.log
gatk --java-options "-Xmx12g -XX:ConcGCThreads=1" HaplotypeCaller -R /media/queen/extradrive1/crotalus_pyrrhus_genome/Crotalus_pyrrhus.hap1.final.fasta --ERC GVCF -I /media/queen/TombBucket/crotalus_genomic_landscape/bam/$indv.bam -O /media/queen/TombBucket/crotalus_genomic_landscape/gvcf/$indv.raw.snps.indels.chr13.g.vcf -L chr13_scaffold_12_1contigs > ./log/runGATKHaplotypeCaller.$indv.log
gatk --java-options "-Xmx12g -XX:ConcGCThreads=1" HaplotypeCaller -R /media/queen/extradrive1/crotalus_pyrrhus_genome/Crotalus_pyrrhus.hap1.final.fasta --ERC GVCF -I /media/queen/TombBucket/crotalus_genomic_landscape/bam/$indv.bam -O /media/queen/TombBucket/crotalus_genomic_landscape/gvcf/$indv.raw.snps.indels.chr14.g.vcf -L chr14_scaffold_16_1contigs > ./log/runGATKHaplotypeCaller.$indv.log
gatk --java-options "-Xmx12g -XX:ConcGCThreads=1" HaplotypeCaller -R /media/queen/extradrive1/crotalus_pyrrhus_genome/Crotalus_pyrrhus.hap1.final.fasta --ERC GVCF -I /media/queen/TombBucket/crotalus_genomic_landscape/bam/$indv.bam -O /media/queen/TombBucket/crotalus_genomic_landscape/gvcf/$indv.raw.snps.indels.chr15.g.vcf -L chr15_scaffold_18_1contigs > ./log/runGATKHaplotypeCaller.$indv.log
gatk --java-options "-Xmx12g -XX:ConcGCThreads=1" HaplotypeCaller -R /media/queen/extradrive1/crotalus_pyrrhus_genome/Crotalus_pyrrhus.hap1.final.fasta --ERC GVCF -I /media/queen/TombBucket/crotalus_genomic_landscape/bam/$indv.bam -O /media/queen/TombBucket/crotalus_genomic_landscape/gvcf/$indv.raw.snps.indels.chr16.g.vcf -L chr16_scaffold_15_1contigs > ./log/runGATKHaplotypeCaller.$indv.log
gatk --java-options "-Xmx12g -XX:ConcGCThreads=1" HaplotypeCaller -R /media/queen/extradrive1/crotalus_pyrrhus_genome/Crotalus_pyrrhus.hap1.final.fasta --ERC GVCF -I /media/queen/TombBucket/crotalus_genomic_landscape/bam/$indv.bam -O /media/queen/TombBucket/crotalus_genomic_landscape/gvcf/$indv.raw.snps.indels.chr17.g.vcf -L chr17_scaffold_17_1contigs > ./log/runGATKHaplotypeCaller.$indv.log
gatk --java-options "-Xmx12g -XX:ConcGCThreads=1" HaplotypeCaller -R /media/queen/extradrive1/crotalus_pyrrhus_genome/Crotalus_pyrrhus.hap1.final.fasta --ERC GVCF -I /media/queen/TombBucket/crotalus_genomic_landscape/bam/$indv.bam -O /media/queen/TombBucket/crotalus_genomic_landscape/gvcf/$indv.raw.snps.indels.chrZ.g.vcf -L chrZ_scaffold_4_1contigs > ./log/runGATKHaplotypeCaller.$indv.log
gatk --java-options "-Xmx12g -XX:ConcGCThreads=1" HaplotypeCaller -R /media/queen/extradrive1/crotalus_pyrrhus_genome/Crotalus_pyrrhus.hap1.final.fasta --ERC GVCF -I /media/queen/TombBucket/crotalus_genomic_landscape/bam/$indv.bam -O /media/queen/TombBucket/crotalus_genomic_landscape/gvcf/$indv.raw.snps.indels.upscaf.g.vcf -L unplaced_scaf_list.list > ./log/runGATKHaplotypeCaller.$indv.log

cd /media/queen/TombBucket/crotalus_genomic_landscape/gvcf/

gatk --java-options "-Xmx12g -XX:ConcGCThreads=1" MergeVcfs I=$indv.raw.snps.indels.chr1.g.vcf  I=$indv.raw.snps.indels.chr2.g.vcf  I=$indv.raw.snps.indels.chr3.g.vcf  I=$indv.raw.snps.indels.chr4.g.vcf  I=$indv.raw.snps.indels.chr5.g.vcf  I=$indv.raw.snps.indels.chr6.g.vcf  I=$indv.raw.snps.indels.chr7.g.vcf  I=$indv.raw.snps.indels.chr8.g.vcf  I=$indv.raw.snps.indels.chr9.g.vcf  I=$indv.raw.snps.indels.chr10.g.vcf  I=$indv.raw.snps.indels.chr11.g.vcf  I=$indv.raw.snps.indels.chr12.g.vcf  I=$indv.raw.snps.indels.chr13.g.vcf  I=$indv.raw.snps.indels.chr14.g.vcf  I=$indv.raw.snps.indels.chr15.g.vcf  I=$indv.raw.snps.indels.chr16.g.vcf  I=$indv.raw.snps.indels.chr17.g.vcf  I=$indv.raw.snps.indels.chrZ.g.vcf  I=$indv.raw.snps.indels.upscaf.g.vcf O=$indv.raw.snps.indels.g.vcf

cd /media/queen/TombBucket/crotalus_genomic_landscape/

bgzip /media/queen/TombBucket/crotalus_genomic_landscape/gvcf/$indv.raw.snps.indels.g.vcf

rm /media/queen/TombBucket/crotalus_genomic_landscape/gvcf/$indv.raw.snps.indels.chr*.g.vcf /media/queen/TombBucket/crotalus_genomic_landscape/gvcf/,$indv.raw.snps.indels.upscaf.g.vcf

# Index out fasta
samtools faidx /media/queen/extradrive1/crotalus_pyrrhus_genome/Crotalus_pyrrhus.hap1.final.fasta

# Create a sequence dictionary for our reference
gatk CreateSequenceDictionary -R /media/queen/extradrive1/crotalus_pyrrhus_genome/Crotalus_pyrrhus.hap1.final.fasta -O /media/queen/extradrive1/crotalus_pyrrhus_genome/Crotalus_pyrrhus.hap1.final.dict

```

####### 3. Run script with `parallel`

```
# Run it 
cd /media/queen/TombBucket/crotalus_genomic_landscape

cat intervallist.txt | while read line; do echo parallel --progress --joblog ./log/logfile.$line.GATKHaplotypeCaller -j 14  --workdir . ./Crotalus_haplotypecaller_bychrom.sh {1} $line :::: list_GATK_v2.txt >> run_Crotalus_haplotypecaller_bychrom.sh; done

# Do a dry run to make sure it looks okay
parallel --dry-run --progress --joblog ./log/logfile.chr1_scaffold_3_3contigs.GATKHaplotypeCaller -j 14 --workdir . ./Crotalus_haplotypecaller_bychrom.sh {1} chr1_scaffold_3_3contigs :::: list_GATK_v2.txt


# Run SNP calling by chromosome
./run_Crotalus_haplotypecaller_bychrom.sh


# Merge the vcfs for each individual 
parallel --progress --joblog ./log/logfile.mergevcfs -j 12 --workdir . ./Crotalus_merge_vcfs.sh :::: list_GATK_v2.txt

# Store the individual chromosome vcfs temporarily
mkdir tmp_stor_chromvcfs

mv *chr*.vcf* ./tmp_stor_chromvcfs
mv *unplaced_scaf_list.list*.vcf* ./tmp_stor_chromvcfs


### Run on our second batch of sequencing
./run_Crotalus_haplotypecaller_bychrom_batch2.sh

# Merge the vcfs for each individual 
parallel --progress --joblog ./log/logfile.mergevcfs -j 12 --workdir . ./Crotalus_merge_vcfs.sh :::: list_GATK_batch2.txt

mv *chr*.vcf* ./tmp_stor_chromvcfs
mv *unplaced_scaf_list.list*.vcf* ./tmp_stor_chromvcfs

# We don't want to include the embryos in this study, so we will move them to another directory

mkdir embyros

mv *EMB* ./embryos/

# Make a final list to use in tabix indexing and genotype gvcfs


ls *.vcf.gz | cut -d . -f 1 >> ../list_GATK_allsamples.txt

cd ../

```

####### 4. Tabix index output gVCFs

```
for i in `cat list_GATK_allsamples.txt`; do tabix -p vcf /media/queen/TombBucket/crotalus_genomic_landscape/gvcf/$i.raw.snps.indels.g.vcf.gz -f; done
```

###### 2. Cohort variant calls using GenotypeGVCFs

###### Format sample list (with paths to gVCF files)

Generate a sample path with the full list.
``` 
ls -1 -d "$PWD/gvcf/"*.vcf.gz > gVCFlist.list 
```

###### Install gatk3

Note: We will create a new conda environment and install gatk3. This allows us to supply a list of .g.vcfs to GenotypeGVCFs. 

```
conda deactivate 

conda create -n gatk3
conda activate gatk3
conda install conda-forge::mamba
mamba install gatk
```

###### Run GenotypeGVCFs on sets of intervals


Generate the command to perform genotyping with GenotypeGVCFs. We escape the quotes (using `\`) so that echo will print them out in the command. 

```
cat intervallist.txt | while read line; do echo "gatk3 \"-Xmx32G\" \"-XX:ConcGCThreads=1\" -T GenotypeGVCFs -R /media/queen/extradrive1/crotalus_pyrrhus_genome/Crotalus_pyrrhus.hap1.final.fasta -L $line -V ./gVCFlist.list -allSites -o /media/queen/TombBucket/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.raw.${line}.vcf.gz > ./log/GenotypeGVCFs.crotalus_genus.allsites.raw.${line}.vcf.gz.log" >> run_Crotalus_GenotypeGVCFs.sh; done 

chmod +x run_Crotalus_GenotypeGvcfs.sh 

nohup sh ./run_Crotalus_GenotypeGvcfs.sh &> ./log/crotalus_genotypyegvcfs.out &

# You can check to make sure there are invdividuals in your vcf while GenotypeGCVFs is running

cd ./vcf/

# The number of samples in the vcf, remove the wc -l if you want a list
query -l crotalus_genus.allsites.raw.chr1_scaffold_3_3contigs.vcf.gz | wc -l

# Get the number of SNPs
zgrep -v "^#" crotalus_genus.allsites.raw.chr1_scaffold_3_3contigs.vcf.gz | wc -l

cd ../

```


------------------------------------------------------------------------------------------
##### 4. Sex identification

We'll use relative read depths on the Z chromosome and autosomes to infer genetic sex of each individual.

We'll extract mapping data for Chromosome 4, taking the median as an 'autosomal median', then compare levels of coverage for individuals to this value across the Z chromosome scaffold.


### ###Set up environment

```
cd /media/queen/extradrive1/crotalus_genomic_landscape/analysis

mkdir sex_identification
cd sex_identification
mkdir mosdepth_results
```

###### 1. Format sample and scaffold input files

Make a list of the samples: `mosdepth_samplelist_batch1.txt`.

```
(cd /media/queen/TombBucket/crotalus_genomic_landscape/bam/ && ls *.bam -1) > mosdepth_samplelist_batch1.txt
```

Make 'genome' files for chromosome 4 and the Z chromosome, with scaffold ID and length:

```
nano chrom.chr4.genome
chr4_scaffold_5_1contigs	126994310


nano chrom.chrZ.genome
chrZ_scaffold_4_1contigs	141807340
```

Converted genome files to BED sliding window files, in 10 kb windows:

```
bedtools makewindows -g chrom.chr4.genome -w 10000 > chrom.chr4.10kb.bed
head chrom.chr4.10kb.bed

bedtools makewindows -g chrom.chrZ.genome -w 10000 > chrom.chrZ.10kb.bed
head chrom.chrZ.10kb.bed
```

###### 2. Write Python script to quantify autosome median and log2 Z/Autosome depth ratio

`identifySex.py`

```
import sys
from numpy import array
from numpy import median
from numpy import mean
from numpy import log2

auto = []
za_norm = []

for line in open(sys.argv[1], 'r'):
	acov = float(line.split()[3])
	auto.append(acov)

auto_med = median(array(auto))

for line in open(sys.argv[2], 'r'):
	zcov = line.split()[3]
	if float(zcov) > 0.00:
		za = log2(float(zcov)/float(auto_med))
		za_norm.append(za)

za_mean = mean(array(za_norm))

sample = str(sys.argv[1])
sample = sample.split('/')[1]
sample = sample.split('.')[0]

if za_mean < -0.2:
	message = 'is likely female'
	print za_mean, sample, message
else:
	message = 'is likely male'
	print za_mean, sample, message
```

###### 3. Write wrapper script to run mosdepth and identifySex.py

`runSexIdentification.sh`

```
list=$1
for i in `cat $list`; do
	mosdepth -t 4 --fast-mode -n -b chrom.chr4.10kb.bed -f /media/queen/extradrive1/crotalus_pyrrhus_genome/Crotalus_pyrrhus.hap1.final.fasta ./mosdepth_results/$i.chr4 /media/queen/TombBucket/crotalus_genomic_landscape/bam/$i
	mosdepth -t 4 --fast-mode -n -b chrom.chrZ.10kb.bed -f /media/queen/extradrive1/crotalus_pyrrhus_genome/Crotalus_pyrrhus.hap1.final.fasta ./mosdepth_results/$i.chrZ /media/queen/TombBucket/crotalus_genomic_landscape/bam/$i
	gunzip ./mosdepth_results/$i.chr4.regions.bed.gz 
	gunzip ./mosdepth_results/$i.chrZ.regions.bed.gz
	python identifySex.py mosdepth_results/$i.chr4.regions.bed mosdepth_results/$i.chrZ.regions.bed
done
```

###### 4. Run script

```
nohup sh runSexIdentification.sh mosdepth_samplelist_batch1.txt > runSexIdentification.log &
```

###### 5. Extract general coverage information from mosdepth results

```
for indv in `cat mosdepth_samplelist_batch1.txt`; do mean=`awk '{ sum += $4; n++ } END { if (n > 0) print sum / n; }' ./mosdepth_results/$indv.chr4.regions.bed`; echo $indv $mean; done
```

------------------------------------------------------------------------------------------
##### 5. Variant filtering

In this section we'll impose various filtering steps to produce high-quality SNPs for downstream analysis.

###### 1. Hard filters in GATK

We'll first flag genotypes to be filtered based on the following specs:
* QD < 2.0
* FS > 60.0
* MQ < 40.0
* MQRankSum < -12.5
* ReadPosRankSum < -8.0

Then we'll set the annotated genotypes as missing in the all-sites VCFs.

Run GATK VariantFiltration to annotate genotypes that don't pass the hard filters.

```
conda activate gatk4 
cd /media/queen/TombBucket/crotalus_genomic_landscape

cd ./vcf

# Create an .sh file with the commands
for i in `ls *.vcf.gz | sed 's/.vcf.gz//'`; do echo "gatk --java-options \"-Xmx32G -XX:ConcGCThreads=1\" VariantFiltration -V ./vcf/${i}.vcf.gz -filter \"QD < 2.0\" --filter-name \"QD2\" -filter \"FS > 60.0\" --filter-name \"FS60\" -filter \"MQ < 40.0\" --filter-name \"MQ40\" -filter \"MQRankSum < -12.5\" --filter-name \"MQRankSum-12.5\" -filter \"ReadPosRankSum < -8.0\" --filter-name \"ReadPosRankSum-8\" -O /media/queen/TombBucket/crotalus_genomic_landscape/vcf/${i}.HardFilter.vcf.gz > ./log/GATK_VariantFilteration.${i}.log" >> ../run_Crotalus_VariantFilter.sh; done

cd .. 

# Make it executable
chmod +x run_Crotalus_VariantFilter.sh

# Add #!/bin/bash to the top of the run_Crotalus_VariantFilter.sh file

# Run with parallel
parallel --progress --joblog ./log/logfile.variantfilter -j 5 < run_Crotalus_VariantFilter.sh
```


Format `./listVCFHardFilter.txt` with paths to interval HardFilter VCFs.

Generate a sample path with the full list.
``` 
ls -1 -d "$PWD/vcf/"*HardFilter.vcf.gz > ./HardFilterVCFlist.list 

```

Run Picard to merge VCFs.
```
nohup gatk MergeVcfs -I ./HardFilterVCFlist.list -O /media/queen/TombBucket/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.HardFilter.vcf.gz > ./log/Picard_MergeVCFs.HardFilter.log &
```

Remove raw interval VCFs taking up disk space.
```
rm /media/queen/TombBucket/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.raw.*
```

Run bcftools to mask indels/hard filtered sites.
```
nohup bcftools filter --threads 32 -e 'TYPE="indel" || FILTER="QD2" || FILTER="FS60" || FILTER="MQ40" || FILTER="MQRankSum-12.5" || FILTER="ReadPosRankSum-8"' --set-GTs . -O z -o /media/queen/TombBucket/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.HardFilter.recode.vcf.gz /media/queen/TombBucket/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.HardFilter.vcf.gz > ../log/Indelfilter.log & 
nohup bcftools index --threads  32 -t /media/queen/TombBucket/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.HardFilter.recode.vcf.gz > ../log/bcftools_index_hardfilter.log & 
```

###### 2. Filtering female heterozygous sites on the sex chromosomes, we only have the Z at this point.

Females are hemizygous ZW, so should not have heterozygous genotype calls on either the Z or W chromosome.

We'll identify any heterozygous calls in known females then conservatively mask these in all individuals.

####### 1. Format scaffold-assigned BED files for parsing chromosomes, these files are tab-delimited and have 3 columns: the chromosome, the start of the chromosome (1), and the end of the chromosome (length of the chromosome)

```
cut -f1,2 /media/queen/extradrive1/crotalus_pyrrhus_genome/Crotalus_pyrrhus.hap1.final.fasta.fai | awk -v OFS='\t' '{print $0,"1"}' | awk -v OFS='\t' '{print$1,$3,$2}' > crotalus_chrom_regions.bed

## Then we split based on autosomes and Z
# Get the Z
grep "chrZ" crotalus_chrom_regions.bed > crotalus_Zchrom.bed

# Get the autosomes (everything but the Z)
grep -v "chrZ" crotalus_chrom_regions.bed > crotalus_autosomes.bed

# Check to confirm it worked
grep "chrZ" crotalus_autosomes.bed
```

####### 2. Format female individual list

This is based on the sex identification procedure above.

```
# Create a text file with the sexes
cut -d " " -f 2,5 /media/queen/extradrive1/crotalus_genomic_landscape/analysis/sex_identification/runSexIdentification.log > crotalus_samplesex.txt

# Remove the embryos
grep -v "EMB" crotalus_samplesex.txt > crotalus_samplesex_noEMB.txt
mv crotalus_samplesex_noEMB.txt crotalus_samplesex.txt

# Isolate males and females
grep -w "female" crotalus_samplesex.txt | cut -f 1 -d " " > crotalus_females.txt
grep -w "male" crotalus_samplesex.txt | cut -f 1 -d " " > crotalus_males.txt


```

####### 3. Parse autosomes and Z chromosome


```
nohup bcftools view --threads 24 -R crotalus_autosomes.bed -O z -o /media/queen/TombBucket/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.HardFilter.recode.auto.vcf.gz /media/queen/TombBucket/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.HardFilter.recode.vcf.gz > ./log/autosome_subset.log &
nohup bcftools view --threads 16 -R crotalus_Zchrom.bed -O z -o /media/queen/TombBucket/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.HardFilter.recode.chrZ.vcf.gz /media/queen/TombBucket/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.HardFilter.recode.vcf.gz > ./log/zchrom_subset.log &
```

####### 4. Extract biallelic SNPs from the Z chromosome VCF

```
nohup bcftools view --threads 16 -m2 -M2 -U -v snps -O z -o /media/queen/TombBucket/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.HardFilter.recode.chrZ.snps.vcf.gz /media/queen/TombBucket/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.HardFilter.recode.chrZ.vcf.gz > ./log/biallelic_zfilter.log &
```

####### 5. Identify female heterozygous sites using Python script

Run `./sexChrFemaleHeterozygous.py` to extract female heterozygous genotype positions on the sex chromosomes.
```
gunzip /media/queen/TombBucket/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.HardFilter.recode.chrZ.snps.vcf.gz 
nohup python sexChrFemaleHeterozygous_v2.py crotalus_females.txt /media/queen/TombBucket/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.HardFilter.recode.chrZ.snps.vcf ./log/crotalus_genus.allsites.HardFilter.recode.chrZ.snps.FemaleZhetSites.txt > ./log/identifyfemalehetsonZ.log &
```

####### 6. Convert output to BED format and index with GATK

```
conda activate gatk4
awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2}' ./log/crotalus_genus.allsites.HardFilter.recode.chrZ.snps.FemaleZhetSites.txt > ./log/crotalus_genus.allsites.HardFilter.recode.chrZ.snps.FemaleZhetSites.bed
nohup gatk IndexFeatureFile -I ./log/crotalus_genus.allsites.HardFilter.recode.chrZ.snps.FemaleZhetSites.bed > ./log/Index_Zhetsites.log &
```

####### 7. Run GATK VariantFiltration to annotate female heterozygous sites and mask with bcftools

Z chromosome:
```
conda activate gatk4
nohup tabix -p vcf /media/queen/TombBucket/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.HardFilter.recode.chrZ.vcf.gz > ./log/tabix_chrZhardfilteredvcf.log &

nohup gatk VariantFiltration -V /media/queen/TombBucket/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.HardFilter.recode.chrZ.vcf.gz --mask ./log/crotalus_genus.allsites.HardFilter.recode.chrZ.snps.FemaleZhetSites.bed --mask-name ZHET -O /media/queen/TombBucket/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.HardFilter.recode.chrZ.filter.vcf.gz > ./log/GATK_VariantFiltration_crotalus_genus.allsites.HardFilter.recode.chrZ.filter.log &

nohup bcftools filter --threads 24 -e 'FILTER="ZHET"' --set-GTs . -O z -o /media/queen/TombBucket/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.final.chrZ.vcf.gz /media/queen/TombBucket/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.HardFilter.recode.chrZ.filter.vcf.gz > ./log/filterZhets.log & 

nohup tabix -p vcf /media/queen/TombBucket/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.final.chrZ.vcf.gz  > ./log/tabix_finalZ.log &
```

###### 3. Repeat masking

Repeat annotations are in `/media/queen/extradrive1/crotalus_pyrrhus_genome/annotation/repeats/repeatmasker/5_full_mask`. We need to create a bed file that has the coordinates of the repeats.

```
cd /media/queen/extradrive1/crotalus_pyrrhus_genome/annotation/repeats/repeatmasker/5_full_mask

grep -v "^#" Crotalus_pyrrhus.hap1.final.full_mask.final.gff3 > Crotalus_pyrrhus_repeats_tmp.gff
cut -f 1,4,5 -d $'\t' Crotalus_pyrrhus_repeats_tmp.gff > Crotalus_pyrrhus_repeats.bed
rm Crotalus_pyrrhus_repeats_tmp.gff

conda activate base 

bedtools sort -i /media/queen/extradrive1/crotalus_pyrrhus_genome/annotation/repeats/repeatmasker/5_full_mask/Crotalus_pyrrhus_repeats.bed

conda activate gatk4
nohup gatk IndexFeatureFile -I /media/queen/extradrive1/crotalus_pyrrhus_genome/annotation/repeats/repeatmasker/Crotalus_pyrrhus_repeats_sorted.bed > ./log/Index_repeats.log &


```

The `Crotalus_pyrrhus_repeats.bed` contains the repeat locations. We'll use that to annotate the repeats. 

####### 1. Annotate repeats with GATK VariantFiltration

```
conda activate gatk4
cd /media/queen/TombBucket/crotalus_genomic_landscape/
 
tabix -p vcf /media/queen/TombBucket/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.HardFilter.recode.auto.vcf.gz

nohup gatk VariantFiltration -V /media/queen/TombBucket/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.HardFilter.recode.auto.vcf.gz --mask /media/queen/extradrive1/crotalus_pyrrhus_genome/annotation/repeats/repeatmasker/Crotalus_pyrrhus_repeats_sorted.bed --mask-name REP -O /media/queen/TombBucket/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.final.auto.tmp-rep.vcf.gz > ./log/GATK_VariantFiltration_repeat-mask.crotalus_genus.auto.log &

nohup gatk VariantFiltration -V /media/queen/TombBucket/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.final.chrZ.vcf.gz --mask /media/queen/extradrive1/crotalus_pyrrhus_genome/annotation/repeats/repeatmasker/Crotalus_pyrrhus_repeats_sorted.bed --mask-name REP -O /media/queen/TombBucket/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.final.chrZ.tmp-rep.vcf.gz > ./log/GATK_VariantFiltration_repeat-mask.crotalus_genus.chrZ.log &

```

####### 2. Recode repeat annotations as missing genotypes and index VCFs

```
nohup bcftools filter --threads 16 -e 'FILTER="REP"' --set-GTs . -O z -o /media/queen/TombBucket/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.final.auto.vcf.gz /media/queen/TombBucket/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.final.auto.tmp-rep.vcf.gz > ./log/crotalus_auto_repeatmasking.log &

nohup tabix -p vcf /media/queen/TombBucket/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.final.auto.vcf.gz > ./log/crotalus_tabix_finalauto.log &

nohup bcftools filter --threads 16 -e 'FILTER="REP"' --set-GTs . -O z -o /media/queen/TombBucket/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.final.chrZ.vcf.gz /media/queen/TombBucket/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.final.chrZ.tmp-rep.vcf.gz > ./log/crotalus_Z_repeatmasking.log &

nohup tabix -p vcf /media/queen/TombBucket/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.final.chrZ.vcf.gz  > ./log/crotalus_tabix_finalZ.log &
```

###### 4. Extract chromosome-specific all-sites VCFs

We will extract chromosome-specific all-sites vcfs for pixy.

####### 1. Set up environment

```
cd /media/queen/TombBucket/crotalus_genomic_landscape/vcf
mkdir chrom-specific-genus
```

####### 2. Parse chromosome-specific autosome VCFs

```
nohup bcftools index --threads 24 -s /media/queen/TombBucket/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.final.auto.vcf.gz | cut -f 1 | while read C; do bcftools view --threads 24 -O z -o /media/queen/TombBucket/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.final.${C}.vcf.gz /media/queen/TombBucket/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.final.auto.vcf.gz "${C}" ; done > ./log/autosome_split_bychrom.log &

```

####### 3. Parse Z chromosome

```
cp ../crotalus_genus.allsites.final.chrZ.vcf.gz ./
```

####### 4. Index VCFs

```
nohup bash -c 'for i in /media/queen/TombBucket/crotalus_genomic_landscape/vcf/chrom-specific-genus/*noCV13*.vcf.gz; do tabix -p vcf $i; done' > ./log/tabix_final_allsitesvcfs.log &
```

Make a directory to store the unplaced scaffolds, move them there.

```

mkdir /media/queen/TombBucket/crotalus_genomic_landscape/vcf/chrom-specific-genus/unplaced_scaf_vcfs

mv *.scaffold* ./unplaced_scaf_vcfs

```

###### 5. Extract biallelic SNPs

####### 1. Extract SNPs

```
cd /media/queen/TombBucket/crotalus_genomic_landscape

nohup bcftools view --threads 16 -m2 -M2 -U -v snps -O z -o /media/queen/TombBucket/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.final.auto.snps.vcf.gz /media/queen/TombBucket/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.final.auto.vcf.gz > ./log/extract_biallelic_auto.log &

nohup bcftools view --threads 16 -m2 -M2 -U -v snps -O z -o /media/queen/TombBucket/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.final.chrZ.snps.vcf.gz /media/queen/TombBucket/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.final.chrZ.vcf.gz > ./log/extract_biallelic_z.log &

```

####### 2. Concatenate autosome and Z chromosome SNP VCFs

```

bcftools concat -O z -o /media/queen/TombBucket/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.final.auto+chrZ.snps.vcf.gz /media/queen/TombBucket/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.final.auto.snps.vcf.gz /media/queen/TombBucket/crotalus_genomic_landscape/vcf/crotalus_genus.allsites.final.chrZ.snps.vcf.gz

```

## Phylogeny, population structure, and demography

We estimated the phylogeny of the Speckled and Western Rattlesnake species complexes using multiple approaches that are listed below, the pipeline to perform this set of analyses is called `phylogeny_structure.md`.

*Software for Phylogeny, population structure, and demography*:

- Concatenated maximum likelihood with [RAxML](https://github.com/amkozlov/raxml-ng)
- Coalescent species tree inference with [SVDquartets](https://www.asc.ohio-state.edu/kubatko.2/software/SVDquartets/)
- Phylogenetic network analysis with [SplitsTree](https://uni-tuebingen.de/en/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/splitstree/)
- Estimated divergence times with [treePL](https://github.com/blackrim/treePL)
- Inferred changes in effective population size with [smc++](https://github.com/popgenmethods/smcpp)
- Inferred genetic structure using principal component analysis using [SNPrelate](https://www.bioconductor.org/packages/release/bioc/html/SNPRelate.html)
- Inferred genetic structure using model-based ancestry using [ADMIXTURE](https://dalexander.github.io/admixture/)
- Processed and visualized results with [R](https://www.r-project.org/)

## Recombination rate and recombination hotspot identification

We estimated recombination maps for *C. pyrrhus*, *C. stephensi*, *C. viridis*, *C. concolor*, and *C. helleri*. We then followed previous approaches developed by [Schield et al. (2020)](https://academic.oup.com/mbe/article/37/5/1272/5700722) and [Hoge et al. (2024)](https://www.science.org/doi/10.1126/science.adj7026) to identify recombination hotspots and coldspots.

*Software for Recombination rate and recombination hotspot identification*:

- Inferred changes in effective population size with [smc++](https://github.com/popgenmethods/smcpp)
- Estimated recombination maps with [pyrho](https://github.com/popgenmethods/pyrho)
- Calculate mean recombination rate in sliding windows with [bedtools](https://bedtools.readthedocs.io/en/latest/)
- Processed and visualized results with [R](https://www.r-project.org/)

## Calculation of introgression and divergence statistics

We estimated the admixture proportion (fd statistic; [Martin et al. (2015)](https://academic.oup.com/mbe/article/32/1/244/2925550?guestAccessKey=) to examine genome-wide patters of introgression. We also estimated absolute (dxy) and relative divergence (Fst) using [pixy](https://pixy.readthedocs.io/en/latest/).

*Software for Calculation of introgression and divergence statistics*:

- Estimated fd with [ABBABABAwindows.py](https://github.com/simonhmartin/genomics_general)
- Estimated divergence measures with [pixy](https://pixy.readthedocs.io/en/latest/)
- Processed and visualized results with [R](https://www.r-project.org/)

## Distinguishing introgression from incomplete lineage sorting

We tested whether our observed patterns may have been produced by incomplete lineage sorting, rather than introgression, using analyses of linkage disequilibrium and absolute sequence divergence (dxy). 

*Software for Distinguishing introgression from incomplete lineage sorting*:

- Estimated linkage disequilibrium with [vcftools](https://vcftools.sourceforge.net/)
- Estimated divergence measures with [pixy](https://pixy.readthedocs.io/en/latest/)
- Processed and visualized results with [R](https://www.r-project.org/)

## Evolutionary and ecological divergence

We calculated the evolutionary and ecological divergence to test the relationship between introgression and types of divergence.

*Software for Evolutionary and ecological divergence*:

- Estimating evolutionary divergence (branch lengths) with [RAxML](https://github.com/amkozlov/raxml-ng)
- Estimating ecological divergence, processing, and visualizing results with [R](https://www.r-project.org/)

## Relationships between introgression, recombination, and divergent selection

We tested the relationship between introgression, recombination rate, exon density, and relative differentiation (Fst). We used the output from previous steps (see above) in these tests.

*Software for Relationships between introgression, recombination, and divergent selection*:

- We tested these relationships in [R](https://www.r-project.org/)

## Identification of species barriers and fine-scale variation in introgression

We identified outlier barrier windows and tested for differential signals of introgression in various genomic regions (e.g., venom, barrier windows) using fd statistics.

*Software for Identification of species barriers and fine-scale variation in introgression*:

- We tested these relationships in [R](https://www.r-project.org/)

