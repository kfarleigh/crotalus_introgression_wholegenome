# crotalus_introgression_wholegenome

*Last updated May 8th, 2026*

**WARNING, LINKS MAY NOT WORK UNTIL AFTER THE MANUSCRIPT HAS BEEN PUBLISHED.**

I am in the process of uploading scripts to this repository to ensure that they are properly annotated, but please reach out if you have any questions.

This repository contains the computational workflow and scripts for [Farleigh et al., (*accepted*)](). Please email Keaka Farleigh (keakafarleigh@virginia.edu; keakafarleigh@gmail.com) if you have any questions. 

## Citation

The citation will be here.

## Genomic data

The genomic data for this project are deposited as parts of BioProject [PRJNA1454467](), [PRJNA593834](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA593834), [PRJNA1150930](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1150930). The genome for this project is deposited at [PRJNA1460337]().

## Structure of this repository.

This repository contains the computational workflow for [Farleigh et al., (*accepted*)](). This ReadMe will tell you which scripts are associated with each analysis. I have separated this ReadMe into each of the methods sections listed in the [supplemental information](). Each script referenced herein is located in the `scripts/` directory.

### Genome assembly and annotation

We assigned chromosome names using a synteny-based approach as implemented in [mashmap](https://github.com/marbl/MashMap). The script to perform this is called `chromosome_assignment.md`.

### Whole genome sequencing and variant calling

We sequenced newly generated libraries on Illumina NovaSeq 6000 lanes to generate 150 bp paired-end reads. We then trimmed our sequence data with [Trimmomatic](https://github.com/usadellab/trimmomatic). We then followed [GATK's best practices workflow](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows) to call and filter SNPs. The pipeline to perform this set of analyses is called `variant_calling.md`.

### Phylogeny, population structure, and demography

We estimated the phylogeny of the Speckled and Western Rattlesnake species complexes using multiple approaches that are listed below, the pipeline to perform this set of analyses is called `phylogeny_structure.md`.

- Concatenated maximum likelihood with [RAxML](https://github.com/amkozlov/raxml-ng)
- Coalescent species tree inference with [SVDquartets](https://www.asc.ohio-state.edu/kubatko.2/software/SVDquartets/)
- Phylogenetic network analysis with [SplitsTree](https://uni-tuebingen.de/en/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/splitstree/)
- Estimated divergence times with [treePL](https://github.com/blackrim/treePL)
- Inferred genetic structure using principal component analysis using [SNPrelate](https://www.bioconductor.org/packages/release/bioc/html/SNPRelate.html)
- Inferred genetic structure using model-based ancestry using [ADMIXTURE](https://dalexander.github.io/admixture/)

### Recombination rate and recombination hotspot identification

### Calculation of introgression and divergence statistics

### Distinguishing introgression from incomplete lineage sorting

### Evolutionary and ecological divergence

### Relationships between introgression, recombination, and divergent selection

### Identification of species barriers and fine-scale variation in introgression



