# Genomic Epidemiology of the 2010 Haiti Cholera Outbreak
 
## Table of Contents
 
- [Project Background](#project-background)
- [Design Rationale and Summary](#design-rationale-and-summary)
- [DAG](#dag)
- [Data Explanation](#data-explanation)
- [Core Dependencies](#core-dependencies)
- [Instructions](#instructions)
- [Expected Results](#expected-results)
- [Directory Structure Overview](#directory-structure-overview)
- [Gen AI Use Statement](#gen-ai-use-statement)
## Project Background
 
For much of human history, epidemics have been poorly understood. Little was known about how they spread or which cases could be confidently linked. Next-generation sequencing has changed that. Within a seemingly large outbreak, genomic data can now reveal distinct clusters of related strains, and the nucleotide-level differences between them can be used to infer who infected whom and how the pathogen moved through a population.

So how do we define an infectious disease outbreak genomically? How can we reconstruct transmission dynamics? Genomic epidemiology seeks to answer these questions using bioinformatic techniques.

While these techniques can be applied to any infectious disease, this project focuses on *Vibrio cholerae*, a water-borne pathogen.
 
> **Fun fact:** *V. cholerae* was the subject of the first effective epidemiological map, constructed by John Snow in 1854 and used to trace the source of an outbreak to a [contaminated well pump in London](https://en.wikipedia.org/wiki/1854_Broad_Street_cholera_outbreak).
>
> **Not-so-fun fact:** *V. cholerae* has caused roughly seven major pandemics since 1817 (the seventh is still ongoing), and causes an estimated 95,000 deaths annually worldwide. [(Ali et al.)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4455997/)
 
In January 2010, Haiti was hit by a devastating earthquake and received international aid in the form of food, water, and aid workers. Nine months later, in October 2010, Haiti confirmed a cholera outbreak. Early studies suggested the outbreak strain was related to strains from South Asia, Bangladesh, or Mozambique, but a genomic epidemiological investigation by [Hendriksen et al.](https://journals.asm.org/doi/10.1128/mbio.00157-11) identified Nepalese aid workers as the likely point of entry. This genomic evidence, together with epidemiological evidence such as the location of the Nepalese camp relative to the first reported cases, strongly supported Nepal as the source.

The investigation by Hendriksen et al. was limited by the technology of its time. The authors aligned sample reads to a known *V. cholerae* reference genome and constructed a maximum parsimony tree. We now know that significant differences between bacterial samples exist in the accessory genome, and these can only be captured through *de novo* assembly, not reference alignment. In addition, increased computational resources have made higher-resolution phylogenetic methods such as maximum likelihood practical, and these have largely replaced maximum parsimony.

**If we redo the 2010 analysis using modern bioinformatic tools and include more countries, can we still expect the Nepalese samples to be closely related to the Haitian ones?**
 
## Design Rationale and Summary
 
> **Note:** This pipeline was designed with historical data in mind (paired-end short-read sequencing). Modern day standard is use both short and long reads in a "hybrid" assembly to achieve the maximum possible accuracy. (See [Ryan Wick's blog](https://rrwick.github.io/) for more details!). Short read bacterial sequencing is still done today, but as long read sequencing becomes more common, this design choice limits the application of the pipeline to modern day outbreaks.
 
### Input
 
- `samplesheet.csv`: sample country and date metadata, SRR accession, and (optionally) `fastq1` and `fastq2` file addresses
- *V. cholerae* reference genome FASTA and GFF files (used for QUAST and Parsnp, but **not** for reference alignment).
### Tools

All tools were available as modules under the [nf-core project](https://nf-co.re/):
 
| Step | Tool | Rationale |
|------|------|-----------|
| Raw read QC | [fastp](https://nf-co.re/modules/fastp) | Chosen for its versatility. A single tool that covers the QC reporting of FastQC and the adapter and quality trimming of Trimmomatic, with additional features such as automatic adapter detection and fast multithreaded processing. |
| Short-read *de novo* assembly | [SKESA](https://nf-co.re/modules/skesa/) | Strategic K-mer Extension for Scrupulous Assemblies. Alongside SPAdes, it is one of the most widely used short-read assemblers. Compared with SPAdes, SKESA is more accurate but produces more fragmented assemblies. Fragmentation loses structural information, but this analysis depends on nucleotide-level differences rather than genome structure, so the trade-off is acceptable. |
| Assembly completeness and contamination QC | [BUSCO](https://nf-co.re/modules/busco_busco) | Uses Benchmarking Universal Single-Copy Orthologs to assess completeness, with duplicated genes indicating possible contamination. The alternative, CheckM, selects its marker gene set by placing each sample taxonomically. Because all samples here are *V. cholerae*, a single lineage-specific BUSCO dataset can be used, which saves space. |
| General assembly QC | [QUAST](https://nf-co.re/modules/quast) | Widely used for reporting assembly quality metrics (N50, total length, number of contigs) and for comparing assemblies against a reference. |
| Core genome phylogeny | [Parsnp](https://nf-co.re/modules/parsnp/) | Performs multiple sequence alignment to identify regions present in all samples, then builds a maximum likelihood tree from the core SNPs (using FastTree). To my knowledge, there is not a ton of differences in phylogenetic tree construction tools. Parsnp happens to be the one I am most familiar with. |
| Resistance gene profiling | [CARD RGI](https://nf-co.re/modules/rgi_main) | Uses the Comprehensive Antibiotic Resistance Database (CARD) and its Resistance Gene Identifier (RGI) to extract a resistance gene profile for each sample. ABRicate, which queries several AMR databases, would normally be preferred, but was not feasible given space and runtime limits. The decision to use CARD over other databases was completely informed and I totally did not just choose it because it's Canadian. |
| Visualization | [ggplot2](https://ggplot2.tidyverse.org/), [ggtree](https://bioconductor.org/packages/ggtree) (R) | Standard libraries for plotting phylogenetic trees and heatmaps. |
 
### Output (found in results/r_plots after runnning)
 
- Phylogenetic tree overlaid with country and sampling date metadata
- Antimicrobial resistance gene presence/absence heatmap
## Directed Acyclic Graph
 
![Pipeline DAG](assets/pipeline.svg)
 
## Data Explanation
 
Rather than using the Hendriksen et al. data directly, global samples were taken from [Mutreja et al.](https://www.nature.com/articles/nature10392), which contains 154 *V. cholerae* isolates collected worldwide. To construct a toy example dataset, single-end samples were filtered out, and one representative sample per country was selected. A Haitian sample and Nepalaese sample deemed closely related by Hendriksen et al. were selected as representatives.

The reference genome used was *V. cholerae* O1 biovar El Tor strain N16961 (isolated in Bangladesh), the El Tor biovar being responsible for the ongoing seventh cholera pandemic.

> **Note:** The reference is **not** used for read alignment or variant calling. Every sample is assembled *de novo*, and the reference serves only as the anchor for Parsnp's core-genome alignment (and as the comparison genome for QUAST).
 
## Core Dependencies
 
This pipeline was tested on ARM64 (an Apple Silicon M1 Mac), but it is expected to work on AMD64 as well, since most nf-core containers are built for that architecture. On ARM64, most modules run by emulating their Docker images. The exceptions are SKESA and the R environment, both of which use conda instead.

- [Docker](https://www.docker.com/)
- [Nextflow](https://www.nextflow.io/)
- [Conda](https://docs.conda.io/)

## Instructions
 
In the project directory, run:
 
```bash
nextflow run .
```
 
That's all! The pipeline will download the necessary container, samples, and databases during its first run.
 
## Expected Results

(If using default samplesheet.csv)
 
![Phylogenetic Tree](assets/plot_tree.png)

![AMR Heatmap](assets/plot_AMR.png)
 
## Directory Structure Overview
 
```
BIOF_501A_project/
├── assets/                  # Static pipeline files (e.g. samplesheet.csv, reference genome FASTA, reference GFF)
├── conf/                    # Configuration files
│   ├── base.config          # Default resource requirements (CPUs, memory, time)
│   └── modules.config       # Per-module options (container, arguments, output directories)
├── databases/               # Reference databases used by pipeline modules
├── inputs/                  # Input data (raw sequencing data per sample)
├── modules/                 # Process definitions for each tool, pulled from nf-core (e.g. SKESA, Parsnp)
├── results/                 # Pipeline outputs (generated, per module)
├── subworkflows/            # Subworkflow definitions pulled from nf-core (e.g. fastq_download)
├── work/                    # Nextflow working directory with intermediate files (generated, per run)
├── main.nf                  # Main workflow script
├── modules.json             # nf-core module/subworkflow version tracking
├── nextflow.config          # Global pipeline configuration
├── pipeline.mmd             # Mermaid source for the pipeline diagram
├── pipeline.svg             # Rendered pipeline diagram
└── README.md                # Project documentation
```
 
## Gen AI Use Statement
 
- All analysis, ideas, and tool choices were made by Ryan Choi.
- Code was also originally written by Ryan Choi, but tidied up with Claude afterwards.
 
