# Comparison of long read sequencing technologies in resolving bacteria and fly genomes

Eric S. Tvedte

2019-11-13

The repository contains Supplementary Data for the manuscript, including Tables, Figures, and Files.

## Table of Contents
1. [Prepare sequencing files for assembly](#ecoli.prep)
2. [Determine read length distributions](#ecoli.read)
3. [E. coli genome assembly using Canu](#ecoli.canu)
4. [E. coli genome assembly using Unicycler](#ecoli.uni)
5. [E. coli BUSCO](#ecoli.busco)
4. [E. coli assembly correctness](#ecoli.correct)
4. [E. coli plasmid analysis](#ecoli.plasmid)



### Prepare sequencing files for assembly <a name="ecoli.prep"></a>
**MinION LIG**  
zcat minion.LIG.raw.fastq.gz | NanoLyse | gzip > minion.LIG.filter.fastq.gz
seqkit sample -s 13 -j 16 -p 0.17218 -o minion.LIG.sample.fastq minion.LIG.filter.fastq.gz  
**PacBio Sequel II**  
seqkit sample -s 13 -j 16 -p 0.0212566 -o pbSequelII.sample.fastq pbSequelII.raw.fastq.gz

### Determine read length distributions <a name="ecoli.read"></a>
**MinION**  
bbtools/readlength.sh in=raw.reads.fastq.gz out=raw.reads.hist.out bin=1000 max=1000000 qin=33
**PacBio**  
bbtools/readlength.sh in=raw.reads.fastq.gz out=raw.reads.hist.out bin=1000 max=1000000
**Generate histograms**  
read_histograms.Rmd  

### E. coli genome assembly using Canu <a name="ecoli.canu"></a>
**MinION**  
canu -p output.prefix -d output.dir genomeSize=4.6m gridEngineThreadsOption="-pe thread THREADS" gridEngineMemoryOption="-l mem_free=MEMORY" gridOptions="-P jdhotopp-lab -q threaded.q" -nanopore-raw raw.nanopore.reads.fastq
**PacBio**  
canu -p output.prefix -d output.dir genomeSize=4.6m gridEngineThreadsOption="-pe thread THREADS" gridEngineMemoryOption="-l mem_free=MEMORY" gridOptions="-P jdhotopp-lab -q threaded.q" -pacbio-raw raw.pacbio.reads.fastq  
**polish, circularize, and rotate**  
pilon_iter.sh canu/assembly.contigs.fasta illuminaPE_R1.fastq.gz illuminaPE_R2.fastq.gz canu/assembly.trimmedReads.fasta  
circlator minimus2 pilon5.fasta circularise.fasta  
circlator fixstart --genes_fa ecoli.dnaA.DNA.fasta circularise.fasta rotated.fasta  
**rename FASTA**  
for f in *_contigs.fasta; do awk '/^>/{print ">ecoli_contig" ++i; next}{print}' < $f > ${f%_c*}_contigs_rn.fasta; done

### E. coli genome assembly using Unicycler <a name="ecoli.uni"></a>

## System requirements

R scripts were run using Windows 10 x64 with RStudio v1.1.463 using this R session:
```
R version 3.5.1 (2018-07-02)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 18362)
```
## Installation
Requirements for installed packages are available in Supplementary Files. 

## 1. supplementary_figures
This folder contains Supplementary Figures reported in the manuscript, Figure S1 - SX

## 2. supplementary_tables
This folder contains Supplementary Tables reported in the manuscript, Table S1 - SX

## 3. scripts
This folder contains Rmd scripts used for data analysis:

FileSX:

Sample input files are provided in sample_data_files. Input and output file paths are hard-coded in the scripts, change these to run the scripts on your local system.

## 4. htmls

This folder contains output html files generated from Rmd files using the R package knitr.

## 5. sample_data_files
This folder contains sample input files for data analysis. 

FileSX: 
