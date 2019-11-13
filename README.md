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
6. [E. coli assembly correctness](#ecoli.correct)
7. [E. coli plasmid analysis](#ecoli.plasmid)



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
**Polish, circularize, and rotate**  
pilon_iter.sh canu/assembly.contigs.fasta illuminaPE_R1.fastq.gz illuminaPE_R2.fastq.gz canu/assembly.trimmedReads.fasta  
circlator minimus2 pilon5.fasta circularise.fasta  
circlator fixstart --genes_fa ecoli.dnaA.DNA.fasta circularise.fasta rotated.fasta  
**Rename FASTA**  
for f in *_contigs.fasta; do awk '/^>/{print ">ecoli_contig" ++i; next}{print}' < $f > ${f%_c*}_contigs_rn.fasta; done

### E. coli genome assembly using Unicycler <a name="ecoli.uni"></a>
**Long read assembly**  
unicycler --mode normal --start_genes ecoli.dnaA.protein.fasta --long long.reads.fastq  -o output.dir -t 32  
**Hybrid assembly**  
unicycler --mode normal --start_genes ecoli.dnaA.protein.fasta --long long.reads.fastq  --short1 short.reads.R1.fastq --short2 short.reads.R2.fastq -o output.dir -t 32  
**Rename FASTA**  
for f in \*contigs.fasta; do awk '/^>/{print ">ecoli_contig" ++i; next}{print}' < $f > ${f%\_c\*})\_contigs_rn.fasta; done

### E. coli BUSCO <a name="ecoli.busco"></a>  
python run_BUSCO.py -f -c 8 -t /local/scratch/etvedte/tmp -i assembly.fasta -o busco_output_dir -l bacteria_odb9 -m geno  

### E. coli assembly correctness <a name="ecoli.correct"></a>  
**De novo assembly of Illumina reads using AbySS**  
abyss-pe -j 16 name=ecoli_illumina_abyss k=116 in='short.reads.R1.fastq short.reads.R2.fastq'  
**De novo assembly of Illumina reads using velvet**  
velveth /velveth/output/dir/ 115 -short -fastq short.reads.R1.fastq  
velvetg /velveth/output/dir/  
**Identify trusted contigs, 100% match between ABySS and velvet**  
mummer -mum -b -s -n -l 5000 ecoli.abyss.contigs.fasta ecoli.velvet.contigs.fasta | grep -i -P "[acgt]{5000,}" | awk '{name += 1; print ">"name"\n"substr(toupper($0), 100, length($0)-200)}' > ecoli.abyss+velvet.trusted.contigs.fasta  
**Build BLAST database for long read assemblies**  
for f in \*contigs_rn.fasta; do makeblastdb -dbtype nucl -in $f -out ${f%\_r\*})\_rn.blastdb -parse_seqids; done
**Estimate chimeric reads using Alvis**  


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
