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
7. [E. coli chimeric reads assessment](#ecoli.chimera)  
8. [E. coli plasmid analysis](#ecoli.plasmid)  
9. [D. ananassae genome assembly using Canu](#dana.canu)  
10. [D. ananassae genome assessment](#dana.eval)  
11. [D. ananassae BUSCO](#dana.busco)  
11. [D. ananassae chromosome map](#dana.chrom)



### Prepare sequencing files for assembly <a name="ecoli.prep"></a>
**MinION LIG**  
zcat minion.LIG.raw.fastq.gz | NanoLyse | gzip > minion.LIG.filter.fastq.gz
seqkit sample -s 13 -j 16 -p 0.17218 -o minion.LIGsample.fastq minion.LIG.filter.fastq.gz  
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
**Identify trusted contigs, >5kb and 100% match between ABySS and velvet**  
mummer -mum -b -s -n -l 5000 ecoli.abyss.contigs.fasta ecoli.velvet.contigs.fasta | grep -i -P "[acgt]{5000,}" | awk '{name += 1; print ">"name"\n"substr(toupper($0), 100, length($0)-200)}' > ecoli.abyss+velvet.trusted.contigs.fasta  
**Build BLAST database for long read assemblies**  
mkdir db  
mv \*\_contigs_rn.fasta db/  
cd db  
for f in \*fasta; do makeblastdb -dbtype nucl -parse_seqids -in $f ; done  
**Retrieve the best BLAST hit for each trusted contig query**  
for f in db/\*fasta; do blastn -db $f -query /local/projects-t3/RDBKO/ecoli.abyss/ecoli.abyss+velvet.5kb+.trusted.contigs.fasta -outfmt 6 | sort -nk1,1 -k12,12gr -k11,11g -k3,3gr | sort -u -nk1,1 > ${f%\_r\*}\_trusted.blast_hits; done  
**Determine percentage correctness, aligned bases, mismatches, gap opens for long read assemblies**  
for f in \*\_trusted_blast_hits; do python3 get_error_rate.py $f > ${f%\_t\*}\_err_data; done  
**Generate summary file**  
for f in \*blast\_hits; do echo ${f%\_t\*} >> ecoli.correctness.names.txt  
for f in \*err_data; do cat $f >> ecoli.correctness.data.txt; done 
paste ecoli.correctness.names.txt ecoli.correctness.data.txt > ecoli.correctness.summary.txt  

### E. coli chimeric reads assessment <a name="ecoli.chimera"></a>  
**Map reads to Unicycler consensus using minimap2. Second round of mapping to cut-and-paste E.coli genome sequence**  
*MinION*  
minimap2 -x map-ont -t 8 ecoli.unicycler.consensus.fasta minion.reads.fastq > mapped.paf  
minimap2 -x map-ont -t 8 ecoli.unicycler.consensus.cut.fasta minion.reads.fastq > mapped.cut.paf  
*PacBio*  
minimap2 -x map-pb -t 8 ecoli.unicycler.consensus.fasta pb.reads.fastq > mapped.paf  
minimap2 -x map-pb -t 8 ecoli.unicycler.consensus.cut.fasta pb.reads.fastq > mapped.cut.paf  
**Filter to retain reads mapped to E. coli genome**  
*FASTA header is ecoli.genome and ecoli.genome.cut, respectively*  
grep ecoli.genome mapped.paf > genome.paf  
grep ecoli.genome.cut mapped.cut.paf > genome.cut.paf  
**Predict chimeras using Alvis**  
alvis/dist/Alvis.jar -type contigAlignment -inputfmt paf -outputfmt svg -chimeras -printChimeras -minChimeraCoveragePC 90 -minChimeraAlignmentPC 10 -in genome.paf -outdir /path/to/outdir/ -out out.prefix  
mv chimeras.txt genome.chimeras.txt  
alvis/dist/Alvis.jar -type contigAlignment -inputfmt paf -outputfmt svg -chimeras -printChimeras -minChimeraCoveragePC 90 -minChimeraAlignmentPC 10 -in genome.cut.paf -outdir /path/to/outdir/ -out out.prefix  
mv chimeras.txt genome.chimeras.cut.txt  
**Determine estimate for percentage chimeras**
*Identify reads assigned as chimeras in both original and cut E.coli genome*
cat genome.chimeras.txt genome.cut.chimeras.txt | awk '{print $1}' - | sort -n | uniq -d > chimeras.candidates.txt  
**Re-map, filter reads to include only chimeras, manual inspection**  
*Map MinION*  
minimap2 -ax map-ont -t 8 ecoli.unicycler.consensus.fasta minion.reads.fastq | samtools sort -o sorted.bam  
*Map PacBio*  
minimap2 -ax map-pb -t 8 ecoli.unicycler.consensus.fasta pb.reads.fastq | samtools sort -o sorted.bam  
*Filter reads*  
java -Xmx10g -jar picard-tools-2.5.0/picard.jar FilterSamReads I=sorted.bam O=chimeras.bam READ_LIST_FILE=chimeras_candidates.txt FILTER=includeReadList
samtools index chimeras.bam  
**Determine estimate for percentage chimeras**  
wc -l chimeras.candidates.txt  
*Total primary reads mapped to genome*  
samtools view sorted.bam -F 4 -F 256 -F 1024 -F 2048 ecoli.genome -c  
percentage chimeras = chimeras candidates / primary reads

### E. coli plasmid analysis <a name="ecoli.plasmid"></a>  
**Produce depth counts for primary reads mapped to genome, pMAR2, p5217**  
samtools view sorted.bam -F 4 -F 256 -F 1024 -F 2048 -bho primary.bam  
samtools depth -aa -m 100000000 primary.bam > primary.depth.txt  
**Count total depth across genome, pMAR2, p5217**  
grep ecoli.genome primary.depth.txt | awk '{total = total + $3}END{print "Total genome depth = "total}' -  
grep pMAR2 primary.depth.txt | awk '{total = total + $3}END{print "Total pMAR2 depth = "total}' -  
grep p5217 primary.depth.txt | awk '{total = total + $3}END{print "Total p5217 depth = "total}' -

### D.ananassae genome assembly using Canu <a name="dana.canu"></a>
**MinION**  
canu -p output.prefix -d output.dir genomeSize=240m gridEngineThreadsOption="-pe thread THREADS" gridEngineMemoryOption="-l mem_free=MEMORY" gridOptions="-P jdhotopp-lab -q threaded.q" -nanopore-raw raw.nanopore.reads.fastq  
**PacBio**  
canu -p output.prefix -d output.dir genomeSize=240m gridEngineThreadsOption="-pe thread THREADS" gridEngineMemoryOption="-l mem_free=MEMORY" gridOptions="-P jdhotopp-lab -q threaded.q" -pacbio-raw raw.pacbio.reads.fastq  
**Hybrid**  
canu -p output.prefix -d output.dir genomeSize=240m gridEngineThreadsOption="-pe thread THREADS" gridEngineMemoryOption="-l mem_free=MEMORY" gridOptions="-P jdhotopp-lab -q threaded.q" -nanopore-raw raw.nanopore.reads.fastq -pacbio-raw raw.pacbio.reads.fastq    
**Genome polishing**  
*Long read*  
pilon_iter.sh canu/assembly.contigs.fasta canu/assembly.trimmedReads.fasta  
*Long + short read*  
pilon_iter.sh canu/assembly.contigs.fasta illuminaPE_R1.fastq.gz illuminaPE_R2.fastq.gz canu/assembly.trimmedReads.fasta  
*PacBio HiFi*  
pilon_iter.sh canu/assembly.contigs.fasta pb.HiFi.ccs.fastq.gz

### D. ananassae BUSCO <a name="dana.busco"></a>  
python run_BUSCO.py -f -c 8 -t /local/scratch/etvedte/tmp -i assembly.fasta -o busco_output_dir -l metazoa_odb9 -m geno  

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
