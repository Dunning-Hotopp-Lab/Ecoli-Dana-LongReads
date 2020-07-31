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
9. [Ecoli DNA modification analysis](#ecoli.dna.mod)  
10. [D. ananassae genome assembly using Canu](#dana.canu)  
11. [D. ananassae genome assessment](#dana.eval)  
12. [Alignment of contigs to D. ananassae polytene map](#dana.map)


### Prepare sequencing files for assembly <a name="ecoli.prep"></a>
**Remove DNA Control Sequence from ONT LIG reads**
```
zcat ont.LIG.raw.fastq.gz | NanoLyse | gzip > ont.LIG.filter.fastq.gz
```

### Determine read length distributions <a name="ecoli.read"></a>
**ONT reads** 
```
bbtools/readlength.sh in=ont.reads.fastq.gz out=ont.reads.hist.out bin=1000 max=1000000 qin=33 
```
**PacBio reads**
```
bbtools/readlength.sh in=pb.reads.fastq.gz out=pb.reads.hist.out bin=1000 max=1000000 
```
**Generate histograms**
```
grep -v '#' reads.hist.out > reads.hist.final.out
scripts/Ecoli_Long_Read_Analysis.Rmd
```
**Subset ONT LIG and PacBio Sequel II reads**
```
seqkit sample -s 13 -j 16 -p 0.17218 -o ont.LIG.sample.fastq ont.LIG.filter.fastq.gz
seqkit sample -s 13 -j 16 -p 0.0212566 -o pb.SQII.sample.fastq pbSQII.raw.fastq.gz
```
### E. coli genome assembly using Canu <a name="ecoli.canu"></a>
**ONT assembly**  
```
canu -p output.prefix -d output.dir genomeSize=4.6m gridEngineThreadsOption="-pe thread THREADS" gridEngineMemoryOption="-l mem_free=MEMORY" gridOptions="-P jdhotopp-lab -q threaded.q" -nanopore-raw nanopore.reads.fastq 
```
**PacBio** 
```
canu -p output.prefix -d output.dir genomeSize=4.6m gridEngineThreadsOption="-pe thread THREADS" gridEngineMemoryOption="-l mem_free=MEMORY" gridOptions="-P jdhotopp-lab -q threaded.q" -pacbio-raw pacbio.reads.fastq
```
**Polish, circularize, and rotate**  
```
Modify pilon_iter.sh for polishing with short reads/long reads/both
scripts/pilon_iter.sh assembly.fasta illumina.reads.R1.fastq illumina.reads.R2.fastq canu/assembly.trimmedReads.fasta  
circlator minimus2 pilon5.fasta circularise.fasta  
circlator fixstart --genes_fa ecoli.dnaA.fasta circularise.fasta rotated.fasta  
```

### E. coli genome assembly using Unicycler <a name="ecoli.uni"></a>
**Long read assembly**  
```
unicycler --mode normal --start_genes ecoli.dnaA.protein.fasta --long long.reads.fastq  -o output.dir -t 32  
```
**Hybrid assembly**  
```
unicycler --mode normal --start_genes ecoli.dnaA.protein.fasta --long long.reads.fastq  --short1 illumina.reads.R1.fastq --short2 illumina.reads.R2.fastq -o output.dir -t 32  
```

### E. coli BUSCO <a name="ecoli.busco"></a>  
```
python run_BUSCO.py -f -c 8 -t /local/scratch/etvedte/tmp -i assembly.fasta -o busco_output_dir -l bacteria_odb9 -m geno  
```

### E. coli assembly correctness <a name="ecoli.correct"></a>  
**De novo assembly of Illumina reads using AbySS**  
```
abyss-pe -j 16 name=ecoli_illumina_abyss k=116 in='illumina.reads.R1.fastq illumina.reads.R2.fastq'  
```
**De novo assembly of Illumina reads using velvet**  
```
velveth /velveth/output/dir/ 115 -short -fastq illumina.reads.R1.fastq  
velvetg /velveth/output/dir/  
```
**Identify trusted contigs, >5kb and 100% match between ABySS and velvet** 
```
mummer -mum -b -s -n -l 5000 ecoli.abyss.contigs.fasta ecoli.velvet.contigs.fasta | grep -i -P "[acgt]{5000,}" | awk '{name += 1; print ">"name"\n"substr(toupper($0), 100, length($0)-200)}' > ecoli.trusted.contigs.fasta  
```
**Build BLAST database for long read assemblies**  
```
mkdir db  
mv \*contigs.fasta db/  
cd db  
for f in \*fasta; do makeblastdb -dbtype nucl -parse_seqids -in $f -out $f; done  
```
**Retrieve the best BLAST hit for each trusted contig query**  
```
for f in db/\*fasta; do blastn -db $f -query /local/projects-t3/RDBKO/ecoli.abyss/ecoli.trusted.contigs.fasta -outfmt 6 | sort -nk1,1 -k12,12gr -k11,11g -k3,3gr | sort -u -nk1,1 > ${f%\_r\*}\_trusted.blast_hits; done  
```
**Determine percentage correctness, aligned bases, mismatches, gap opens for long read assemblies**  
```
for f in \*\_trusted_blast_hits; do python3 get_error_rate.py $f > ${f%\_t\*}\_err_data; done  
```
**Generate summary file**  
```
for f in \*blast\_hits; do echo ${f%\_t\*} >> ecoli.correctness.names.txt  
for f in \*err_data; do cat $f >> ecoli.correctness.data.txt; done 
paste ecoli.correctness.names.txt ecoli.correctness.data.txt > ecoli.correctness.summary.txt  
```

### E. coli chimeric reads assessment <a name="ecoli.chimera"></a>  
**Map reads to Unicycler consensus using minimap2. Second round of mapping to cut-and-paste E.coli genome sequence**  
*ONT reads*  
```
minimap2 -x map-ont -t 8 ecoli.unicycler.consensus.fasta minion.reads.fastq > mapped.paf  
minimap2 -x map-ont -t 8 ecoli.unicycler.consensus.cut.fasta minion.reads.fastq > mapped.cut.paf  
```
*PacBio reads*  
```
minimap2 -x map-pb -t 8 ecoli.unicycler.consensus.fasta pb.reads.fastq > mapped.paf  
minimap2 -x map-pb -t 8 ecoli.unicycler.consensus.cut.fasta pb.reads.fastq > mapped.cut.paf  
```
**Filter to retain reads mapped to E. coli genome**  
*FASTA header is ecoli.genome and ecoli.genome.cut, respectively*  
```
grep ecoli.genome mapped.paf > genome.paf  
grep ecoli.genome.cut mapped.cut.paf > genome.cut.paf  
```
**Predict chimeras using Alvis**  
```
alvis/dist/Alvis.jar -type contigAlignment -inputfmt paf -outputfmt svg -chimeras -printChimeras -minChimeraCoveragePC 90 -minChimeraAlignmentPC 10 -in genome.paf -outdir /path/to/outdir/ -out out.prefix  
mv chimeras.txt genome.chimeras.txt  
alvis/dist/Alvis.jar -type contigAlignment -inputfmt paf -outputfmt svg -chimeras -printChimeras -minChimeraCoveragePC 90 -minChimeraAlignmentPC 10 -in genome.cut.paf -outdir /path/to/outdir/ -out out.prefix  
mv chimeras.txt genome.chimeras.cut.txt  
```
**Determine estimate for percentage chimeras**
*Identify reads assigned as chimeras in both original and cut E.coli genome*  
```
cat genome.chimeras.txt genome.cut.chimeras.txt | awk '{print $1}' - | sort -n | uniq -d > chimeras.candidates.txt  
```
**Re-map to produce BAM output, filter reads to include only chimeras, manual inspection**  
*ONT reads*  
```
minimap2 -ax map-ont -t 8 ecoli.unicycler.consensus.fasta minion.reads.fastq | samtools sort -o sorted.bam  
```
*PacBio reads*  
```
minimap2 -ax map-pb -t 8 ecoli.unicycler.consensus.fasta pb.reads.fastq | samtools sort -o sorted.bam  
```
*Filter reads*  
```
java -Xmx10g -jar picard-tools-2.5.0/picard.jar FilterSamReads I=sorted.bam O=chimeras.bam READ_LIST_FILE=chimeras_candidates.txt FILTER=includeReadList
samtools index chimeras.bam
```
**Determine estimate for percentage chimeras**  
*Total number of chimeric reads*
```
wc -l chimeras.candidates.txt  
```
*Total primary reads mapped to genome*  
```
samtools view sorted.bam -F 4 -F 256 -F 1024 -F 2048 ecoli.genome -c  
```

### E. coli plasmid analysis <a name="ecoli.plasmid"></a>  
**Produce depth counts for primary reads mapped to genome, pMAR2, p5217**  
*Map long reads*
```
minimap2 -ax map-ont -t 8 ecoli.unicycler.consensus.fasta ont.reads.fastq | samtools sort -o ont.sorted.bam 
minimap2 -ax map-pb -t 8 ecoli.unicycler.consensus.fasta pb.reads.fastq | samtools sort -o pb.sorted.bam  
```
*Filter long read BAM files to retain primary reads*  
```
samtools view sorted.bam -F 4 -F 256 -F 2048 -bho primary.bam  
```
*Calculcate sequencing depth*  
```
samtools depth -aa -m 100000000 primary.bam > primary.depth.txt 
```
*Map short read data*
```
bwa mem -k 23 -t 8 ecoli.unicycler.consensus.fasta illumina.reads.R1.fastq illumina.reads.R2.fastq | samtools sort -o sorted.illumina.bam 
```
*Mark Duplicates*  
```
picard.jar MarkDuplicates I=sorted.illumina.bam O=dupsmarked.illumina.bam M=dupsmarked.illumina.metrics VALIDATION_STRINGENCY=SILENT AS=true CREATE_INDEX=true
```
*Filter to retain primary reads*  
```
samtools view dupsmarked.illumina.bam -F 4 -F 256 -F 1024 -F 2048 -bho primary.illumina.bam  
```
*Calculcate sequencing depth*  
```
samtools depth -aa -m 100000000 primary.illumina.bam > primary.illumina.depth.txt  
```
**Count total depth across genome, pMAR2, p5217**  
```
grep ecoli.genome primary.depth.txt | awk '{total = total + $3}END{print "Total genome depth = "total}' -  
grep pMAR2 primary.depth.txt | awk '{total = total + $3}END{print "Total pMAR2 depth = "total}' -  
grep p5217 primary.depth.txt | awk '{total = total + $3}END{print "Total p5217 depth = "total}' -  
```
**Plasmid qPCR visualization** 
```
scripts/Ecoli_Long_Read_Analysis.Rmd
```

### Ecoli DNA modification analysis <a name="ecoli.dna.mod"></a>  

**PacBio RS II reads**
```
pbmm2 align ecoli.unicycler.consensus.fasta subreads.bam ecoli.PB.mapped_sorted.bam --sort -j 16 -J 8
samtools index ecoli.PB.mapped_sorted.bam
pbindex ecoli.PB.mapped_sorted.bam
ipdSummary ecoli.PB.mapped_sorted.bam --reference ecoli.unicycler.consensus.fasta --gff pb.basemods.gff --csv pb.basemods.csv --pvalue 0.001 --numWorkers 16 --identify m4C,m6A,m5C_TET

motifMaker find -f ecoli.unicycler.consensus.fasta -g pb.basemods.gff -o pb.motifs.csv

motifMaker reprocess -f ecoli.unicycler.consensus.fasta -g pb.basemods.gff -m pb.motifs.csv -o pb.motifs.gff
```
**PacBio Sequel II reads**
```
smrtlink_8.0.0.80529/smrtcmds/bin/dataset create --type SubreadSet --name ecoli.subreadset subreadset.xml pb.SQII.subreads.bam
smrtlink_8.0.0.80529/smrtcmds/bin/dataset create --type ReferenceSet --name ecoli.referenceset referenceset.xml ecoli.unicycler.consensus.fasta

smrtlink_8.0.0.80529/smrtcmds/bin/pbcromwell run pb_basemods -e subreadset.xml -e referenceset.xml -t kineticstools_compute_methyl_fraction=True -t kineticstools_identify_mods=m4C,m6A,m5C_TET -t run_find_motifs=True

smrtlink_8.0.0.80529/smrtcmds/bin/pbcromwell run pb_basemods -e subreadset.xml -e referenceset.xml -t kineticstools_compute_methyl_fraction=True -t kineticstools_identify_mods=m4C,m6A,m5C_TET -t run_find_motifs=True -t motif_min_score=125
```

**ONT reads**
```
tombo preprocess annotate_raw_with_fastqs --overwrite --fast5-basedir fast5_dir --fastq-filenames ont.reads.fastq --processes 16

tombo resquiggle fast5_dir ecoli.unicycler.consensus.fasta --processes 16 --num-most-common-errors 5

tombo detect_modifications de_novo --fast5-basedirs fast5_dir --statistics-file-basename ONT.denovo --processes 16

tombo text_output browser_files --fast5-basedirs fast5_dir --statistics-filename ONT.denovo.tombo.stats --file-types dampened_fraction --browser-file-basename ONT.denovo

tombo plot motif_with_stats --fast5-basedirs fast5_dir --statistics-filename ONT.denovo.tombo.stats --genome-fasta ecoli.unicycler.consensus.fasta --motif GATC --plot-standard-model --num-statistics 10000 --num-regions 1 --pdf-filename ONT.denovo.plot.GATC.pdf

tombo plot motif_with_stats --fast5-basedirs fast5_dir --statistics-filename ONT.denovo.tombo.stats --genome-fasta ecoli.unicycler.consensus.fasta --motif CCWGG --plot-standard-model --num-statistics 10000 --num-regions 1 --pdf-filename ONT.denovo.plot.CCWGG.pdf

tombo plot roc --statistics-filenames RAPID.5mC.tombo.stats RAPID.6mA.tombo.stats LIG.5mC.tombo.stats LIG.6mA.tombo.stats --motif-descriptions CCWGG:2:"CCWGG RAPID" GATC:2:"GATC RAPID" CCWGG:2:"CCWGG LIG"::GATC:2:"GATC LIG" --genome-fasta ecoli.unicycler.consensus.fasta
```

### D.ananassae genome assembly using Canu <a name="dana.canu"></a>
**MinION**  
```
canu -p output.prefix -d output.dir genomeSize=240m gridEngineThreadsOption="-pe thread THREADS" gridEngineMemoryOption="-l mem_free=MEMORY" gridOptions="-P jdhotopp-lab -q threaded.q" -nanopore-raw raw.nanopore.reads.fastq  
```
**PacBio**  
```
canu -p output.prefix -d output.dir genomeSize=240m gridEngineThreadsOption="-pe thread THREADS" gridEngineMemoryOption="-l mem_free=MEMORY" gridOptions="-P jdhotopp-lab -q threaded.q" -pacbio-raw raw.pacbio.reads.fastq  
```
**Hybrid**  
```
canu -p output.prefix -d output.dir genomeSize=240m gridEngineThreadsOption="-pe thread THREADS" gridEngineMemoryOption="-l mem_free=MEMORY" gridOptions="-P jdhotopp-lab -q threaded.q" -nanopore-raw raw.nanopore.reads.fastq -pacbio-raw raw.pacbio.reads.fastq  
```
**Genome polishing**  
```
pilon_iter.sh canu/assembly.contigs.fasta illuminaPE_R1.fastq.gz illuminaPE_R2.fastq.gz canu/assembly.trimmedReads.fasta  
```

### D. ananassae genome assessment <a name="dana.eval"></a>  
assemblathon_stats.pl contigs.fasta -csv -graph -genome_size 240000000
**D. ananassae QUAST-LG** 
```
use python-3.5
echo "/home/etvedte/scripts/quast-5.0.2/quast.py ../dana.hybrid.LIG+SQII.pilon.l_contigs.rn.fasta ../dana.hybrid.RAPID+SQII.pilon.l_contigs.fasta -r /local/projects-t3/RDBKO/nonIGS_dana/caf1/GCA_000005115.1_dana_caf1_genomic_scaffolds.fna -o quast_test -t 8 --large -k --est-ref-size 240000000 --fragmented" | qsub -P jdhotopp-lab -l mem_free=20G -q threaded.q -pe thread 8 -N quast.LG -cwd -V

echo "/home/etvedte/scripts/quast-5.0.2/quast.py /local/projects-t3/RDBKO/dana.postassembly/dana.minion.RAPID.raw_contigs.fasta -r /local/projects-t3/RDBKO/nonIGS_dana/caf1/GCA_000005115.1_dana_caf1_genomic.redux.fasta --features gene:/local/projects-t3/RDBKO/nonIGS_dana/caf1/GCA_000005115.1_dana_caf1_genomic.gff -o quast_minion_features -t 24 --large -e --est-ref-size 240000000 --fragmented" | qsub -P jdhotopp-lab -l mem_free=20G -q threaded.q -pe thread 24 -N quast.LG -cwd -V


echo "/home/etvedte/scripts/quast-5.0.2/quast.py /local/projects-t3/RDBKO/dana.postassembly/dana.minion.RAPID.raw_contigs.fasta /local/projects-t3/RDBKO/dana.postassembly/dana.minion.LIG.raw_contigs.fasta /local/projects-t3/RDBKO/dana.postassembly/dana.pb.sqII.raw_contigs.fasta /local/projects-t3/RDBKO/dana.postassembly/dana.hybrid.RAPID+SQII.raw_contigs.fasta /local/projects-t3/RDBKO/dana.postassembly/dana.hybrid.LIG+SQII.raw_contigs.fasta /local/projects-t3/RDBKO/dana.postassembly/dana.hybrid.RS+SQII.raw_contigs.fasta /local/projects-t3/RDBKO/nonIGS_dana/Miller2018/Dana.pass.minimap2.racon.x3.pilon.x3.fasta /local/projects-t3/RDBKO/nonIGS_dana/caf1/GCA_000005115.1_dana_caf1_genomic.redux.fasta -r /local/projects-t3/RDBKO/nonIGS_dana/caf1/GCA_000005115.1_dana_caf1_genomic.redux.fasta --features gene:/local/projects-t3/RDBKO/nonIGS_dana/caf1/GCA_000005115.1_dana_caf1_genomic.gff -o quast_all_assemblies -t 24 --large -m 0 --fragmented --split-scaffolds" | qsub -P jdhotopp-lab -l mem_free=20G -q threaded.q -pe thread 24 -N quast.LG -cwd -V


echo "kat comp -n -t 16 -o minion.RAPID.v.illumina /local/projects-t3/RDBKO/sequencing/RANDD_LIG_Dana_20190405_merged_pass.fastq.gz /local/projects-t3/RDBKO/sequencing/cHI_Dana_2_15_19_ILLUMINA_DATA/RANDD_20190322_K00134_IL100123454_MX29_L004_R1.fastq" | qsub -P jdhotopp-lab -l mem_free=10G -q threaded.q -pe thread 16 -N kat.comp.reads -cwd -V

```

non-chromosomal contigs
```
echo "/home/etvedte/scripts/quast-5.0.2/quast.py /local/projects-t3/RDBKO/dana.postassembly/dana.minion.RAPID.raw_contigs.nonchr.fasta /local/projects-t3/RDBKO/dana.postassembly/dana.minion.LIG.raw_contigs.nonchr.fasta /local/projects-t3/RDBKO/dana.postassembly/dana.pb.sqII.raw_contigs.nonchr.fasta /local/projects-t3/LGT/Dananassae_2020/dana.quickmerge/flye+canu.FREEZE.custom.params/pilon.long.bases/dana.assembly.FREEZE.plusMITO.6.1.20.nonchr.fasta /local/projects-t3/RDBKO/nonIGS_dana/Miller2018/Dana.pass.minimap2.racon.x3.pilon.x3.rh.nonchr.fasta /local/projects-t3/RDBKO/nonIGS_dana/caf1/GCA_000005115.1_dana_caf1_genomic_scaffolds.nonchr.fna -o quast_all_nonchr -t 24 --large -m 0 --split-scaffolds" | qsub -P jdhotopp-lab -l mem_free=20G -q threaded.q -pe thread 24 -N quast.LG -cwd -V
```

**KAT**
```
use kat-2.4.0
echo "kat comp -t 16 -o minion.RAPID.assembly /local/projects-t3/RDBKO/sequencing/cHI_Dana_2_15_19_ILLUMINA_DATA/RANDD_20190322_K00134_IL100123454_MX29_L004_R1.fastq /local/projects-t3/RDBKO/dana.postassembly/dana.minion.RAPID.pilon.l_contigs.fasta" | qsub -P jdhotopp-lab -l mem_free=10G -q threaded.q -pe thread 16 -N kat.comp -cwd -V
kat plot spectra-cn -x 150 -p pdf -o Miller2018_cHi_illumina-main.mx.spectra-cn.pdf Miller2018_cHi_illumina-main.mx

**D. ananassae BUSCO**  
```
python run_BUSCO.py -f -c 8 -t /local/scratch/etvedte/tmp -i assembly.fasta -o busco_output_dir -l metazoa_odb9 -m geno  
```

### Alignment of contigs to D. ananassae polytene map <a name="dana.map"></a>  
*Initial BLAST against polished contigs to determine contig orientation*  
```
makeblastdb -in polished.contigs.fasta -out polished.contigs.fasta -dbtype nucl -parse_seqids  
blastn -query mapping_loci.fasta -db polished.contigs.fasta -max_target_seqs 10 -max_hsps 10 -outfmt 6 > initial.blast.out  
```
*After manual inspection, reverse complementation if necessary* 
```
samtools faidx polished.contigs.fasta fwd.contig.name >> polished.contigs.correct.fasta  
samtools faidx -i polished.contigs.fasta rev.contig.name >> polished.contigs.correct.fasta  
```
*Second round of BLAST. Can limit target sequences if alignment lengths of BLAST matches are similar to queries. Also used custom output columns*  
```
makeblastdb -in polished.contigs.correct.fasta -out polished.contigs.correct.fasta -dbtype nucl -parse_seqids  
blastn -query mapping_loci.fasta -db polished.contigs.correct.fasta -max_target_seqs 1 -max_hsps 1 -outfmt "6 qseqid sseqid pident length sstart send evalue slen" > final.blast.out  
```

**filter out chromosomes, redo QUAST stats
/usr/local/packages/bbtools/filterbyname.sh in=dana.minion.RAPID.raw_contigs.fasta out=dana.minion.RAPID.raw_contigs.nonchr.fasta names=dana.minion.RAPID.raw_contigs.chr.list include=f


**Rename FASTA**  
```
for f in *_contigs.fasta; do awk '/^>/{print ">ecoli_contig" ++i; next}{print}' < $f > ${f%_c*}_contigs_rn.fasta; done
```


#Dana basemod
tombo text_output signif_sequence_context --fast5-basedirs fast5_dir --statistics-filename ONT.denovo.tombo.stats --num-regions 1000 --num-bases 50

echo -e "/usr/local/packages/meme-4.12.0/bin/meme -oc RANDD_RAPID_Ecoli.tombo.stats.dam.meme -dna -mod zoops -nmotifs 50 tombo_results.significant_regions.fasta" | qsub -P jdhotopp-lab -l mem_free=5G -N meme -cwd


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
