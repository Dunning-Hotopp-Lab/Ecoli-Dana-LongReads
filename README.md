# Comparison of long read sequencing technologies in interrogating bacteria and fly genomes

Eric S. Tvedte

2020-09-14

## Table of Contents
1. [Prepare sequencing files for assembly](#ecoli.prep)
2. [Read length distributions](#ecoli.read)
3. [E. coli genome assembly](#ecoli.asm)
4. [Evaluation of E. coli genome assemblies and chimeras](#ecoli.eval)
5. [Ecoli DNA modification analysis](#ecoli.dna.mod)  
6. [E. coli plasmid quantification](#ecoli.plasmid)  
7. [D. ananassae genome assembly](#dana.canu)
8. [Dana.UMIGS genome assembly](#dana.umigs)
9. [Anchoring D. ananassae assembly contigs](#dana.anchor)
10. [Evaluation of D. ananassae genome assemblies](#dana.eval)  
11. [Detection of DNA modification in D. ananassae](#dana.dna.mod)
12. [R scripts for data visualization](#data.vis)


### Prepare sequencing files for assembly <a name="ecoli.prep"></a>
**Basecalling R9 ONT reads**
```
guppy_basecaller --input_path fast5_dir --save_path output_dir --config dna_r9.4.1_450bps_hac.cfg --fast5_out --post_out --qscore_filtering --min_qscore 7 --records_per_fastq 10000000 -x "cuda:5 cuda:6 cuda:7" > guppy.4.2.2.log
```
**Basecalling R10 ONT reads**
```
guppy_basecaller --input_path fast5_dir --save_path output_dir --config dna_r10.3_450bps_hac.cfg --fast5_out --post_out --qscore_filtering --min_qscore 7 --records_per_fastq 10000000 -x "cuda:5 cuda:6 cuda:7" > guppy.4.2.2.log
```
**Remove DNA Control Sequence from ONT LIG reads**
```
cat ONT.LIG.fastq | NanoLyse --reference DCS.fasta | gzip > ONT.LIG.filterDCS.fastq.gz 
```

### Read length distributions <a name="ecoli.read"></a>
**Determine read length distributions** 
```
bbtools/readlength.sh in=reads.fastq.gz out=reads.hist.out bin=1000 max=1000000 qin=33
grep -v '#' reads.hist.out > reads.hist.final.out
```
**Subsample read sets to 100X**
```
seqkit sample -j 8 -p 0.0339 -o ONT.LIG.100X.fastq ONT.LIG.filterDCS.fastq.gz
seqkit sample -j 8 -p 0.3677 -o Ecoli.PB.RSII.100X.fastq Ecoli.PB.RSII.fastq.gz
seqkit sample -j 8 -p 0.0078 -o Ecoli.PB.CLR.100X.fastq Ecoli.PB.CLR.fastq.gz
seqkit sample -j 8 -p 0.0216 -o Ecoli.PB.HiFi.100X.fastq Ecoli.PB.HiFi.fastq.gz
```

**Generate histograms in R**  
Perform analysis for all reads and subsampled reads  
Input: binned read count frequencies from bbtools, set as in.path
Output: read length histograms
```r
library(ggplot2)
library(gtable)
library(grid)
library(scales)

all.read.data = data.frame()
replicate.data = data.frame()

#PacBio Sequel II
in.path = "pb.reads.hist.final.out"
data <- read.table(in.path, header = F, sep = '\t')
colnames(data) <- c("length", "reads", "pct_reads", "cum_reads", "cum_pct_reads", "bases", "pct_bases", "cum_bases", "cum_pct_bases")
data[,"rlib"] <- "PB_SEQUELII"
data$length[1] <- 1
data[,"length.mean"] <- mean(data$length)
data[,"length.max"] <- max(data$length)

all.read.data <- rbind(all.read.data,data)

#ONT RAPID
in.path = "ont.reads.hist.final.out"
data <- read.table(in.path, header = F, sep = '\t')
colnames(data) <- c("length", "reads", "pct_reads", "cum_reads", "cum_pct_reads", "bases", "pct_bases", "cum_bases", "cum_pct_bases")
data[,"rlib"] <- "ONT_RAPID"
data$length[1] <- 1
data[,"length.mean"] <- mean(data$length)
data[,"length.max"] <- max(data$length)

all.read.data <- rbind(all.read.data,data)

#convert to numeric
all.read.data$cum_pct_bases <- as.character(all.read.data$cum_pct_bases)
all.read.data$cum_pct_bases <- gsub("%", "", all.read.data$cum_pct_bases)
all.read.data$cum_pct_bases <- as.numeric(all.read.data$cum_pct_bases)
gsub("%", "", as.character(factor(all.read.data$cum_pct_bases)))

#convert to numeric
all.read.data$pct_bases <- as.character(all.read.data$pct_bases)
all.read.data$pct_bases <- gsub("%", "", all.read.data$pct_bases)
all.read.data$pct_bases <- as.numeric(all.read.data$pct_bases)
gsub("%", "", as.character(factor(all.read.data$pct_bases)))

#convert to numeric
all.read.data$cum_pct_reads <- as.character(all.read.data$cum_pct_reads)
all.read.data$cum_pct_reads <- gsub("%", "", all.read.data$cum_pct_reads)
all.read.data$cum_pct_reads <- as.numeric(all.read.data$cum_pct_reads)
gsub("%", "", as.character(factor(all.read.data$cum_pct_reads)))

#convert to numeric
all.read.data$pct_reads <- as.character(all.read.data$pct_reads)
all.read.data$pct_reads <- gsub("%", "", all.read.data$pct_reads)
all.read.data$pct_reads <- as.numeric(all.read.data$pct_reads)
gsub("%", "", as.character(factor(all.read.data$pct_reads)))

#add small numbers so log transformation can be used
all.read.data$cum_pct_bases <- (100 - all.read.data$cum_pct_bases) + 0.001
all.read.data$cum_pct_reads <- 100 - all.read.data$cum_pct_reads

top_theme = theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_text(size = 8), axis.text.y = element_text(size = 8), legend.position = "none")

p1 <- ggplot(all.read.data, aes(x=all.read.data$length, y=reads, colour = rlib)) + 
    geom_line() +
    xlab("read length (bp)") +
    ylab("read count") +
    scale_x_continuous(trans = "log10", limits = c(1000,1000000)) +
    #scale_x_continuous(limits = c(0,250000)) +
    scale_y_continuous(label = comma, trans = "log10", limits = c(1,1000000)) +
    scale_color_manual(values=c('#1b9e77','#1b9e7770', '#d95f02','#d95f0270', '#e7298a','#7570b3')) +
    geom_vline(aes(xintercept=length.max, color=rlib), linetype="dashed") +
    theme_bw() +  
    top_theme

p2 <- ggplot(all.read.data, aes(x=all.read.data$length, y=pct_reads, colour = rlib)) + 
    geom_line() +
    xlab("read length (bp)") +
    ylab("read count (% of total)") +
    scale_x_continuous(trans = "log10", limits = c(1000,1000000)) +
    #scale_x_continuous(limits = c(0,250000)) +
    scale_y_continuous(limits = c(0,20)) +
    scale_color_manual(values=c('#1b9e77','#1b9e7770', '#d95f02','#d95f0270', '#e7298a','#7570b3')) +
    geom_vline(aes(xintercept=length.max, color=rlib), linetype="dashed") +
    theme_bw() +  
    top_theme

p3 <- ggplot(all.read.data, aes(x=length, y=bases, colour = rlib)) + 
    geom_line() +
    xlab("read length (bp)") +
    ylab("bases sequenced") +
    scale_x_continuous(label = comma, trans = "log10", limits = c(1000,1000000)) +
    #scale_x_continuous(limits = c(0,250000)) +
    theme(axis.ticks = element_blank(), axis.text.y = element_blank()) + 
    geom_vline(aes(xintercept=length.max, color=rlib), linetype="dashed") +
    #scale_color_brewer(palette = "Dark2") +
    scale_color_manual(values=c('#1b9e77','#1b9e7770', '#d95f02','#d95f0270', '#e7298a','#7570b3')) +
    theme_bw() + 
    theme(axis.title.y = element_text(size = 8), axis.text.y = element_text(size = 8), legend.position = "bottom")

p4 <- ggplot(all.read.data, aes(x=length, y=pct_bases, colour = rlib)) + 
    geom_line() +
    xlab("read length (bp)") +
    ylab("bases sequenced (% of total)") +
    scale_x_continuous(label = comma, trans = "log10", limits = c(1000,1000000)) +
    #scale_x_continuous(limits = c(0,250000)) +
    scale_y_continuous(limits = c(0,20)) +
    theme(axis.ticks = element_blank(), axis.text.y = element_blank()) + 
    geom_vline(aes(xintercept=length.max, color=rlib), linetype="dashed") +
    #scale_color_brewer(palette = "Dark2") +
    scale_color_manual(values=c('#1b9e77','#1b9e7770', '#d95f02','#d95f0270', '#e7298a','#7570b3')) +
    theme_bw() + 
    theme(axis.title.y = element_text(size = 8), axis.text.y = element_text(size = 8), legend.position = "bottom")

#optional pdf output
pdf("X:/RDBKO/raw.read.composition/Ecoli.readcomp.full.pdf", width = 6.5, height = 5)
   g1 <- ggplotGrob(p1)
   g2 <- ggplotGrob(p2)
   g3 <- ggplotGrob(p3)
   g4 <- ggplotGrob(p4)
  
   g5 <- rbind(g1, g3, size = "first")
   g6 <- rbind(g2, g4, size = "first")
  
   g5$widths <- unit.pmax(g1$widths, g4$widths)
   g6$widths <- unit.pmax(g2$widths, g5$widths)
   
   g <- cbind(g5, g6, size = "first")
   grid.newpage()
   grid.draw(g)

#optional pdf close 
dev.off()

```

### E. coli genome assembly <a name="ecoli.asm"></a>

**Canu ONT assembly**  
```
canu -p output.prefix -d output.dir genomeSize=4.6m corOutCoverage=1000 gridEngineThreadsOption="-pe thread THREADS" gridEngineMemoryOption="-l mem_free=MEMORY" gridOptions="-P jdhotopp-lab -q threaded.q" -nanopore ont.reads.fastq 
```
**Canu PacBio RSII/CLR assembly** 
```
canu -p output.prefix -d output.dir genomeSize=4.6m gridEngineThreadsOption="-pe thread THREADS" gridEngineMemoryOption="-l mem_free=MEMORY" gridOptions="-P jdhotopp-lab -q threaded.q" -pacbio pacbio.reads.fastq
```
**HiCanu PacBio HiFi assembly** 
```
canu -p output.prefix -d output.dir genomeSize=4.6m gridEngineThreadsOption="-pe thread THREADS" gridEngineMemoryOption="-l mem_free=MEMORY" gridOptions="-P jdhotopp-lab -q threaded.q" -pacbio-hifi pacbio.hifi.reads.fastq
```

**Polish, circularize, and rotate**  
```
bwa index genome.fasta
bwa mem -t 16 genome.fasta fwd.R1.fastq rev.R2.fastq | samtools view -@ 16 -bS - |samtools sort -@ 16 -o ShortRead.bam -
pilon-1.22.jar --genome genome.fasta --frags ShortRead.bam --changes --threads 16 --minmq 10 --fix bases --output pilon       
circlator minimus2 pilon.fasta circularise.fasta  
circlator fixstart --genes_fa ecoli.dnaA.fasta circularise.fasta rotated.fasta  
```
**Flye ONT assembly**
```
flye -t 24 --plasmids -o output.dir --nano-raw ont.reads.fastq
```
**Flye PacBio RSII/CLR assembly**
```
flye -t 24 --plasmids -o output.dir --pacbio-raw pacbio.reads.fastq
```
**Flye PacBio HiFi assembly**
```
flye -t 24 --plasmids -o output.dir --pacbio-hifi pacbio.hifi.reads.fastq
```
**Unicycler hybrid assembly**  
```
unicycler --mode normal --start_genes ecoli.dnaA.protein.fasta --long long.reads.fastq --short1 fwd.R1.fastq --short2 rev.R2.fastq -o output.dir -t 32  
```

### E. coli genome assembly quality <a name="ecoli.eval"></a>   
**BUSCO: conserved genes** 
```
busco --config config.ini -i contigs.fasta -c 8 -o output_dir -l bacteria -m geno  
```
**QUAST: consensus identity, missassemblies**
```
quast.py -m 0 -t 8 contigs.1.fasta contigs.2.fasta contigs.3.fasta -o output_dir -r Ecoli.UMIGS.fasta
```
**Alvis: chimeric reads**  
*Map ONT reads*  
```
minimap2 -x map-ont -t 8 Ecoli.UMIGS.fasta ont.reads.fastq > mapped.paf  
minimap2 -x map-ont -t 8 Ecoli.UMIGS.cut.fasta ont.reads.fastq > mapped.cut.paf  
```
*Map PacBio reads*  
```
minimap2 -x map-pb -t 8 Ecoli.UMIGS.fasta pb.reads.fastq > mapped.paf  
minimap2 -x map-pb -t 8 Ecoli.UMIGS.cut.fasta pb.reads.fastq > mapped.cut.paf  
```
*Filter to retain reads mapped to E. coli genome*  
```
grep ecoli.genome mapped.paf > genome.paf  
grep ecoli.genome.cut mapped.cut.paf > genome.cut.paf  
```
*Predict chimeras using Alvis*  
```
alvis/dist/Alvis.jar -type contigAlignment -inputfmt paf -outputfmt svg -chimeras -printChimeras -minChimeraCoveragePC 90 -minChimeraAlignmentPC 10 -in genome.paf -outdir /path/to/outdir/ -out out.prefix  
mv chimeras.txt genome.chimeras.txt  
alvis/dist/Alvis.jar -type contigAlignment -inputfmt paf -outputfmt svg -chimeras -printChimeras -minChimeraCoveragePC 90 -minChimeraAlignmentPC 10 -in genome.cut.paf -outdir /path/to/outdir/ -out out.prefix  
mv chimeras.txt genome.chimeras.cut.txt  
```
*Identify reads assigned as chimeras in both original and cut-and-paste E.coli genome*  
```
cat genome.chimeras.txt genome.cut.chimeras.txt | awk '{print $1}' - | sort -n | uniq -d > chimeras.candidates.txt  
```
*Re-map to produce BAM output, filter reads to include only chimeras, manual inspection*  
*ONT reads*  
```
minimap2 -ax map-ont -t 8 Ecoli.UMIGS.fasta ont.reads.fastq | samtools sort -o sorted.bam
minimap2 -ax map-pb -t 8 Ecoli.UMIGS.fasta pb.reads.fastq | samtools sort -o sorted.bam 
picard.jar FilterSamReads I=sorted.bam O=chimeras.bam READ_LIST_FILE=chimeras_candidates.txt FILTER=includeReadList
samtools index chimeras.bam
```
*Determine estimate for percentage chimeras*  
```
wc -l chimeras.candidates.txt #total chimeras mapped to E.coli genome
samtools view sorted.bam -F 4 -F 256 -F 2048 ecoli.genome -c  #total primary reads mapped to E.coli genome
```
### E. coli plasmid quantification <a name="ecoli.plasmid"></a>  
**Minimap2: map long reads to Ecoli.UMIGS**  
```
minimap2 -ax map-ont -t 8 Ecoli.UMIGS.fasta ont.reads.fastq | samtools sort -o ont.sorted.bam 
minimap2 -ax map-pb -t 8 Ecoli.UMIGS.fasta pb.reads.fastq | samtools sort -o pb.sorted.bam  
```
*Filter long read BAM files to retain primary reads, calculate sequencing depth*  
```
samtools view sorted.bam -F 4 -F 256 -F 2048 -bho primary.bam 
samtools depth -aa -m 100000000 primary.bam > primary.depth.txt
```
**BWA-MEM: map short reads to Ecoli.UMIGS**
```
bwa mem -k 23 -t 8 Ecoli.UMIGS.fasta fwd.R1.fastq rev.R2.fastq | samtools sort -o sorted.illumina.bam 
```
*Mark duplicates, filter to retain primary reads, calculcate sequencing depth*  
```
picard.jar MarkDuplicates I=sorted.illumina.bam O=dupsmarked.illumina.bam M=dupsmarked.illumina.metrics VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=TRUE AS=TRUE CREATE_INDEX=TRUE  
samtools view dupsmarked.illumina.bam -F 4 -F 256 -F 2048 -bho primary.illumina.bam  
samtools depth -aa -m 100000000 primary.illumina.bam > primary.illumina.depth.txt
```
**Count total depth across genome, pMAR2, p5217**  
```
grep ecoli.genome primary.depth.txt | awk '{total = total + $3}END{print "Total genome depth = "total}' -  
grep pMAR2 primary.depth.txt | awk '{total = total + $3}END{print "Total pMAR2 depth = "total}' -  
grep p5217 primary.depth.txt | awk '{total = total + $3}END{print "Total p5217 depth = "total}' -  
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

**Plot PacBio modification QV versus ONT dampened fraction**
*PacBio*
```
grep 'm6A' cromwell_out/cromwell-executions/pb_basemods/cb4f28f5-fff9-426a-8010-a03646fb69de/call-reprocess_motifs/execution/motifs.gff | awk '{print $1"\t"$4"\t"$6}' > pb.sqII.basemods.m6A.tsv
```
*ONT*
```
tombo detect_modifications alternative_model --alternate-bases 6mA --fast5-basedirs fast5_dir --statistics-file-basename ONT.6mA.alt --processes 16

tombo text_output browser_files --fast5-basedirs fast5_dir --statistics-filename ONT.6mA.alt.tombo.stats --file-types dampened_fraction --browser-file-basename ONT.denovo

wig2bed < ONT.6mA.alt.dampened_fraction_modified_reads.plus.wig > ONT.6mA.alt.dampened_fraction_modified_reads.plus.bed  
awk '{print $1"\t"$3"\t"$5}' ONT.6mA.alt.dampened_fraction_modified_reads.plus.bed > ONT.6mA.alt.dampened_fraction_modified_reads.plus.final.tsv  
wig2bed < ONT.6mA.alt.dampened_fraction_modified_reads.minus.wig > ONT.6mA.alt.dampened_fraction_modified_reads.minus.bed  
awk '{print $1"\t"$3"\t"$5}' ONT.6mA.alt.dampened_fraction_modified_reads.minus.bed > ONT.6mA.alt.dampened_fraction_modified_reads.minus.final.tsv  
```
*Visualization*
```
scripts/Ecoli_Long_Read_Analysis.Rmd
```



### D. ananassae genome assembly <a name="dana.canu"></a>
**Canu**
*ONT*
```
canu -p output.prefix -d output.dir genomeSize=240m gridEngineThreadsOption="-pe thread THREADS" gridEngineMemoryOption="-l mem_free=MEMORY" gridOptions="-P jdhotopp-lab -q threaded.q" -nanopore ont.reads.fastq  
```
*PacBio CLR*  
```
canu -p output.prefix -d output.dir genomeSize=240m gridEngineThreadsOption="-pe thread THREADS" gridEngineMemoryOption="-l mem_free=MEMORY" gridOptions="-P jdhotopp-lab -q threaded.q" -pacbio pacbio.reads.fastq  
```
*Hybrid ONT and CLR*
```
canu -p output.prefix -d output.dir genomeSize=240m gridEngineThreadsOption="-pe thread THREADS" gridEngineMemoryOption="-l mem_free=MEMORY" gridOptions="-P jdhotopp-lab -q threaded.q" -pacbio pacbio.reads.fastq -nanopore ont.reads.fastq
```
*PacBio HiFi*
```
canu -p output.prefix -d output.dir genomeSize=240m gridEngineThreadsOption="-pe thread THREADS" gridEngineMemoryOption="-l mem_free=MEMORY" gridOptions="-P jdhotopp-lab -q threaded.q" -pacbio-hifi pacbio.hifi.reads.fastq
```
**Flye**
*ONT*
```
flye --asm-coverage 40 --genome-size 240m --nano-raw ont.reads.fastq -t 24 -o output_dir  
```
*PacBio CLR*  
```
flye --asm-coverage 40 --genome-size 240m --pacbio-raw pb.reads.fastq -t 24 -o output_dir  
```
*PacBio HiFi*
```
flye --asm-coverage 40 --genome-size 240m --pacbio-hifi pb.hifi.reads.fastq -t 24 -o output_dir
```

**Genome polishing (performed on PacBio Sequel II + ONT LIG hybrid assembly)**  
```
Modify pilon_iter.sh for polishing with short reads/long reads/both
scripts/pilon_iter.sh canu/assembly.contigs.fasta illuminaPE_R1.fastq.gz illuminaPE_R2.fastq.gz canu/assembly.trimmedReads.fasta  
```
### Dana.UMIGS genome assembly <a name="dana.umigs"></a>
**Canu hybrid**
```
canu -p output.prefix -d output.dir genomeSize=240m corOutCoverage=80 gridEngineThreadsOption="-pe thread THREADS" gridEngineMemoryOption="-l mem_free=MEMORY" gridOptions="-P jdhotopp-lab -q threaded.q" -pacbio-raw pacbio.sequelII.reads.fastq.gz -nanopore-raw ont.LIG.reads.fastq.gz
```
**Polish with Arrow (2X)**  
```
pbmm2 align canu.assembly.contigs.fasta pacbio.sequelII.subreads.bam canu.assembly.mapped.pb.sqII_sorted.bam --sort -j 16 -J 8

```
**polish with arrow using Sequel II data**
```
echo "/usr/local/packages/smrttools/install/current/bundles/smrttools/smrtcmds/bin/pbmm2 align Circularized_assembly_1_dana.illumina.mito.NOVOPlasty.fasta /local/projects-t3/RDBKO/sequencing/cHI_Dana_2_15_19_PACBIO_DATA/RANDD_20190301_S64018_PL100122512-1_C01.subreads.bam Circularized_assembly_1_dana.illumina.mito.NOVOPlasty.sqII.arrow1_sorted.bam --sort -j 16 -J 8" | qsub -P jdhotopp-lab -l mem_free=50G -N pbmm2.align -q threaded.q -pe thread 16 -cwd -V
```

**Flye Sequel II**
```
flye -g 240m -t 24 -o flye_assembly_dir --asm-coverage 60 --pacbio-raw pacbio.sequelII.reads.fastq.gz
```

**Map PacBio HiFi data**  



**purge haplotigs from assembly**
```
echo "minimap2 -xmap-pb /local/projects-t3/LGT/Dananassae_2020/dana.postassembly/arrow/sqII.rd2/dana.hybrid.80X.arrow.rd2.contigs.fasta /local/projects-t3/RDBKO/sequencing/Dana.Hawaii.pbSequelII.raw.fastq.gz | gzip -c - > dana.hybrid.80X.arrow.rd2.mappedsqII.paf.gz" | qsub -P jdhotopp-lab -l mem_free=10G -N minimap2 -cwd

/local/projects-t3/RDBKO/scripts/purge_dups/bin/split_fa /local/projects-t3/LGT/Dananassae_2020/dana.postassembly/arrow/sqII.rd2/dana.hybrid.80X.arrow.rd2.contigs.fasta > /local/projects-t3/LGT/Dananassae_2020/dana.postassembly/arrow/sqII.rd2/dana.hybrid.80X.arrow.rd2.contigs.split

/home/etvedte/scripts/purge_dups/bin/pbcstat dana.hybrid.80X.contigs.arrow.polished.mappedhifi.paf.gz
/home/etvedte/scripts/purge_dups/bin/calcuts PB.stat > cutoffs 2> calcuts.log
/home/etvedte/scripts/purge_dups/scripts/hist_plot.py PB.stat hist.out.pdf
```

**Convert bases to upper case**
*By default arrow outputs regions with no consensus as lower case, i.e. 'acgt'. In order to properly annotate repetitive regions as lower case, all bases must be converted to upper case. Note that this could have been accomplished in arrow using the parameter --noEvidenceConsensusCall reference* 
```
awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' dana.hybrid.80X.arrow.rd2.contigs.FREEZE.fasta > dana.hybrid.80X.arrow.rd2.contigs.FREEZE.fmt.fasta
```

### D. ananassae post-assembly processing <a name="dana.post"></a>
**purge_dups**
```
for f in /local/projects-t3/RDBKO/dana.flye/*/*contigs.fasta; do echo "minimap2 -xmap-pb $f /local/projects-t3/RDBKO/sequencing/cHI_Dana_2_15_19_PACBIO_DATA_HiFi/cHI_Dana_2_15_19/PACBIO_DATA/RANDD_20191011_S64018_PL100122512-3_A01.ccs.fastq.gz > /local/projects-t3/RDBKO/dana.postassembly/purge_dups/$(basename ${f%_s*}).paf" | qsub -P jdhotopp-lab -l mem_free=5G -N minimap2 -cwd -V; done
/home/etvedte/scripts/purge_dups/bin/pbcstat PB.HiFi.canu.contigs.fasta.paf
/home/etvedte/scripts/purge_dups/scripts/hist_plot.py PB.stat PB.stat.hist.jpeg
/home/etvedte/scripts/purge_dups/bin/calcuts -d 1 -l5 -m195 -u300 PB.stat > cutoffs 2>calcuts.log

for f in /local/projects-t3/RDBKO/dana.flye/*/*contigs.fasta; do /home/etvedte/scripts/purge_dups/bin/split_fa $f > /local/projects-t3/RDBKO/dana.postassembly/purge_dups/$(basename ${f%_s*}).split; done
for f in *split; do echo "minimap2 -xasm5 -DP $f $f > $f.self.paf" | qsub -P jdhotopp-lab -l mem_free=5G -N minimap2 -cwd -V; done

/home/etvedte/scripts/purge_dups/bin/purge_dups -2 -T cutoffs -c PB.base.cov PB.CLR.canu.raw.contigs.fasta.split.self.paf > dups.bed 2> purge_dups.log
/home/etvedte/scripts/purge_dups/bin/get_seqs -e dups.bed /local/projects-t3/RDBKO/dana.postassembly/PB.CLR.canu.raw.contigs.fasta
```
**purge_haplotigs** 
```
for f in /local/projects-t3/RDBKO/dana.flye/*/*contigs.fasta; do echo minimap2 -t 4 -ax map-pb $f /local/projects-t3/RDBKO/sequencing/cHI_Dana_2_15_19_PACBIO_DATA_HiFi/cHI_Dana_2_15_19/PACBIO_DATA/RANDD_20191011_S64018_PL100122512-3_A01.ccs.fastq.gz --secondary=no | samtools sort -m 5G -o /local/projects-t3/RDBKO/dana.postassembly/purge_haplotigs/$(basename ${f%_s*}).HiFi.aligned.bam -T /local/scratch/etvedte; done
echo "/home/etvedte/scripts/purge_haplotigs/bin/purge_haplotigs hist -b PB.CLR.canu.raw.contigs.fasta.HiFi.aligned.bam -g /local/projects-t3/RDBKO/dana.postassembly/PB.CLR.canu.raw.contigs.fasta -t 4 -d 400" | qsub -P jdhotopp-lab -l mem_free=5G -q threaded.q -pe thread 4 -N purge.hist -cwd -V
/home/etvedte/scripts/purge_haplotigs/bin/purge_haplotigs cov -l 5 -m 195 -h 300 -i PB.CLR.canu.raw.contigs.fasta.HiFi.aligned.bam.gencov -j 80 -s 80
echo "/home/etvedte/scripts/purge_haplotigs/bin/purge_haplotigs purge -g /local/projects-t3/RDBKO/dana.postassembly/PB.CLR.canu.raw.contigs.fasta -b PB.CLR.canu.raw.contigs.fasta.HiFi.aligned.bam -c coverage_stats.csv -d -t 8" | qsub -P jdhotopp-lab -l mem_free=10G -q threaded.q -pe thread 8 -N ph.purge -cwd -V
echo "/home/etvedte/scripts/purge_haplotigs/bin/purge_haplotigs clip -p PB.CLR.flye.curated.fasta -h curated.haplotigs.fasta -t 4" | qsub -P jdhotopp-lab -l mem_free=5G -q threaded.q -pe thread 4 -N ph.clip -cwd -V
```
**purge_haplotigs final ** 
```
minimap2 -t 4 -ax map-pb asm.fasta pb.hifi.fastq.gz --secondary=no | samtools sort -m 5G -o aligned.bam
purge_haplotigs hist -b aligned.bam -g asm.fasta -t 4 -d 400
purge_haplotigs cov -l 5 -m 195 -h 300 -i aligned.bam.gencov -j 80 -s 80
purge_haplotigs purge -g asm.fasta -b aligned.bam -c coverage_stats.csv -d -t 8
purge_haplotigs clip -p curated.fasta -h curated.haplotigs.fasta -t 4
```
## Anchoring D. ananassae assembly contigs <a name="dana.anchor"></a>
**Major chromosome arm contigs (X, 2, 3)**
```
blastn -query chrom.map.loci.fasta -db asm.fasta -max_target_seqs 5 -max_hsps 5 -outfmt "6 qseqid sseqid pident length sstart send evalue slen" > initial.blast.out #inspect output to generate list of chromosome arm sequences
seqkit grep -f chr.arm.contigs.list asm.fasta > chr.arm.contigs.fasta
blastn -query chrom.map.loci.fasta -db chr.arm.contigs.fasta -max_target_seqs 1 -max_hsps 1 -outfmt "6 qseqid sseqid pident length sstart send evalue slen" > final.blast.out
```

**Chromosome 4 contigs and LGT contigs**
*Chromosome 4*
```

```

*LGT*
```
nucmer -l 1000 --prefix asm.LGT wAna.genome.fasta asm.fasta
delta-filter asm.LGT.delta -q > asm.LGT.filter
show-coords -r asm.LGT.filter -T > asm.LGT.filter.coords
tail -n +5 asm.LGT.delta.coords | awk '{print $13}' | sort -n | uniq > LGT.contigs.list #retrieves set of LGT contigs



tail -n +5 finalpass.coords | awk '{print $9"\t"$3"\t"$4}' > finalpass.bed
Rscript /home/etvedte/scripts/fixbed.R finalpass.bed fixed.bed
bedtools coverage -a fixed.bed -b fixed.bed -hist | grep all | awk '{total = total + $3}END{print "LGT segment length = "total}' #estimates alignment block length of wAna to LGT


```

**Chromosome Y contigs**
```

```

### D. ananassae genome assessment <a name="dana.eval"></a>  

**D. ananassae QUAST-LG** 
```
quast.py ONT.RAPID.raw_contigs.fasta ONT.LIG.raw_contigs.fasta PB.SequelII.raw_contigs.fasta Hybrid.RAPID+SequelII.raw_contigs.fasta Hybrid.LIG+SequelII.raw_contigs.fasta Hybrid.RSII+SequelII.raw_contigs.fasta Hybrid.LIG+SequelII.polished.short_contigs.fasta Hybrid.LIG+SequelII.polished.long_contigs.fasta Hybrid.LIG+SequelII.polished_both_contigs.fasta Hybrid.LIG+SequelII.polished_hifi_contigs.fasta Dana.Miller2018.fasta Dana.caf1.fasta Dana.UMIGS_contigs.fasta -r Dana.UMIGS_contigs.fasta -o quast_output_dir -t 24 --large -m 0 --conserved-genes-finding --split-scaffolds
```

non-chromosomal contigs
```
echo "/home/etvedte/scripts/quast-5.0.2/quast.py /local/projects-t3/RDBKO/dana.postassembly/dana.minion.RAPID.raw_contigs.nonchr.fasta /local/projects-t3/RDBKO/dana.postassembly/dana.minion.LIG.raw_contigs.nonchr.fasta /local/projects-t3/RDBKO/dana.postassembly/dana.pb.sqII.raw_contigs.nonchr.fasta /local/projects-t3/LGT/Dananassae_2020/dana.quickmerge/flye+canu.FREEZE.custom.params/pilon.long.bases/dana.assembly.FREEZE.plusMITO.6.1.20.nonchr.fasta /local/projects-t3/RDBKO/nonIGS_dana/Miller2018/Dana.pass.minimap2.racon.x3.pilon.x3.rh.nonchr.fasta /local/projects-t3/RDBKO/nonIGS_dana/caf1/GCA_000005115.1_dana_caf1_genomic_scaffolds.nonchr.fna -o quast_all_nonchr -t 24 --large -m 0 --split-scaffolds" | qsub -P jdhotopp-lab -l mem_free=20G -q threaded.q -pe thread 24 -N quast.LG -cwd -V
```

### Alignment of contigs to D. ananassae polytene map <a name="dana.map"></a>  
*Initial BLAST against polished contigs to determine contig orientation*  
```
makeblastdb -in contigs.fasta -out contigs.fasta -dbtype nucl -parse_seqids  
blastn -query dana.polytene.map.loci.fasta -db contigs.fasta -max_target_seqs 10 -max_hsps 10 -outfmt 6 > initial.blast.out  
```
*Extract chromosome contigs, reverse complement if necessary* 
```
samtools faidx contigs.fasta fwd.contig.name >> chr.contigs.fasta  
samtools faidx -i contigs.fasta rev.contig.name >> chr.contigs.fasta  
```
*Second round of BLAST using chromosome contigs and custom output columns*  
```
makeblastdb -in chr.contigs.fasta -out chr.contigs.fasta -dbtype nucl -parse_seqids  
blastn -query dana.polytene.map.loci.fasta -db chr.contigs.fasta -max_target_seqs 1 -max_hsps 1 -outfmt "6 qseqid sseqid pident length sstart send evalue slen" > final.blast.out  
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
