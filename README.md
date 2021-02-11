# Comparison of long read sequencing technologies in interrogating bacteria and fly genomes

Eric S. Tvedte

2021-01-28

## Table of Contents
1. [Prepare sequencing files for assembly](#ecoli.prep)
2. [Read length distributions](#ecoli.read)
3. [E. coli genome assembly](#ecoli.asm)
4. [Evaluation of E. coli genome assemblies](#ecoli.eval)
5. [E. coli DNA modification analysis](#ecoli.dna.mod)  
6. [E. coli plasmid composition analysis](#ecoli.plasmid)  
7. [D. ananassae genome assembly](#dana.canu)
8. [Dana.UMIGS genome assembly](#dana.umigs)
9. [Anchoring D. ananassae assembly contigs](#dana.anchor)
10. [Evaluation of D. ananassae genome assemblies](#dana.eval)  
11. [D. ananassae DNA modification analysis](#dana.dna.mod)
12. [R scripts for data visualization](#data.vis)


### Prepare sequencing files for assembly <a name="ecoli.prep"></a>
**Basecalling R9 ONT reads**
```
guppy_basecaller --input_path fast5_dir --save_path output_dir --config dna_r9.4.1_450bps_hac.cfg --fast5_out --post_out --qscore_filtering --min_qscore 7 --records_per_fastq 10000000 -x "cuda:5 cuda:6 cuda:7" > guppy.4.2.2.log
```
**Basecalling R10 ONT reads**
```
guppy_basecaller --input_path fast5_dir --save_path output_dir --config dna_r10_450bps_hac.cfg --fast5_out --post_out --qscore_filtering --min_qscore 7 --records_per_fastq 10000000 -x "cuda:5 cuda:6 cuda:7" > guppy.4.2.2.log
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
### E. coli genome assembly <a name="ecoli.asm"></a>

**Canu ONT assembly**  
```
canu -p output.prefix -d output.dir genomeSize=4.6m corOutCoverage=1000 gridEngineThreadsOption="-pe thread THREADS" gridEngineMemoryOption="-l mem_free=MEMORY" gridOptions="-P jdhotopp-lab -q threaded.q" -nanopore ont.reads.fastq 
```
**Canu PacBio RSII/CLR assembly** 
```
canu -p output.prefix -d output.dir genomeSize=4.6m corOutCoverage=1000 gridEngineThreadsOption="-pe thread THREADS" gridEngineMemoryOption="-l mem_free=MEMORY" gridOptions="-P jdhotopp-lab -q threaded.q" -pacbio pacbio.reads.fastq
```
**HiCanu PacBio HiFi assembly** 
```
canu -p output.prefix -d output.dir genomeSize=4.6m corOutCoverage=1000 gridEngineThreadsOption="-pe thread THREADS" gridEngineMemoryOption="-l mem_free=MEMORY" gridOptions="-P jdhotopp-lab -q threaded.q" -pacbio-hifi pacbio.hifi.reads.fastq
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

### Evaluation of E. coli genome assemblies <a name="ecoli.eval"></a>   
**BUSCO: conserved genes** 
```
busco --config config.ini -i contigs.fasta -c 8 -o output_dir -l bacteria -m geno  
```
**QUAST: consensus identity, missassemblies**
```
quast.py -m 0 -t 8 ONT.RAPID.canu.raw.genome.contig.fasta ONT.RAPID.canu.polished.genome.contig.fasta ONT.RAPID.flye.genome.contig.fasta ONT.LIG.canu.raw.genome.contig.fasta ONT.LIG.canu.polished.genome.contig.fasta ONT.LIG.flye.genome.contig.fasta PB.RSII.canu.raw.genome.contig.fasta PB.RSII.canu.polished.genome.contig.fasta PB.RSII.flye.genome.contig.fasta PB.CLR.canu.raw.genome.contig.fasta PB.CLR.canu.polished.genome.contig.fasta PB.CLR.flye.genome.contig.fasta PB.HiFi.canu.raw.genome.contig.fasta PB.HiFi.canu.polished.genome.contig.fasta PB.HiFi.flye.genome.contig.fasta -r Ecoli.UMIGS.2X.fasta -o output_dir
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
### E. coli plasmid composition analysis <a name="ecoli.plasmid"></a>  
**Minimap2: map long reads to Ecoli.UMIGS**  
```
minimap2 -ax map-ont -t 8 Ecoli.UMIGS.fasta ont.reads.fastq | samtools sort -o ont.sorted.bam 
minimap2 -ax map-pb -t 8 Ecoli.UMIGS.fasta pb.reads.fastq | samtools sort -o pb.sorted.bam  
```
**BWA-MEM: map short reads to Ecoli.UMIGS**
```
bwa mem -k 23 -t 8 Ecoli.UMIGS.fasta illumina.R1.fastq illumina.R2.fastq | samtools sort -o sorted.illumina.bam  
```
**PICARD: remove Illumina duplicates** 
```
picard.jar MarkDuplicates I=sorted.illumina.bam O=dupsmarked.illumina.bam M=dupsmarked.illumina.metrics VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=TRUE AS=TRUE CREATE_INDEX=TRUE  
```
**Calculate mapping rate and sequencing depth**  
```
samtools view sorted.bam -F 4 -F 256 -F 2048 -c aln.bam # mapping rate
samtools depth -aa -m 100000000 aln.bam > aln.depth.txt # sequencing depth
```
**Count total depth across genome, pMAR2, p5217**  
```
grep ecoli.genome aln.depth.tx | awk '{total = total + $3}END{print "Total genome depth = "total}' -  # total depth of ecoli genome
grep pMAR2 aln.depth.tx | awk '{total = total + $3}END{print "Total pMAR2 depth = "total}' -  # total depth of pMAR2
grep p5217 aln.depth.tx | awk '{total = total + $3}END{print "Total p5217 depth = "total}' -  # total depth of p5217
```

### Ecoli DNA modification analysis <a name="ecoli.dna.mod"></a>  

**PacBio RS II reads**
```
samtools faidx Ecoli.UMIGS.fasta ecoli.genome > Ecoli.UMIGS.genome.fasta
pbmm2 align Ecoli.UMIGS.genome.fasta subreads.bam RSII.mapped_sorted.bam --sort -j 16 -J 8
samtools index RSII.mapped_sorted.bam
pbindex RSII.mapped_sorted.bam
ipdSummary RSII.mapped_sorted.bam --reference Ecoli.UMIGS.fasta --gff pb.basemods.gff --csv pb.basemods.csv --pvalue 0.001 --numWorkers 16 --identify m4C,m6A,m5C_TET
motifMaker find -f Ecoli.UMIGS.fasta -g pb.basemods.gff -o pb.motifs.csv
motifMaker reprocess -f Ecoli.UMIGS.fasta -g pb.basemods.gff -m pb.motifs.csv -o pb.motifs.gff
```
**PacBio Sequel II CLR reads**
```
smrtcmds/bin/dataset create --type SubreadSet --name ecoli.subreadset subreadset.xml SQII.CLR.subreads.bam
smrtcmds/bin/dataset create --type ReferenceSet --name ecoli.referenceset referenceset.xml Ecoli.UMIGS.genome.fasta
smrtcmds/bin/pbcromwell run pb_basemods -e subreadset.xml -e referenceset.xml -t kineticstools_compute_methyl_fraction=True -t kineticstools_identify_mods=m4C,m6A,m5C_TET -t run_find_motifs=True
smrtcmds/bin/pbcromwell run pb_basemods -e subreadset.xml -e referenceset.xml -t kineticstools_compute_methyl_fraction=True -t kineticstools_identify_mods=m4C,m6A,m5C_TET -t run_find_motifs=True -t motif_min_score=125
```

**ONT reads**
```
tombo preprocess annotate_raw_with_fastqs --overwrite --fast5-basedir fast5_dir --fastq-filenames ont.reads.fastq --processes 16
tombo resquiggle fast5_dir Ecoli.UMIGS.genome.fasta --processes 16 --num-most-common-errors 5
tombo detect_modifications de_novo --fast5-basedirs fast5_dir --statistics-file-basename ONT.denovo --processes 16
tombo text_output browser_files --fast5-basedirs fast5_dir --statistics-filename ONT.denovo.tombo.stats --file-types dampened_fraction --browser-file-basename ONT.denovo
tombo plot motif_with_stats --fast5-basedirs fast5_dir --statistics-filename ONT.denovo.tombo.stats --genome-fasta ecoli.unicycler.consensus.fasta --motif GATC --plot-standard-model --num-statistics 10000 --num-regions 1 --pdf-filename ONT.denovo.plot.GATC.pdf
tombo plot motif_with_stats --fast5-basedirs fast5_dir --statistics-filename ONT.denovo.tombo.stats --genome-fasta ecoli.unicycler.consensus.fasta --motif CCWGG --plot-standard-model --num-statistics 10000 --num-regions 1 --pdf-filename ONT.denovo.plot.CCWGG.pdf
tombo plot roc --statistics-filenames RAPID.5mC.tombo.stats RAPID.6mA.tombo.stats LIG.5mC.tombo.stats LIG.6mA.tombo.stats --motif-descriptions CCWGG:2:"CCWGG RAPID" GATC:2:"GATC RAPID" CCWGG:2:"CCWGG LIG"::GATC:2:"GATC LIG" --genome-fasta Ecoli.UMIGS.fasta
```

**Plot PacBio modification QV versus ONT dampened fraction**  
*PacBio*
```
grep 'm6A' pb.motifs.gff | awk '{print $1"\t"$4"\t"$6}' > pb.basemods.m6A.tsv
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

### D. ananassae genome assembly <a name="dana.canu"></a>
**Canu ONT assembly**  
```
canu -p output.prefix -d output.dir genomeSize=240m gridEngineThreadsOption="-pe thread THREADS" gridEngineMemoryOption="-l mem_free=MEMORY" gridOptions="-P jdhotopp-lab -q threaded.q" -nanopore ont.reads.fastq 
```
**Canu PacBio CLR assembly** 
```
canu -p output.prefix -d output.dir genomeSize=240m gridEngineThreadsOption="-pe thread THREADS" gridEngineMemoryOption="-l mem_free=MEMORY" gridOptions="-P jdhotopp-lab -q threaded.q" -pacbio pacbio.reads.fastq
```
**Canu PacBio ONT-CLR assembly** 
```
canu -p output.prefix -d output.dir genomeSize=240m gridEngineThreadsOption="-pe thread THREADS" gridEngineMemoryOption="-l mem_free=MEMORY" gridOptions="-P jdhotopp-lab -q threaded.q" -pacbio pacbio.reads.fastq -nanopore ont.reads.fastq
```
**HiCanu PacBio HiFi assembly** 
```
canu -p output.prefix -d output.dir genomeSize=240m gridEngineThreadsOption="-pe thread THREADS" gridEngineMemoryOption="-l mem_free=MEMORY" gridOptions="-P jdhotopp-lab -q threaded.q" -pacbio-hifi pacbio.hifi.reads.fastq
```
**Flye ONT assembly**
```
flye --asm-coverage 40 --genome-size 240m --nano-raw ont.reads.fastq -t 24 -o output_dir
```
**Flye PacBio CLR assembly**
```
flye --asm-coverage 40 --genome-size 240m --pacbio-raw pb.reads.fastq -t 24 -o output_dir
```
**Flye PacBio HiFi assembly**
```
flye --asm-coverage 40 --genome-size 240m --pacbio-hifi pb.hifi.reads.fastq -t 24 -o output_dir
```
**Polish ONT assemblies**  
```
bwa index ont.asm.fasta
bwa mem -t 16 ont.asm.fasta fwd.R1.fastq rev.R2.fastq | samtools view -@ 16 -bS - | samtools sort -@ 16 -o ShortRead.bam -
pilon-1.22.jar --genome genome.fasta --frags ShortRead.bam --changes --threads 16 --minmq 10 --fix bases --output pilon       
```
### Dana.UMIGS genome assembly <a name="dana.umigs"></a>
**Extract chromosome arm contigs**
```
samtools faidx pb.clr.canu.contigs.fasta tig00000049 > chrXL.fasta
samtools faidx pb.clr.canu.contigs.fasta tig00000036 > chrXR.fasta
samtools faidx pb.clr.flye.contigs.fasta contig_269 > chr2L.fasta
samtools faidx pb.clr.flye.contigs.fasta contig_298 > chr2R.fasta
samtools faidx pb.clr.canu.contigs.fasta tig00000025 > chr3L.fasta
samtools faidx ont.flye.contigs.fasta scaffold_346 > chr3R.fasta
cat chr* pb.clr.canu.non.chr.contigs.fasta > pb.clr.merge.fasta
```
**Merge and polish**
```
seqkit seq -m 50000 pb.hifi.contigs.fasta > pb.hifi.50kbp.contigs.fasta
quickmerge/merge_wrapper.py pb.clr.merge.fasta pb.hifi.50kbp.contigs.fasta
pbmm2 align merged.asm.fasta pb.clr.subreads.bam aln.bam --sort -j 16 -J 8 
arrow -j 16 --algorithm arrow --noEvidenceConsensusCall reference --referenceFilename merged.asm.fasta aln.bam -o arrow.polished.fasta -o arrow.polished.gff
minimap2 -t 16 -ax map-pb arrow.polished.fasta pb.hifi.reads.fastq | samtools view -@ 24 -bS - | samtools sort -@ 24 -o  LongRead.bam -
pilon-1.22.jar --genome genome.fasta --unpaired LongRead.bam --changes --threads 16 --minmq 10 --fix bases --output pilon
```
**Purge haplotigs**
```
minimap2 -t 4 -ax map-pb pilon.fasta pb.hifi.fastq.gz --secondary=no | samtools sort -m 5G -o aln.bam
purge_haplotigs hist -b aln.bam -g pilon.fasta -t 4 -d 400
purge_haplotigs cov -l 5 -m 195 -h 300 -i aln.bam.gencov -j 80 -s 80
purge_haplotigs purge -g pilon.fasta -b aln.bam -c coverage_stats.csv -d -t 8
```

## Anchoring D. ananassae assembly contigs <a name="dana.anchor"></a>
**Major chromosome arm contigs (X, 2, 3)**
```
blastn -query chrom.map.loci.fasta -db asm.fasta -max_target_seqs 5 -max_hsps 5 -outfmt "6 qseqid sseqid pident length sstart send evalue slen" > initial.blast.out #inspect output to generate list of chromosome arm sequences
seqkit grep -f chr.arm.contigs.list asm.fasta > chr.arm.contigs.fasta
blastn -query chrom.map.loci.fasta -db chr.arm.contigs.fasta -max_target_seqs 1 -max_hsps 1 -outfmt "6 qseqid sseqid pident length sstart send sstrand evalue slen" > final.blast.out
```

**Heterochromatic contigs**  
*Chromosome 4*
```
minimap2 -cx asm5 /local/projects-t3/RDBKO/dana.postassembly/Dana.UMIGS.unpurged.contigs.fasta Leung2017_chr4_scaffolds.final.fasta --secondary=no > chr4.aln.paf
awk '{print $6}' chr4.aln.paf | sort -n | uniq > chr4.candidate.list
seqkit grep -f chr4.candidate.list Dana.UMIGS.contigs.fasta > chr4.candidate.contigs.fasta 

nucmer -l 1000 --prefix chr4 Leung2017_chr4_scaffolds.final.fasta chr4.candidate.contigs.fasta
mummerplot --color --png --prefix chr4 chr4.delta #manually inspect and refine results
```
*LGT*
```
nucmer -l 1000 --prefix asm.LGT wAna.genome.fasta asm.fasta
delta-filter asm.LGT.delta -q > asm.LGT.filter
show-coords -r asm.LGT.filter -T > asm.LGT.filter.coords
tail -n +5 asm.LGT.delta.coords | awk '{print $13}' | sort -n | uniq > LGT.contigs.list #retrieves set of LGT contigs
```

*Chromosome Y*
```
seqkit sample -n 107000000  
seqkit sample -j 8 -p 0.0216 -o Ecoli.PB.HiFi.100X.fastq Ecoli.PB.HiFi.fastq.gz
bwa mem -k 23 -t 8 Ecoli.UMIGS.fasta fwd.R1.fastq rev.R2.fastq | samtools sort -o sorted.illumina.bam  
picard.jar MarkDuplicates I=sorted.illumina.bam O=dupsmarked.illumina.bam M=dupsmarked.illumina.metrics VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=TRUE AS=TRUE CREATE_INDEX=TRUE
samtools merge -@ 4 female.aln.bam female.1.aln.bam female.2.aln.bam
samtools merge -@ 4 male.aln.bam male.1.aln.bam male.2.aln.bam
samtools depth -Q 10 female.aln.bam male.aln.bam > fem.v.male.aln.out
perl Chang_Larracuente_frame_depth.pl fem.v.male.aln.out
```

### Evaluation of D. ananassae genome assemblies <a name="dana.eval"></a>  
**QUAST - full assemblies** 
```
quast.py ONT.canu.contigs.fasta ONT.flye.contigs.fasta PB.CLR.canu.contigs.fasta PB.CLR.flye.contigs.fasta PB.HiFi.canu.contigs.fasta PB.HiFi.flye.contigs.fasta Hybrid.ONT+CLR.canu.contigs.fasta Dana.Miller2018.fasta Dana.D12GC2007caf1.fasta -r Dana.UMIGS.chr.contigs.fasta -o output_dir -t 24 --large -m 0 --split-scaffolds
```
**KAT**
```
kat comp -t 8 -o output_prefix 'Illumina.R1.fastq Illumina.R2.fastq' asm.fasta 
kat plot spectra-cn -m 200 kat.mx 
```
**BUSCO**
```
busco --config config.ini -i asm.fasta -o output_dir -l arthropoda -m geno 
```
**QUAST - chromosome contigs**
```
minimap2 -cx asm5 asm.fasta Dana.UMIGS.chromosome.fasta --secondary=no > chr.aln.paf  
awk '$11>50000' chr.aln.paf | awk '{print $6}' | sort -n | uniq > chr.contig.list  
seqkit grep -f chr.contig.list asm.fasta > chr.contig.fasta  
quast.py ONT.canu.chr.contigs.fasta ONT.flye.chr.contigs.fasta PB.CLR.canu.chr.contigs.fasta PB.CLR.flye.chr.contigs.fasta PB.HiFi.canu.chr.contigs.fasta PB.HiFi.flye.chr.contigs.fasta Hybrid.ONT+CLR.canu.chr.contigs.fasta Dana.Miller2018.chr.contigs.fasta Dana.D12GC2007caf1.chr.contigs.fasta -r Dana.UMIGS.chr.contigs.fasta -o output_dir -t 24 --large -m 0 --split-scaffolds
```
**NUCmer - chromosome contigs**
```
nucmer -l 500 --maxmatch --prefix chr.align Dana.UMIGS.chromosome.fasta chr.contig.fasta
mummerplot --color --postscript --small --prefix chr.align -Q chr.contig.orientation.list -R Dana.UMIGS.chromosome.fasta
```
**Estimating LGT in D. ananassae assemblies**
```
nucmer -l 1000 --prefix firstpass wAna.genome.fasta asm.fasta
show-coords -r -T firstpass.delta > firstpass.coords
tail -n +5 finalpass.coords | awk '{print $9}' | sort -n | uniq > LGT.contigs.list
seqkit grep -f LGT.contigs.list asm.fasta > LGT.contigs.fasta
nucmer -l 1000 --prefix finalpass wAna.genome.fasta LGT.contigs.fasta
tail -n +5 finalpass.coords | awk '{print $9"\t"$3"\t"$4}' > finalpass.bed
fixbed.R finalpass.bed fixed.bed
bedtools coverage -a fixed.bed -b fixed.bed -hist | grep all | awk '{total = total + $3}END{print "LGT segment length = "total}' #estimates alignment block length of wAna to LGT
```
**D. ananassae heterochromatin sequencing bias**
```
minimap2 -t 4 -ax map-pb asm.fasta pb.fastq.gz --secondary=no | samtools sort -m 5G -o pb.aln.bam
minimap2 -t 4 -ax map-ont asm.fasta ont.fastq.gz --secondary=no | samtools sort -m 5G -o ont.aln.bam
purge_haplotigs hist -b aln.bam -g Dana.UMIGS.fasta -t 4
samtools view aln.bam -L eu.bed > eu.aln.bam
samtools view aln.bam -L het.bed > het.aln.bam
purge_haplotigs hist -b eu.aln.bam -g Dana.UMIGS.fasta -t 4
purge_haplotigs hist -b het.aln.bam -g Dana.UMIGS.fasta -t 4
```

## D. ananassae DNA modification analysis <a name="dana.dna.mod"></a>
```
tombo preprocess annotate_raw_with_fastqs --overwrite --fast5-basedir fast5_dir --fastq-filenames ont.reads.fastq --processes 16
tombo resquiggle fast5_dir Dana.UMIGS.fasta --processes 16 --num-most-common-errors 5
tombo detect_modifications de_novo --fast5-basedirs fast5_dir --statistics-file-basename ONT.denovo --processes 16
tombo text_output signif_sequence_context --fast5-basedirs fast5_dir --statistics-filename ONT.denovo.tombo.stats --num-regions 1000 --num-bases 50
meme -oc Dana.tombo.meme -dna -mod zoops -nmotifs 50 tombo.significant_regions.fasta
```

## R scripts for data visualization <a name="data.vis"></a>
R scripts and input data used in data visualization are available on Figshare

## System requirements

R scripts were run using Windows 10 x64 with RStudio v1.1.463 using this R session:
```
R version 3.5.1 (2018-07-02)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 18362)
```
