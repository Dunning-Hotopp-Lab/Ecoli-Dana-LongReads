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
11. [Alignment of contigs to D. ananassae polytene map](#dana.map)


### Prepare sequencing files for assembly <a name="ecoli.prep"></a>
**MinION LIG**
```
zcat minion.LIG.raw.fastq.gz | NanoLyse | gzip > minion.LIG.filter.fastq.gz
seqkit sample -s 13 -j 16 -p 0.17218 -o minion.LIGsample.fastq minion.LIG.filter.fastq.gz
```
**PacBio Sequel II**
```
seqkit sample -s 13 -j 16 -p 0.0212566 -o pbSequelII.sample.fastq pbSequelII.raw.fastq.gz
```
### Determine read length distributions <a name="ecoli.read"></a>
**MinION** 
```
bbtools/readlength.sh in=raw.reads.fastq.gz out=raw.reads.hist.out bin=1000 max=1000000 qin=33 
```
**PacBio**
```
bbtools/readlength.sh in=raw.reads.fastq.gz out=raw.reads.hist.out bin=1000 max=1000000 
```
**Generate histograms**
```
read_histograms.Rmd
```

### E. coli genome assembly using Canu <a name="ecoli.canu"></a>
**MinION**  
```
canu -p output.prefix -d output.dir genomeSize=4.6m gridEngineThreadsOption="-pe thread THREADS" gridEngineMemoryOption="-l mem_free=MEMORY" gridOptions="-P jdhotopp-lab -q threaded.q" -nanopore-raw raw.nanopore.reads.fastq 
```
**PacBio** 
```
canu -p output.prefix -d output.dir genomeSize=4.6m gridEngineThreadsOption="-pe thread THREADS" gridEngineMemoryOption="-l mem_free=MEMORY" gridOptions="-P jdhotopp-lab -q threaded.q" -pacbio-raw raw.pacbio.reads.fastq
```
**Polish, circularize, and rotate**  
```
pilon_iter.sh canu/assembly.contigs.fasta illumina.reads.R1.fastq illumina.reads.R2.fastq canu/assembly.trimmedReads.fasta  
circlator minimus2 pilon5.fasta circularise.fasta  
circlator fixstart --genes_fa ecoli.dnaA.DNA.fasta circularise.fasta rotated.fasta  
```
**Rename FASTA**  
```
for f in *_contigs.fasta; do awk '/^>/{print ">ecoli_contig" ++i; next}{print}' < $f > ${f%_c*}_contigs_rn.fasta; done
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
**Rename FASTA**  
```
for f in \*contigs.fasta; do awk '/^>/{print ">ecoli_contig" ++i; next}{print}' < $f > ${f%\_c\*})\_contigs_rn.fasta; done
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
mummer -mum -b -s -n -l 5000 ecoli.abyss.contigs.fasta ecoli.velvet.contigs.fasta | grep -i -P "[acgt]{5000,}" | awk '{name += 1; print ">"name"\n"substr(toupper($0), 100, length($0)-200)}' > ecoli.abyss+velvet.trusted.contigs.fasta  
```
**Build BLAST database for long read assemblies**  
```
mkdir db  
mv \*\_contigs_rn.fasta db/  
cd db  
for f in \*fasta; do makeblastdb -dbtype nucl -parse_seqids -in $f ; done  
```
**Retrieve the best BLAST hit for each trusted contig query**  
```
for f in db/\*fasta; do blastn -db $f -query /local/projects-t3/RDBKO/ecoli.abyss/ecoli.abyss+velvet.5kb+.trusted.contigs.fasta -outfmt 6 | sort -nk1,1 -k12,12gr -k11,11g -k3,3gr | sort -u -nk1,1 > ${f%\_r\*}\_trusted.blast_hits; done  
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
*MinION*  
```
minimap2 -x map-ont -t 8 ecoli.unicycler.consensus.fasta minion.reads.fastq > mapped.paf  
minimap2 -x map-ont -t 8 ecoli.unicycler.consensus.cut.fasta minion.reads.fastq > mapped.cut.paf  
```
*PacBio*  
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
**Re-map, filter reads to include only chimeras, manual inspection**  
*Map MinION*  
```
minimap2 -ax map-ont -t 8 ecoli.unicycler.consensus.fasta minion.reads.fastq | samtools sort -o sorted.bam  
```
*Map PacBio*  
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
minimap2 -ax map-ont -t 8 ecoli.unicycler.consensus.fasta minion.reads.fastq | samtools sort -o sorted.bam 
minimap2 -ax map-pb -t 8 ecoli.unicycler.consensus.fasta pb.reads.fastq | samtools sort -o sorted.bam  
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

### Ecoli DNA mod

**PacBio**
samtools index ecoli.genome.mapped.pbmm2.RSII_sorted.bam
/usr/local/packages/smrttools/install/current/bundles/smrttools/smrtcmds/bin/pbindex ecoli.genome.mapped.pbmm2.RSII_sorted.bam
echo "/usr/local/packages/smrttools/install/current/bundles/smrttools/smrtcmds/bin/ipdSummary ecoli.genome.mapped.pbmm2.RSII_sorted.bam --reference /local/projects-t3/RDBKO/ecoli.postassembly/mmap2/ecoli.genome.fasta --gff rsII.basemods.ALL.gff --csv rsII.basemods.ALL.csv --pvalue 0.001 --numWorkers 16 --identify m4C,m6A,m5C_TET" | qsub -P jdhotopp-lab -l mem_free=20G -q threaded.q -pe thread 16 -N smrt.ipdsummary -cwd

echo "/usr/local/packages/smrttools/install/current/bundles/smrttools/smrtcmds/bin/motifMaker find -f /local/projects-t3/RDBKO/ecoli.postassembly/mmap2/ecoli.genome.fasta -g rsII.basemods.ALL.gff -o rsII.motifs.csv" | qsub -P jdhotopp-lab -l mem_free=10G -cwd -N motif.find

echo "/usr/local/packages/smrttools/install/current/bundles/smrttools/smrtcmds/bin/motifMaker reprocess -f /local/projects-t3/RDBKO/ecoli.postassembly/mmap2/ecoli.genome.fasta -g rsII.basemods.ALL.gff -m rsII.motifs.csv -o rsII.motifs.gff" | qsub -P jdhotopp-lab -l mem_free=10G -cwd -N motif.reprocess

```
echo "/local/projects-t3/RDBKO/scripts/smrtlink_8.0.0.80529/smrtcmds/bin/bamsieve --percentage 2.16 --seed 13 /local/projects-t3/RDBKO/sequencing/E2348_69_2_25_19_PACBIO_DATA/RANDD_20190405_S64018_PL100122513-1_C01.subreads.bam /local/projects-t3/RDBKO/sequencing/E2348_69_2_25_19_PACBIO_DATA/RANDD_20190405_S64018_PL100122513-1_C01.300X.subset.subreads.bam" | qsub -P jdhotopp-lab -l mem_free=50G -N pb_bamsieve -cwd

/local/projects-t3/RDBKO/scripts/smrtlink_8.0.0.80529/smrtcmds/bin/dataset create --type SubreadSet --name dana.sequelII /path/to/subreadset.xml  /path/to/subreads.bam
/local/projects-t3/RDBKO/scripts/smrtlink_8.0.0.80529/smrtcmds/bin/dataset create --type ReferenceSet --name dana.chr2R /path/to/referenceset.xml /path/to/chr2R.fasta


```
echo "/local/projects-t3/RDBKO/scripts/smrtlink_8.0.0.80529/smrtcmds/bin/pbcromwell run pb_basemods -e /local/projects-t3/RDBKO/sequencing/E2348_69_2_25_19_PACBIO_DATA/RANDD_20190405_S64018_PL100122513-1_C01.subreadset.xml -e /local/projects-t3/RDBKO/ecoli.postassembly/mmap2/ecoli.unicycler.genome+plasmids.referenceset.xml -t kineticstools_compute_methyl_fraction=True -t kineticstools_identify_mods=m4C,m6A,m5C_TET -t run_find_motifs=True" | qsub -P jdhotopp-lab -l mem_free=50G -N pb_basemods -cwd


```
PATH=/local/aberdeen2rw/julie/Matt_dir/packages/miniconda3/bin:"$PATH"

qsub -P jdhotopp-lab -q threaded.q  -pe thread 32 -l mem_free=10G -N guppy -cwd /local/projects-t3/RDBKO/scripts/MinIONGuppyBasecall_v3.1.5.sh -fast5_dir /local/projects-t3/RDBKO/sequencing/RANDD_LIG_Ecoli_MS100122513/20190405_1744_MN23690_FAK57346_fe910659/fast5/ -output_dir /local/projects-t3/RDBKO/ecoli.epi/RANDD_LIG_Ecoli_MS100122513 -config dna_r9.4.1_450bps_fast.cfg

echo "multi_to_single_fast5 -i workspace -s single_fast5 -t 16" | qsub -P jdhotopp-lab -q threaded.q -pe thread 16 -l mem_free=50G -N multi_to_single_fast5 -cwd -V

echo "/local/aberdeen2rw/julie/Matt_dir/packages/miniconda3/bin/tombo resquiggle single_fast5 /local/projects-t3/RDBKO/ecoli.postassembly/mmap2/ecoli.unicycler.consensus.fasta --processes 16 --num-most-common-errors 5" | qsub -P jdhotopp-lab -l mem_free=20G -q threaded.q -pe thread 16 -N tombo.resquiggle -cwd -V

model
echo "/local/aberdeen2rw/julie/Matt_dir/packages/miniconda3/bin/tombo detect_modifications alternative_model --fast5-basedirs single_fast5 --statistics-file-basename RAPID.alt --alternate-bases 5mC 6mA --processes 16" | qsub -P jdhotopp-lab -l mem_free=20G -q threaded.q -pe thread 16 -N tombo.detect -cwd -V

denovo
echo "/local/aberdeen2rw/julie/Matt_dir/packages/miniconda3/bin/tombo detect_modifications de_novo --fast5-basedirs single_fast5 --statistics-file-basename RAPID.denovo --processes 16" | qsub -P jdhotopp-lab -l mem_free=20G -q threaded.q -pe thread 16 -N tombo.detect -cwd -V

echo "/local/aberdeen2rw/julie/Matt_dir/packages/miniconda3/bin/tombo text_output browser_files --fast5-basedirs single_fast5 --statistics-filename RAPID.denovo.tombo.stats --file-types dampened_fraction --browser-file-basename RAPID.denovo" | qsub -P jdhotopp-lab -q threaded.q -pe thread 8 -l mem_free=50G -N tombo.textoutput.browser -cwd -V

echo "/local/aberdeen2rw/julie/Matt_dir/packages/miniconda3/bin/tombo text_output signif_sequence_context --fast5-basedirs single_fast5 --statistics-filename RAPID.denovo.tombo.stats --num-regions 1000 --num-bases 50" | qsub -P jdhotopp-lab -q threaded.q -pe thread 8 -l mem_free=50G -N tombo.signif.context -cwd -V

echo -e "/usr/local/packages/meme-4.12.0/bin/meme -oc RANDD_RAPID_Ecoli.tombo.stats.dam.meme -dna -mod zoops -nmotifs 50 tombo_results.significant_regions.fasta" | qsub -P jdhotopp-lab -l mem_free=5G -N meme -cwd
```
motif-based plotting
```
echo "/local/aberdeen2rw/julie/Matt_dir/packages/miniconda3/bin/tombo plot motif_with_stats --fast5-basedirs single_fast5 --statistics-filename RAPID.denovo.tombo.stats --genome-fasta /local/projects-t3/RDBKO/ecoli.postassembly/mmap2/ecoli.genome.fasta --motif GATC --plot-standard-model --num-statistics 10000 --num-regions 1 --pdf-filename RAPID.denovo.plot.GATC.pdf" | qsub -P jdhotopp-lab -l mem_free=20G -N tombo.plot -cwd -V
echo "/local/aberdeen2rw/julie/Matt_dir/packages/miniconda3/bin/tombo plot motif_with_stats --fast5-basedirs single_fast5 --statistics-filename RAPID.denovo.tombo.stats --genome-fasta /local/projects-t3/RDBKO/ecoli.postassembly/mmap2/ecoli.genome.fasta --motif CCWGG --plot-standard-model --num-statistics 10000 --num-regions 1 --pdf-filename RAPID.denovo.plot.CCWGG.pdf" | qsub -P jdhotopp-lab -l mem_free=20G -N tombo.plot -cwd -V
```
ROC curve - needs to be done on command line
```
/local/aberdeen2rw/julie/Matt_dir/packages/miniconda3/bin/tombo plot roc --statistics-filenames RAPID.alt.5mC.tombo.stats RAPID.alt.6mA.tombo.stats RAPID.denovo.tombo.stats --motif-descriptions CCWGG:2:"dcm 5mC Alt. Model" GATC:2:"dam 6mA Alt. Model" CCWGG:2:"dcm 5mC De Novo"::GATC:2:"dam 6mA De Novo" --genome-fasta /local/projects-t3/RDBKO/ecoli.postassembly/mmap2/ecoli.unicycler.consensus.fasta
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
**KAT**
```
use kat-2.4.0
echo "kat comp -t 16 -o minion.RAPID.assembly /local/projects-t3/RDBKO/sequencing/cHI_Dana_2_15_19_ILLUMINA_DATA/RANDD_20190322_K00134_IL100123454_MX29_L004_R1.fastq /local/projects-t3/RDBKO/dana.postassembly/dana.minion.RAPID.pilon.l_contigs.fasta" | qsub -P jdhotopp-lab -l mem_free=10G -q threaded.q -pe thread 16 -N kat.comp -cwd -V
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
