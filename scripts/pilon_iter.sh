#this script was adapted from:
#Chang and Larracuente 2019. Heterochromatin-enriched assemblies reveal the sequence and organization of the Drosophila melanogaster Y chromosome
#https://github.com/LarracuenteLab/mel.heterochromatin.Y.assembly
 

#genome in $1, shortreadR1 in $2, shortreadR2 in $3, longread in $4

mkdir -p pilon.long+short
cd pilon.long+short
cp $1 pilon0.fasta

#number of loops is hardcoded to 5 runs of pilon
for i in {1..5}; do
	if [[ $i == '1' ]]; then
		GENOME=pilon0.fasta
	else
		x=$(($i - 1))
		GENOME=pilon$x.fasta
	fi
	if [[ ! -f pilon$1.fasta ]]; then
		echo "Working on $GENOME"
		#map Illumina reads to the genome 
		#comment out based on need for short read polishing
		bwa index $GENOME $GENOME
		bwa mem -t 24 $GENOME $2 $3 |samtools view -@ 24 -bS - |samtools sort -@ 24 -o  ShortRead.bam -
		#change ont or pb based on lr 
		#comment out based on need for long read polishing
		minimap2 -t 24 -ax map-ont $GENOME $4 |samtools view -@ 24 -bS - |samtools sort -@ 24 -o  LongRead.bam -
		samtools index ShortRead.bam
		samtools index LongRead.bam
		#run pilon
		java -Xmx200g -jar /usr/local/packages/pilon-1.22/pilon-1.22.jar --genome $GENOME --frags ShortRead.bam --unpaired LongRead.bam --changes --threads 24 --minmq 10 --fix bases --output pilon$i
		rm *.bam *.bai *.sa *amb *ann *pac *bwt

	else
		echo "$GENOME has already been run, moving onto next iteration"
	fi
done
