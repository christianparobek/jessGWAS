### This is a BWA alignment pipeline
### For Pf whole-genome assembly
### Adapted September 2014
### Features:
###	Paired-end alignments using BWA-MEM    
###	Variant-calling using GATK


##########################################################################
###################### REFERENCE GENOME PREPARATION ######################
##########################################################################

## INDEX REFERENCE SEQUENCE FOR BWA
#bwa index PlasmoDB-9.3_Pfalciparum3D7_Genome.fasta

## INDEX REFERENCE SEQUENCE FOR BOWTIE2
#bowtie2-build PlasmoDB-9.3_Pfalciparum3D7_Genome.fasta Pf3D7_9.3

## INDEX REFERENCE SEQUENCE FOR SAMTOOLS... necessary for the mpileup step
#samtools faidx PvSal1_v10.0.fasta

## INDEX REFERENCE SEQUENCE FOR GATK
#java -jar /nas02/apps/picard-1.88/picard-tools-1.88/CreateSequenceDictionary.jar R=index PlasmoDB-9.3_Pfalciparum3D7_Genome.fasta O=index PlasmoDB-9.3_Pfalciparum3D7_Genome.dict

##########################################################################
###################### BOWTIE2 ALIGNMENT & CLEANING ######################
##########################################################################

#reads=/proj/julianog/sequence_reads/htsf_seq_backups/2014_09_19_Jess_PPQ_Pools/trimmed
#ref=/proj/julianog/refs/Pf3D7_v9.3
#masks=/proj/julianog/refs/Pf3D7_v9.3/Pf3D7_9.3/masks
#picard=/nas02/apps/picard-1.88/picard-tools-1.88


#for name in `cat filenames_trimmed.txt`
#do

#### ALIGN PAIRED-END READS WITH BOWTIE2
#bowtie2 --threads 8 \
#	-x $ref/Pf3D7_9.3 \
#	-1 $reads/$name\_R1.fastq.gz \
#	-2 $reads/$name\_R2.fastq.gz \
#	-S aln_bt2/$name.sam \
#	--maxins 750 \
#	--rg-id $name --rg PL:illumina --rg LB:$name --rg SM:$name

### SORT, AND COMPRESS SAM FILES INTO BAM
#java -jar /nas02/apps/picard-1.88/picard-tools-1.88/SortSam.jar I=aln_bt2/$name.sam O=aln_bt2/$name.sorted.bam SO=coordinate

### MARK DUPLICATES
#java -jar /nas02/apps/picard-1.88/picard-tools-1.88/MarkDuplicates.jar I=aln_bt2/$name.sorted.bam O=aln_bt2/$name.dedup.bam METRICS_FILE=aln_bt2/$name.dedup.metrics REMOVE_DUPLICATES=True

### INDEX BAM FILE PRIOR TO REALIGNMENT
#java -jar /nas02/apps/picard-1.88/picard-tools-1.88/BuildBamIndex.jar INPUT=aln_bt2/$name.dedup.bam

### IDENTIFY WHAT REGIONS NEED TO BE REALIGNED 
#java -jar /nas02/apps/biojars-1.0/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $ref/PlasmoDB-9.3_Pfalciparum3D7_Genome.fasta -I aln_bt2/$name.dedup.bam -o aln_bt2/$name.realigner.intervals -nt 8

### PERFORM THE ACTUAL REALIGNMENT
#java -jar /nas02/apps/biojars-1.0/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T IndelRealigner -R $ref/PlasmoDB-9.3_Pfalciparum3D7_Genome.fasta -I aln_bt2/$name.dedup.bam -targetIntervals aln_bt2/$name.realigner.intervals -o aln_bt2/$name.realn.bam

#done

##########################################################################
###################### BWA-MEM ALIGNMENT & CLEANING ######################
##########################################################################

reads=/proj/julianog/sequence_reads/htsf_seq_backups/2014_09_19_Jess_PPQ_Pools
ref=/proj/julianog/refs/Pf3D7_v9.3
masks=/proj/julianog/refs/Pf3D7_v9.3/Pf3D7_9.3/masks
picard=/nas02/apps/picard-1.88/picard-tools-1.88


for name in `cat filenames.txt`
do

## ALIGN PAIRED-END READS WITH BWA-MEM
bwa mem -M \
	-t 8 \
	-v 2 \
	-R "@RG\tID:$name\tPL:illumina\tLB:$name\tSM:$name" \
	-k 30 \
	$ref/PlasmoDB-9.3_Pfalciparum3D7_Genome.fasta \
	$reads/$name\_R1\_001.fastq.gz \
	$reads/$name\_R2\_001.fastq.gz \
	> aln_bwa/$name.sam
		# -M marks shorter split hits as secondary
		# -t indicates number of threads
		# -v 2 is verbosity ... warnings and errors only
		# -k 30 sets min seed length to 30

## SORT, AND COMPRESS SAM FILES INTO BAM
java -jar /nas02/apps/picard-1.88/picard-tools-1.88/SortSam.jar I=aln_bwa/$name.sam O=aln_bwa/$name.sorted.bam SO=coordinate

## MARK DUPLICATES
java -jar /nas02/apps/picard-1.88/picard-tools-1.88/MarkDuplicates.jar I=aln_bwa/$name.sorted.bam O=aln_bwa/$name.dedup.bam METRICS_FILE=aln_bwa/$name.dedup.metrics REMOVE_DUPLICATES=True

## INDEX BAM FILE PRIOR TO REALIGNMENT
java -jar /nas02/apps/picard-1.88/picard-tools-1.88/BuildBamIndex.jar INPUT=aln_bwa/$name.dedup.bam

## IDENTIFY WHAT REGIONS NEED TO BE REALIGNED 
java -jar /nas02/apps/biojars-1.0/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $ref/PlasmoDB-9.3_Pfalciparum3D7_Genome.fasta -I aln_bwa/$name.dedup.bam -o aln_bwa/$name.realigner.intervals -nt 8

## PERFORM THE ACTUAL REALIGNMENT
java -jar /nas02/apps/biojars-1.0/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T IndelRealigner -R $ref/PlasmoDB-9.3_Pfalciparum3D7_Genome.fasta -I aln_bwa/$name.dedup.bam -targetIntervals aln_bwa/$name.realigner.intervals -o aln_bwa/$name.realn.bam

done

##########################################################################
############################## POPOOLATION2 ##############################
##########################################################################

## GENERATE MPILEUP
#samtools mpileup -B -Q 0 -f /proj/julianog/refs/Pf3D7_v9.3/PlasmoDB-9.3_Pfalciparum3D7_Genome.fasta aln/high.realn.bam aln/low.realn.bam -l intervals/intervals2include.bed > poo2/hi_lo.mpileup

## CONVERT MPILEUP TO SYNC FORMAT
#java -jar ~/bin/popoolation2_1201/mpileup2sync.jar --input poo2/hi_lo.mpileup --output poo2/hi_lo.sync --fastq-type sanger --min-qual 20 --threads 8

## SUBSAMPLING
#perl ~/bin/popoolation2_1201/subsample-synchronized.pl --input poo2/hi_lo.sync --output poo2/hi_lo_subsampled.sync --target-coverage 40 --max-coverage 250 --method withreplace

####################################################
############### SNPWISE CALCULATIONS ###############
####################################################

## CALCULATE SNPWISE FISHER'S EXACT TEST
#perl ~/bin/popoolation2_1201/fisher-test.pl --input poo2/hi_lo_subsampled.sync --output poo2/hi_lo_snp.fet --min-count 2 --min-coverage 4 --max-coverage 250 --suppress-noninformative

## CONVERT SNPWISE FET TO IGV FORMAT
#perl ~/bin/popoolation2_1201/export/pwc2igv.pl --input poo2/hi_lo_snp.fet --output poo2/hi_lo_snp_fet.igv

## CALCULATE SNPWISE FST
#perl ~/bin/popoolation2_1201/fst-sliding.pl --window-size 1 --step-size 1 --suppress-noninformative --input poo2/hi_lo_subsampled.sync --min-covered-fraction 1.0 --min-coverage 15 --max-coverage 40 --min-count 3 --output poo2/hi_lo_snp.fst --pool-size 40

## CONVERT SNPWISE FST FILE TO IGV FORMAT
#perl ~/bin/popoolation2_1201/export/pwc2igv.pl --input poo2/hi_lo_snp.fst --output poo2/hi_lo_snp_fst.igv

## CALCULATE FST BY WINDOW
#perl ~/bin/popoolation2_1201/fst-sliding.pl --window-size 1000 --step-size 500 --suppress-noninformative --input poo2/hi_lo_subsampled.sync --min-covered-fraction 1.0 --min-coverage 15 --max-coverage 40 --min-count 3 --output poo2/hi_lo_snp.fst --pool-size 40

####################################################
############### GENEWISE CALCULATIONS ##############
####################################################

## CREATE GENEWISE SYNC FILE
#perl ~/bin/popoolation2_1201/create-genewise-sync.pl --input poo2/hi_lo_subsampled.sync --gtf /proj/julianog/refs/Pf3D7_v9.3/PlasmoDB-9.3_Pfalciparum3D7.gtf --output poo2/hi_lo_gene.sync

## CALCULATE GENEWISE FISHER'S EXACT
#perl ~/bin/popoolation2_1201/fisher-test.pl --input poo2/hi_lo_gene.sync --output poo2/hi_lo_gene.fet --min-count 2 --min-coverage 4 --max-coverage 250 --suppress-noninformative

## CONVERT GENEWISE FET TO IGV FORMAT
#perl ~/bin/popoolation2_1201/export/pwc2igv.pl --input poo2/hi_lo_gene.fet --output poo2/hi_lo_gene_fet.igv

## CALCULATE GENEWISE FST
#perl ~/bin/popoolation2_1201/fst-sliding.pl --window-size 1 --step-size 1 --suppress-noninformative --input poo2/hi_lo_gene.sync --min-covered-fraction 1.0 --min-coverage 15 --max-coverage 40 --min-count 3 --output poo2/hi_lo_gene.fst --pool-size 40

## CONVERT GENEWISE FST FILE TO IGV FORMAT
#perl ~/bin/popoolation2_1201/export/pwc2igv.pl --input poo2/hi_lo_gene.fst --output poo2/hi_lo_gene_fst.igv

##########################################################################
############################## EXTRA TOOLS ###############################
##########################################################################

## CALCULATE COVERAGE
#bedtools genomecov -ibam aln/$name.sorted.bam -max 10 | grep genome > $name.cov
