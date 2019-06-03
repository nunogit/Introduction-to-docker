#!/bin/bash

# This contains all the lines of script needed to generate ASE variant counts using
#  the ASEReadCounter function of GATK 3.8.1
# Newer versions of ASEReadCounter (GATK >4.0) apparently stopped working. 
#  Hence until fixed, will rely on older version (November 2017, so not that old).

# 0   Preliminaries
#     ORIGINAL FASTA REQUIRED TO WHICH BAM FILES ARE ALLIGNED
#     Needs to be converted to 2bit (for GC correction) and a DICT file
#      needs to be created (for GATK):
#./faToTwoBit genome.fa genome.2bit
#java -jar picard.jar CreateSequenceDictionary REFERENCE=genome.fa OUTPUT=genome.dict

# Need to change VCF file to match chromosome notation of BAM and FASTA
#awk '{gsub(/^chrM/,"chrMT"); print}' 1000G_phase1.snps.high_confidence.hg19.vcf | \
#	awk '{gsub(/^chr/,""); print}' > 1000G_phase1.snps.high_confidence.chrRecode.hg19.vcf
#bgzip 1000G_phase1.snps.high_confidence.chrRecode.hg19.vcf
#tabix -p vcf 1000G_phase1.snps.high_confidence.chrRecode.hg19.vcf.gz
	
# Start of loop:
for f in $PWD/*.bam ; do
  FILE=`basename ${f%%.*}`
  
# 1.  GC correction: (effectiveGenomeSize: https://goo.gl/x9VuPg)
#     Takes about 1-3 hours in total depending on BAM size
computeGCBias -b ${FILE}.bam --effectiveGenomeSize 2736124973 \
 -g genome.2bit --fragmentLength 75 --GCbiasFrequenciesFile gc.${FILE}.txt \
 --numberOfProcessors 10
correctGCBias -b ${FILE}.bam --effectiveGenomeSize 2736124973 \
 -g genome.2bit --GCbiasFrequenciesFile gc.${FILE}.txt -o ${FILE}.gc.bam \
 --numberOfProcessors 10
rm -f gc.${FILE}.txt

# 2.a Reorder BAM file to match exactly with FASTA reference file contigs. This is ok as long as
#     the main chromosomes (i.e. autosomes and M, X, Y) contigs are identical between the 
#     FASTA reference file and the BAM files. 
java -jar picard.jar ReorderSam INPUT=${FILE}.gc.bam \
 OUTPUT=${FILE}.r.bam CREATE_INDEX=true REFERENCE=genome.fa \
 ALLOW_INCOMPLETE_DICT_CONCORDANCE=true ALLOW_CONTIG_LENGTH_DISCORDANCE=true
rm -f ${FILE}.gc.ba*
 
 # The error it gives is due to more contigs in FASTA file compared to BAM.
 #   See: https://goo.gl/YE1k6D
 
# 3.  Apparently the below is needed to generate additional read group characteristics
#     and generate a new index file. The former is required by GATK:
java -jar picard.jar AddOrReplaceReadGroups \
       I=${FILE}.r.bam \
       O=${FILE}.rg.bam \
       RGID=dcm \
       RGLB=dcm \
       RGPL=illumina \
       RGPU=dcm \
       RGSM=dcm \
       CREATE_INDEX=true
rm -f ${FILE}.r.ba*

# 4. ASE variant calling:

# Using GATK 3.8.1: (45 minutes)
java -jar GenomeAnalysisTK.jar -R genome.fa -T ASEReadCounter \
 -I ${FILE}.rg.bam -o $PWD/ASE_results/${FILE}_ASEresults.csv -sites \
 1000G_phase1.snps.high_confidence.chrRecode.hg19.vcf.gz \
 --minDepthOfNonFilteredBase 10 --minMappingQuality 20 --minBaseQuality 5 \
 -U ALLOW_N_CIGAR_READS --outputFormat CSV
rm -f ${FILE}.rg.ba*

done

echo $SECONDS

# END
