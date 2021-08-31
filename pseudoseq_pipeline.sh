#!/bin/bash

#$ -S /bin/bash
#$ -N preGWAS_pseudo
#$ -l h_vmem=25G
#$ -l tmem=25G
#$ -l h_rt=20:0:0
#$ -j y
#$ -t 1-3432
#$ -l tscratch=15G
#$ -o /

hostname
date

WORK_DIR=.;
REF=$WORK_DIR/reference/Mycobacterium_tuberculosis_h37rv.ASM19595v2.dna_sm.toplevel.fa;

SAMPLE_LIST=$WORK_DIR/sample_list.txt;

ADAPTERS=$WORK_DIR/software/trimmomatic/trimmomatic.fa;
MASK_FILE=$WORK_DIR/reference/mtb.h37rv.hypervariable.bed;
DR_FILE=$WORK_DIR/reference/tb_profiler_db.complete.bed;

GATK=$WORK_DIR/software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar
BEDTOOLS=$WORK_DIR/software/bedtools-2.25.0/bin/bedtools;
PICARD=$WORK_DIR/software/picard-2.21.9/picard.jar;

IN_DIR=$WORK_DIR/fastq

OUT_DIR=$WORK_DIR/vcf;
mkdir -p $OUT_DIR;

ALN_DIR=$WORK_DIR/pseudosequence;
mkdir -p $ALN_DIR;

CONT_DIR=$WORK_DIR/contigs
mkdir -p $CONT_DIR


TEMP_DIR=$WORK_DIR/tmp


###################
#    GETA META    #
###################

SAMPLE=$(awk "NR==${SGE_TASK_ID}" $SAMPLE_LIST)

echo "Analyzing sample ${SAMPLE}"


###################
#    GET FASTQ    #
###################

FWD=$IN_DIR/${SAMPLE}_1.fastq.gz;
REV=$IN_DIR/${SAMPLE}_2.fastq.gz;

trimmomatic PE $FWD $REV $TEMP_DIR/${SAMPLE}_R1.fastq.gz /dev/null $TEMP_DIR/${SAMPLE}_R2.fastq.gz /dev/null -phred33 ILLUMINACLIP:$ADAPTERS:1:30:11 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:10 AVGQUAL:20

FWD=$TEMP_DIR/${SAMPLE}_R1.fastq.gz;
REV=$TEMP_DIR/${SAMPLE}_R2.fastq.gz;


##################
# First assembly #
##################


echo "Assembly for $SAMPLE using SPADES"

KMER=21,33,45,55,65,75,81,101,111,121

echo "Running spades for $SAMPLE"
spades.py -o $TEMP_DIR -1 $FWD -2 $REV --isolate -t 2 -m 250 -k $KMER;

echo "Moving contigs from $TEMP_DIR to $CONT_DIR"
scp $TEMP_DIR/contigs.fasta $CONT_DIR/$SAMPLE.contigs.fasta;


###################
# Mapping contigs #
###################

CONTIGS=$TEMP_DIR/$SAMPLE.contigs.fasta

echo "Mapping contigs to H37Rv for $SAMPLE"

minimap2 -ax asm20 -R "@RG\tID:$SAMPLE\tPU:$SAMPLE\tPL:ILLUMINA\tCN:Sanger\tSM:$SAMPLE" $REF $CONTIGS | samtools view -bS - | samtools sort -@ 4 -T $TEMP_DIR/$SAMPLE.temp -O bam -o $TEMP_DIR/$SAMPLE.algn_contig.bam;


#############################
# Sample-specific reference #
#############################

echo "Creating $SAMPLE specific reference with pseudosequence"
# Creating pseudosequence for sample specific reference
# To help mapping:
  # Failing or low coverage areas get reference base
  # Deletions and gaps get reference base
bcftools mpileup -a AD,ADF,ADR,DP,SP -A -x -E -pm3 -d 1000 -L 1000 -f $REF -q 0 -Q 0 -Ou $TEMP_DIR/$SAMPLE.algn_contig.bam | bcftools call -mv -M -A -f GQ,GP -Ov | /SAN/ballouxlab/tb_lima/dev/vcf2pseudoseq.py -v - -o $TEMP_DIR/$SAMPLE.pseudo.contig.fa -H max --fail_as_ref -w -r $REF;

# Get Chromosome into header
sed -i "s/>/&Chromosome /" $TEMP_DIR/$SAMPLE.pseudo.contig.fa;

# Just in case, change any non-IUPAC into N
cat <(head -n1 $TEMP_DIR/$SAMPLE.pseudo.contig.fa) <(tail -n+2 $TEMP_DIR/$SAMPLE.pseudo.contig.fa | tr 'MRWSYKVHDB-' 'N') > $TEMP_DIR/$SAMPLE.pseudo.contig.fa.temp;
mv $TEMP_DIR/$SAMPLE.pseudo.contig.fa.temp $TEMP_DIR/$SAMPLE.pseudo.contig.fa;


##################################
# Mapp reads to sample-reference #
##################################

echo "Mapping reads into sample-specific reference"

bwa index $TEMP_DIR/$SAMPLE.pseudo.contig.fa;
java -jar $PICARD CreateSequenceDictionary R=$TEMP_DIR/$SAMPLE.pseudo.contig.fa O=$TEMP_DIR/$SAMPLE.pseudo.contig.dict
samtools faidx $TEMP_DIR/$SAMPLE.pseudo.contig.fa;

HEADER=$(zcat $FWD | head -n1 | cut -d' ' -f2-);

if [ $(echo $HEADER | tr ':' '\n' | wc -l) -eq 5 ]; then
  FLOWCELL=$(echo $HEADER | cut -d' ' -f2 | cut -d':' -f1);
  LANE=$(echo $HEADER | cut -d' ' -f2 | cut -d':' -f2);

  rg_id=$FLOWCELL.$LANE;
  pu=$FLOWCELL.$LANE.$SAMPLE;
  pl=ILLUMINA;
  sm=$SAMPLE;
else
  FLOWCELL=$(echo $HEADER | cut -d' ' -f2 | cut -d':' -f3);
  LANE=$(echo $HEADER | cut -d' ' -f2 | cut -d':' -f4);

  rg_id=$FLOWCELL.$LANE;
  pu=$FLOWCELL.$LANE.$SAMPLE;
  pl=ILLUMINA;
  sm=$SAMPLE;
fi



echo "Aligning raw reads back to $SAMPLE specific reference"
bwa mem -M -t 2 -v 1 -R "@RG\tID:$rg_id\tPU:$pu\tPL:$pl\tSM:$sm" $TEMP_DIR/$SAMPLE.pseudo.contig.fa $FWD $REV | samtools view -bS - | samtools sort -@ 4 -T $TEMP_DIR/$SAMPLE.temp -O bam -o $TEMP_DIR/$SAMPLE.toContig.bam;

echo "Marking duplicates"
java -Xms4G -Xmx6G -jar $PICARD MarkDuplicates I=$TEMP_DIR/$SAMPLE.toContig.bam O=$TEMP_DIR/$SAMPLE.toContig.markDup.bam M=$TEMP_DIR/$SAMPLE.marked_dup_metrics.txt TAGGING_POLICY=All MAX_RECORDS_IN_RAM=250000;

echo "Indel re-alignment"
java -jar $GATK -T RealignerTargetCreator -nt 4 -R $TEMP_DIR/$SAMPLE.pseudo.contig.fa -I $TEMP_DIR/$SAMPLE.toContig.markDup.bam -o $TEMP_DIR/$SAMPLE.IndelRealigner.intervals;
java -jar $GATK -T IndelRealigner -R $TEMP_DIR/$SAMPLE.pseudo.contig.fa -I $TEMP_DIR/$SAMPLE.toContig.markDup.bam -targetIntervals $TEMP_DIR/$SAMPLE.IndelRealigner.intervals -o $TEMP_DIR/$SAMPLE.toContig.markDup.reAlign.bam;

samtools index $TEMP_DIR/$SAMPLE.toContig.markDup.reAlign.bam;


####################
# Calling variants #
####################


echo "Calling variants"

bcftools mpileup -a AD,INFO/AD,ADF,INFO/ADF,ADR,INFO/ADR,DP,SP -A -x -E -pm3 -d 10000 -L 10000 -f $TEMP_DIR/$SAMPLE.pseudo.contig.fa -q 0 -Q 0 -Ou $TEMP_DIR/$SAMPLE.toContig.markDup.reAlign.bam | bcftools call -m -M -A -f GQ,GP -Oz | bcftools norm -f $TEMP_DIR/$SAMPLE.pseudo.contig.fa -Oz  | bcftools filter -m+ -s'MinMQ' -e 'INFO/MQ < 20' | bcftools filter -m+ -s'QUAL' -e 'QUAL < 20' | bcftools filter -m+ -s'PosBias' -e 'INFO/RPB < 0.001' | bcftools filter -m+ -s'MqBias' -e 'INFO/MQB < 0.001' | bcftools filter -m+ -s'StrandBias' -e 'FMT/SP > 50' | bcftools filter -m+ -s'MinDP' -e "INFO/DP < 20" -Oz -o $TEMP_DIR/$SAMPLE.sfilt.vcf.gz;
tabix -p vcf $TEMP_DIR/$SAMPLE.sfilt.vcf.gz;

echo "Adding FT tag"
$WORK_DIR/dev/addFT.py $TEMP_DIR/$SAMPLE.sfilt.vcf.gz $TEMP_DIR/$SAMPLE.sfilt.withFT.vcf;

bcftools view $TEMP_DIR/$SAMPLE.sfilt.withFT.vcf -Oz -o $TEMP_DIR/$SAMPLE.vcf.gz;

echo "Moving output to final folder"
scp $TEMP_DIR/$SAMPLE.vcf.gz $OUT_DIR/$SAMPLE.vcf.gz

echo "DONE CALLING VARIANTS"


####################
#  Pseudosequence  #
####################

VCF_FILE=$TEMP_DIR/$SAMPLE.vcf.gz;

echo "Get pseudosequence for sample $SAMPLE"

$WORK_DIR/dev/vcf2pseudoseq.py -v $VCF_FILE -o $TEMP_DIR/$SAMPLE.fa -H iupac --fill_with_N -d -w -r $REF -f "FORMAT/FT == 0" -f "FILTER != PASS"

echo "Mask for sample $SAMPLE"

# Remove 100bp around hypervariable sites
bedtools maskfasta -fi $TEMP_DIR/$SAMPLE.fa -bed <(sed "s/Chromosome/$SAMPLE/" $MASK_FILE | awk 'BEGIN{FS=OFS="\t"} $2-=101, $3+=100') -fo $TEMP_DIR/$SAMPLE.masked.fa;

# Remove drug resistance genes
bedtools maskfasta -fi $TEMP_DIR/$SAMPLE.fa -bed <(sed "s/Chromosome/$SAMPLE/" $DR_FILE | awk 'BEGIN{FS=OFS="\t"} $2-=101, $3+=100') -fo $TEMP_DIR/$SAMPLE.noDR.fa;

# Remove drug resistance genes + hypervariable
bedtools maskfasta -fi $TEMP_DIR/$SAMPLE.masked.fa -bed <(sed "s/Chromosome/$SAMPLE/" $DR_FILE | awk 'BEGIN{FS=OFS="\t"} $2-=101, $3+=100') -fo $TEMP_DIR/$SAMPLE.masked.noDR.fa;


echo "Moving output files for sample $SAMPLE"

mv $TEMP_DIR/$SAMPLE.fa $ALN_DIR;
mv $TEMP_DIR/$SAMPLE.masked.fa $ALN_DIR;
mv $TEMP_DIR/$SAMPLE.noDR.fa $ALN_DIR;
mv $TEMP_DIR/$SAMPLE.masked.noDR.fa $ALN_DIR;

date
