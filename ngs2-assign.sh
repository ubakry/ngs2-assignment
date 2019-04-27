###########################################
#  - NGS 2 Course - Assignment            #
#  - Bash Script                          #
#  - April 26,2019                        #
#  - Copyright: Usama Bakry               #
#  - Nile University                      #
###########################################
#!/bin/bash

# Downloading reference from gencode website
echo "[   PROCESS   ]     Downloading reference from gencode website..."
mkdir -p ../ref/
wget -P ../ref/ http://genomedata.org/rnaseq-tutorial/fasta/GRCh38/chr22_with_ERCC92.fa
echo "[     OK      ]      Reference is downloaded."

# Exploring the samples names
ls -tral ../data

# Renaming the samples names
mv ../data/SRR8797509_1.part_001.part_001.fastq.gz ../data/S1_L001_R1_001.fastq.gz
mv ../data/SRR8797509_2.part_001.part_001.fastq.gz ../data/S1_L001_R2_001.fastq.gz

mv ../data/shuffled_SRR8797509_1.part_001.part_001.fastq.gz ../data/S2_L001_R1_001.fastq.gz
mv ../data/shuffled_SRR8797509_2.part_001.part_001.fastq.gz ../data/S2_L001_R2_001.fastq.gz

# Indexing the reference genome
echo "[   PROCESS   ]     Indexing the reference genome..."
GENOME_DIR="../ref/star-index/"
mkdir -p $GENOME_DIR && cd $GENOME_DIR
ln -s ../chr22_with_ERCC92.fa .
cd -
STAR --runMode genomeGenerate --genomeDir $GENOME_DIR --genomeFastaFiles ../ref/chr22_with_ERCC92.fa  --runThreadN 4
echo "[     OK      ] Indexing is done" 

# Unzipping samples fastq files
cd ../data/
gunzip -k *.fastq.gz
cd -

# Add [Read group information] and align all reads
PASS_1_DIR=../GATK_results/star_res/1pass
mkdir -p $PASS_1_DIR
for R1 in ../data/*_R1_001.fastq;do
    mkdir -p $PASS_1_DIR/$(basename $R1 _R1_001.fastq)
    cd $PASS_1_DIR/$(basename $R1 _R1_001.fastq)
    R2=$(echo $R1 | sed 's/_R1_/_R2_/')
    echo $R1 $R2
    STAR --genomeDir ../../../../ref/star-index/ --readFilesIn ../../../$R1 ../../../$R2 --runThreadN 4
    cd -
done

# Indexing the reference genome
echo "[   PROCESS   ]     Indexing the reference genome..."
for S in {1..2};do
    GENOME_DIR="../ref/star-index-2/S"$S
    mkdir -p $GENOME_DIR && cd $GENOME_DIR
    ln -s ../../chr22_with_ERCC92.fa .
    cd -
    STAR --runMode genomeGenerate --genomeDir $GENOME_DIR --genomeFastaFiles ../ref/chr22_with_ERCC92.fa --sjdbFileChrStartEnd ../GATK_results/star_res/1pass/S${S}_L001/SJ.out.tab --sjdbOverhang 75 --runThreadN 4
    echo "[     OK      ] Indexing is done for sample $S" 
done

# Add [Read group information] and align all reads
PASS_2_DIR=../GATK_results/star_res/2pass
mkdir -p $PASS_2_DIR
for R1 in ../data/*_R1_001.fastq;do
    mkdir -p $PASS_2_DIR/$(basename $R1 _R1_001.fastq)
    cd $PASS_2_DIR/$(basename $R1 _R1_001.fastq)
    R2=$(echo $R1 | sed 's/_R1_/_R2_/')
    echo $R1 $R2
    echo ../../../../ref/star-index-2/$(basename $R1 | cut -d"_" -f1)
    STAR --genomeDir ../../../../ref/star-index-2/$(basename $R1 | cut -d"_" -f1) --readFilesIn ../../../$R1 ../../../$R2 --runThreadN 4
    cd -
done


# Add [Read group information]
mkdir -p ../GATK_results/picard_res
for i in S1_L001 S2_L001;do
    SM=$(basename $i | cut -d"_" -f1)                                          ##sample ID
    LB=$i                                        ##library ID
    PL="Illumina"                                                           ##platform (e.g. illumina, solid)
    RGID=$(cat ../data/${i}_R1_001.fastq | head -n1 | sed 's/ /_/g' | cut -d "_" -f1)       ##read group identifier 
    PU=$RGID.$LB                                                            ##Platform Unit
    echo -e "@RG\tID:$RGID\tSM:$SM\tPL:$PL\tLB:$LB\tPU:$PU"

    picard AddOrReplaceReadGroups I=../GATK_results/star_res/2pass/$i/Aligned.out.sam O=../GATK_results/picard_res/${SM}_rg_added_sorted.bam SO=coordinate RGID=$RGID RGLB=$LB RGPL=$PL RGPU=$PU RGSM=$SM 

    picard MarkDuplicates I=../GATK_results/picard_res/${SM}_rg_added_sorted.bam O=../GATK_results/picard_res/${SM}_dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=../GATK_results/picard_res/${SM}_output.metrics 
done

cd ../ref
samtools fqidx chr22_with_ERCC92.fa 
gatk CreateSequenceDictionary -R chr22_with_ERCC92.fa -O chr22_with_ERCC92.dict
cd -

# Split'N'Trim and reassign mapping qualities
mkdir -p ../GATK_results/gatk_res
for i in S1 S2;do
    gatk SplitNCigarReads -R ../ref/chr22_with_ERCC92.fa -I ../GATK_results/picard_res/${i}_dedupped.bam -O ../GATK_results/gatk_res/${i}_split.bam
done

wget -P ../ref ftp://ftp.ensembl.org/pub/release-96/variation/vcf/homo_sapiens/homo_sapiens-chr22.vcf.gz

cd ../ref
gunzip -k homo_sapiens-chr22.vcf.gz 
grep "^#" homo_sapiens-chr22.vcf > chr22.vcf
grep "^22" homo_sapiens-chr22.vcf | sed 's/^22/chr22/' >> chr22.vcf
gatk IndexFeatureFile -F chr22.vcf
cd -

for sample in S1 S2;do
  name=${sample%.dedup.bam}

  gatk --java-options "-Xmx2G" BaseRecalibrator \
-R ../refdog_chr5.fa -I $sample --known-sites ../ref/chr22.vcf \
-O $name.report

  gatk --java-options "-Xmx2G" ApplyBQSR \
-R dog_chr5.fa -I $sample -bqsr $name.report \
-O $name.bqsr.bam --add-output-sam-program-record --emit-original-quals
done

