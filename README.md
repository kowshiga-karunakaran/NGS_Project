# NGS_Project
Repository for NGS project done during Nyberman internship

# install sra toolkit in Linux


# create directory and safe data
mkdir Project_SRA_data

# install FASTQC
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
unzip fastqc_v0.12.1.zip

# Quality check of sra data - FASTQC
cd Project_SRA_data
/workspace/NGS_Project/FastQC/fastqc ERR11468775_1.fastq.gz ERR11468775_2.fastq.gz ERR11468776_1.fastq.gz ERR11468776_2.fastq.gz ERR11468777_1.fastq.gz ERR11468777_2.fastq.gz 

# install FASTP
cd ..
wget http://opengene.org/fastp/fastp
chmod a+x ./fastp

# Read Trimming (poor quality reads or adaptors)
cd Project_SRA_data
/workspace/NGS_Project/fastp -i ERR11468775_1.fastq.gz -I ERR11468775_2.fastq.gz -o ERR11468775_1_trimmed.fastq.gz -O ERR11468775_2_trimmed.fastq.gz
/workspace/NGS_Project/fastp -i ERR11468776_1.fastq.gz -I ERR11468776_2.fastq.gz -o ERR11468776_1_trimmed.fastq.gz -O ERR11468776_2_trimmed.fastq.gz
/workspace/NGS_Project/fastp -i ERR11468777_1.fastq.gz -I ERR11468777_2.fastq.gz -o ERR11468777_1_trimmed.fastq.gz -O ERR11468777_2_trimmed.fastq.gz

# Downloading Reference Genome
mkdir project_ref_genoome
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr11.fa.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr12.fa.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr16.fa.gz
gunzip chr11.fa.gz
gunzip chr12.fa.gz
gunzip chr16.fa.gz

# Install bwa
git clone https://github.com/lh3/bwa.git
cd bwa; make

# build index
cd project_ref_genome
/workspace/NGS_Project/bwa/bwa index chr11.fa
/workspace/NGS_Project/bwa/bwa index chr12.fa
/workspace/NGS_Project/bwa/bwa index chr16.fa

# Read mapping
mkdir project_read_mapping
/workspace/NGS_Project/bwa/bwa mem -t 5 project_ref_genome/chr11.fa Project_SRA_data/ERR11468775_1_trimmed.fastq.gz Project_SRA_data/ERR11468775_2_trimmed.fastq.gz > project_read_mapping/P_sample1.sam
/workspace/NGS_Project/bwa/bwa mem -t 5 project_ref_genome/chr11.fa Project_SRA_data/ERR11468776_1_trimmed.fastq.gz Project_SRA_data/ERR11468776_2_trimmed.fastq.gz > project_read_mapping/P_sample2.sam
/workspace/NGS_Project/bwa/bwa mem -t 5 project_ref_genome/chr11.fa Project_SRA_data/ERR11468777_1_trimmed.fastq.gz Project_SRA_data/ERR11468777_2_trimmed.fastq.gz > project_read_mapping/P_sample3.sam

# install Docker - gatk
docker pull broadinstitute/gatk:latest
docker run -it -v $PWD:/data/ broadinstitute/gatk:latest

# sorting sam files
gatk SortSam -I project_read_mapping/P_sample1.sam -O P_bam_files/P_sample1_sorted.bam -SORT_ORDER coordinate
gatk SortSam -I project_read_mapping/P_sample2.sam -O P_bam_files/P_sample2_sorted.bam -SORT_ORDER coordinate
gatk SortSam -I project_read_mapping/P_sample3.sam -O P_bam_files/P_sample3_sorted.bam -SORT_ORDER coordinate

# Mark Duplicates
gatk MarkDuplicates I=P_bam_files/P_sample1_sorted.bam O=P_bam_files/sampple1_marked_duplicates.bam M=P_bam_files/sample1_marked_dup_metrics.txt
gatk MarkDuplicates I=P_bam_files/P_sample2_sorted.bam O=P_bam_files/sampple2_marked_duplicates.bam M=P_bam_files/sample2_marked_dup_metrics.txt
gatk MarkDuplicates I=P_bam_files/P_sample3_sorted.bam O=P_bam_files/sampple3_marked_duplicates.bam M=P_bam_files/sample3_marked_dup_metrics.txt

# BQSR
# download the file 
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
gatk IndexFeatureFile --input read_mapping/Homo_sapiens_assembly38.dbsnp138.vcf

# create ref genome index using samtools
samtools faidx project_ref_genome/chr11.fa
samtools faidx project_ref_genome/chr12.fa
samtools faidx project_ref_genome/chr16.fa
## create dict file
gatk CreateSequenceDictionary -R project_ref_ref_gemone/chr11.fa
gatk CreateSequenceDictionary -R project_ref_ref_gemone/chr12.fa
gatk CreateSequenceDictionary -R project_ref_ref_gemone/chr16.fa

# add read group info
gatk AddOrReplaceReadGroups I= P_bam_files/sample1_marked_duplicates.bam  O= P_mapping_qc/sam1_rg.sam RGID=sam1 RGLB=sam1 RGPL=ILLUMINA RGPU=unit1 RGSM=20 
gatk AddOrReplaceReadGroups I= P_bam_files/sample2_marked_duplicates.bam  O= P_mapping_qc/sam2_rg.sam RGID=sam2 RGLB=sam2 RGPL=ILLUMINA RGPU=unit1 RGSM=20 
gatk AddOrReplaceReadGroups I= P_bam_files/sample3_marked_duplicates.bam  O= P_mapping_qc/sam3_rg.sam RGID=sam3 RGLB=sam3 RGPL=ILLUMINA RGPU=unit1 RGSM=20 

# extract records from chr1 alone using bed file
samtools view -L P_bam_files/P_bed_files.bed P_bam_files/P_sample1_sorted.bam -o P_bam_files/sample1_chr11_bamfile.bam
samtools view -L P_bam_files/P_bed_files.bed P_bam_files/P_sample2_sorted.bam -o P_bam_files/sample2_chr11_bamfile.bam
samtools view -L P_bam_files/P_bed_files.bed P_bam_files/P_sample3_sorted.bam -o P_bam_files/sample3_chr11_bamfile.bam

# BQSR
gatk BaseRecalibrator -I P_mapping_qc/sam1_rg.sam -R project_ref_genome/chr11.fa --known-sites read_mapping/Homo_sapiens_assembly38.dbsnp138.vcf -O P_recal/sample1_recal_data.table
gatk BaseRecalibrator -I P_mapping_qc/sam1_rg.sam -R project_ref_genome/chr11.fa --known-sites read_mapping/Homo_sapiens_assembly38.dbsnp138.vcf -O P_recal/sample2_recal_data.table
gatk BaseRecalibrator -I P_mapping_qc/sam1_rg.sam -R project_ref_genome/chr11.fa --known-sites read_mapping/Homo_sapiens_assembly38.dbsnp138.vcf -O P_recal/sample3_recal_data.table

gatk ApplyBQSR -R project_ref_genome/chr11.fa -I P_mapping_qc/sam1_rg.sam --bqsr-recal-file P_recal/sample1_recal_data.table -O P_recal/sample1_recal.bam
gatk ApplyBQSR -R project_ref_genome/chr11.fa -I P_mapping_qc/sam2_rg.sam --bqsr-recal-file P_recal/sample2_recal_data.table -O P_recal/sample2_recal.bam
gatk ApplyBQSR -R project_ref_genome/chr11.fa -I P_mapping_qc/sam3_rg.sam --bqsr-recal-file P_recal/sample3_recal_data.table -O P_recal/sample3_recal.bam

samtools index P_recal/sample1_recal.bam
samtools index P_recal/sample1_recal.bam
samtools index P_recal/sample1_recal.bam

# Somatic variant calling - Mutect2
gatk Mutect2 -I P_recal/sample1_recal.bam -R project_ref_genome/chr11.fa -O P_somatic_variants/sample1.vcf.gz
gatk Mutect2 -I P_recal/sample2_recal.bam -R project_ref_genome/chr11.fa -O P_somatic_variants/sample2.vcf.gz
gatk Mutect2 -I P_recal/sample3_recal.bam -R project_ref_genome/chr11.fa -O P_somatic_variants/sample3.vcf.gz





