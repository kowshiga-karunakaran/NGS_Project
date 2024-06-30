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
/home/hp/FastQC/fastqc ERR11468775_1.fastq.gz ERR11468775_2.fastq.gz ERR11468776_1.fastq.gz ERR11468776_2.fastq.gz ERR11468777_1.fastq.gz ERR11468777_2.fastq.gz 

# install FASTP
cd ..
wget http://opengene.org/fastp/fastp
chmod a+x ./fastp

# Read Trimming (poor quality reads or adaptors)
cd Project_SRA_data
/home/hp/fastp -i ERR11468775_1.fastq.gz -I ERR11468775_2.fastq.gz -o ERR11468775_1_trimmed.fastq.gz -O ERR11468775_2_trimmed.fastq.gz
/home/hp/fastp -i ERR11468776_1.fastq.gz -I ERR11468776_2.fastq.gz -o ERR11468776_1_trimmed.fastq.gz -O ERR11468776_2_trimmed.fastq.gz
/home/hp/fastp -i ERR11468777_1.fastq.gz -I ERR11468777_2.fastq.gz -o ERR11468777_1_trimmed.fastq.gz -O ERR11468777_2_trimmed.fastq.gz

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
/home/hp/bwa/bwa index chr11.fa
/home/hp/bwa/bwa index chr12.fa
/home/hp/bwa/bwa index chr16.fa

# Read mapping
mkdir project_read_mapping
/home/hp/bwa/bwa mem -t 5 project_ref_genome/chr11.fa Project_SRA_data/ERR11468775_1_trimmed.fastq.gz Project_SRA_data/ERR11468775_2_trimmed.fastq.gz > project_read_mapping/P_sample1.sam
/home/hp/bwa/bwa mem -t 5 project_ref_genome/chr11.fa Project_SRA_data/ERR11468776_1_trimmed.fastq.gz Project_SRA_data/ERR11468776_2_trimmed.fastq.gz > project_read_mapping/P_sample2.sam
/home/hp/bwa/bwa mem -t 5 project_ref_genome/chr11.fa Project_SRA_data/ERR11468777_1_trimmed.fastq.gz Project_SRA_data/ERR11468777_2_trimmed.fastq.gz > project_read_mapping/P_sample3.sam

# install Docker - gatk
docker pull broadinstitute/gatk:latest
docker run -it -v $PWD:/data/ broadinstitute/gatk:latest

# 



