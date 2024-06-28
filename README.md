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


