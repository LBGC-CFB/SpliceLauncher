# SpliceLauncher  

---

The RNAseq pipeline tool to study the alternative splicing

> **Cite as:** SpliceLauncher: a tool for detection, annotation and relative quantification of alternative junctions from RNAseq data. *Raphaël Leman, Grégoire Davy, Antoine Rousselain, Alexandre Atkinson, Sophie Krieger*


## Repository contents  

---

* dataTest: example of input files
* refData: reference files that SpliceLauncher needs
* scripts: complementary scripts to run SpliceLauncher pipeline

## install and configure SpliceLauncher pipeline

---

The SpliceLauncher pipeline needs to start from fastq files  
The tools:

* STAR  
* samtools  
* BEDtools  

The reference files:

* the human genome assembly  
* the exons annotations  
* the transcripts information  

An example of these two last files is provide in the [refData folder](https://github.com/raphaelleman/SpliceLauncher/tree/master/refData "tittle")

### install STAR, samtools, BEDtools  

---

To use SpliceLauncher pipeline please install the three tools: STAR, samtools, BEDtools

### STAR  

#### these following instruction were from the [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf "Title")  

Get the g++ compiler for linux

    sudo apt-get update
    sudo apt-get install g++
    sudo apt-get install make

Download the [latest release](https://github.com/alexdobin/STAR/releases "Title") from and uncompress it 

    # Get latest STAR source from releases
    wget https://github.com/alexdobin/STAR/archive/2.7.0c.tar.gz
    tar -xzf 2.7.0c.tar.gz
    cd STAR-2.7.0c
    
    # Alternatively, get STAR source using git
    git clone https://github.com/alexdobin/STAR.git

Compile under Linux

    # Compile
    cd STAR/source
    make STAR

### Samtools  

Download the samtools package at: http://www.htslib.org/download/  

Configure samtools for linux:

    cd samtools-1.x    # and similarly for bcftools and htslib
    ./configure --prefix=/where/to/install
    make
    make install

For more information, please see the [samtools manual](http://www.htslib.org/doc/samtools.html "tittle")

### BEDtools

Installation of BEDtools for linux:

    wget https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz
    tar -zxvf bedtools-2.25.0.tar.gz
    cd bedtools2
    make

For more information, please see the [BEDtools tutorial](http://quinlanlab.org/tutorials/bedtools/bedtools.html "tittle")

## install R libraries  

---

Open the R console:

    install.packages("WriteXLS")

    install.packages("optparse")

    install.packages("Cairo")

**NB:** SpliceLauncher pipeline requires also the perl compiler but not particular perl libraries

## Download SpliceLauncher pipeline  

---

    git clone https://github.com/raphaelleman/SpliceLauncher
    cd ./SpliceLauncher

## Make the reference files  

---

### STAR genome
To create STAR genome you will need :

* genome fasta file
* RefSeq annot GTF file
* RefSeq annot BED file

#### get genome fasta file

From [UCSC](http://hgdownload.soe.ucsc.edu/downloads.html#human "tittle"), with hg19 example:

    #the ftp URL depends on your assembly genome choice
    mkdir ./fastaGenome
    cd ./fastaGenome
    `wget --timestamping ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz -O chromFa.tar.gz`
    tar xvzf ./chromFa.tar.gz
    cd ..

#### get RefSeq annot GTF file

From [UCSC table browser](https://genome.ucsc.edu/cgi-bin/hgTables "tittle"), with hg19 example:

Options setting:

+ *clade:* Mammal
+ *genome:* Human
+ *assembly:* Feb. 2009 (GRCh37/hg19) #you can choice here your genome assembly
+ *group:* Genes and Gene Predictions
+ *track:* NCBI RefSeq
+ *table:* RefSeq All (ncbiRefSeq)
+ *output format:* GTF - gene transfer format (limited)
+ *output file:* RefSeqAnnot.gtf.gz

After download the file, uncompress it by `gunzip RefSeqAnnot.gtf.gz`

#### get RefSeq annot BED file

From [UCSC table browser](https://genome.ucsc.edu/cgi-bin/hgTables "tittle"), with hg19 example:

Options setting:

+ *clade:* Mammal
+ *genome:* Human
+ *assembly:* Feb. 2009 (GRCh37/hg19) #you can choice here your genome assembly
+ *group:* Genes and Gene Predictions
+ *track:* NCBI RefSeq
+ *table:* RefSeq All (ncbiRefSeq)
+ *output format:* BED - browser extensible data
+ *output file:* RefSeqAnnot.bed.gz
+ *Create one BED record per:* Whole gene

After download the file, uncompress it by `gunzip RefSeqAnnot.bed.gz`

#### create STAR genome

**First step, we create the intron coordinates file from BED RefSeq file**
the generated file is an sjdb file that will use by STAR to get the junction reads

    cd /path/to/SpliceLauncher/
    Rscript ./scripts/generateRefSeqsjdb.r -i /path/to/RefSeqAnnot.bed -o ./RefSeqAnnot.sjdb

**Second step, we create the genome STAR files**

    mkdir ./genomeSTAR
    /path/to/STAR \ 
     --runMode genomeGenerate \
     --runThreadN 5 #define here the number of thread to use \
     --genomeDir ./genomeSTAR \
     --genomeFastaFiles ./fastaGenome/*.fa \
     --sjdbFileChrStartEnd ./RefSeqAnnot.sjdb \
     --sjdbGTFfile /path/to/RefSeqAnnot.gtf \
     --sjdbOverhang 99

### create exon BED annotations

This file is provide in this repository at */refData/refExons.bed*  

At this step you can also define a list of transcripts to get only the splicing junctions on these transcripts (especially for targeted RNAseq). An list example is provide in *./dataTest/transcriptsToSelect.txt*.
If you wish to generate new exon BED annotations, you can use the RefSeqAnnot.bed file generated in the **STAR genome** section.  
the command is:

    Rscript ./scripts/generateExonBEDRef.r -i /path/to/RefSeqAnnot.bed \
     -o ./RefExons.bed \
     -t /path/to/transcriptsList.txt #optional, the list of selected trasncript (example of list in ./dataTest/transcriptsToSelect.txt)

### create the transcripts information

This file is provide in this repository at */refData/RefSpliceLauncher.txt*  

If you wish to generate new transcripts information file, **at first step** you need to download the RefSeq annotation database from [UCSC table browser](https://genome.ucsc.edu/cgi-bin/hgTables "tittle"), with hg19 example:


Options setting:

+ *clade:* Mammal
+ *genome:* Human
+ *assembly:* Feb. 2009 (GRCh37/hg19) #you can choice here your genome assembly
+ *group:* Genes and Gene Predictions
+ *track:* NCBI RefSeq
+ *table:* RefSeq All (ncbiRefSeq)
+ *output format:* all fields from selected table
+ *output file:* RefSeqAnnot.txt.gz

After download the file, uncompress it by `gunzip RefSeqAnnot.txt.gz`

**In the second step**, convert this file in the transcript information file by the command:

    Rscript ./scripts/generateSpliceLauncherRef.r -i /path/to/RefSeqAnnot.txt -o ./RefSpliceLauncher.txt


## Run SpliceLauncher Pipeline

---

The pipeline works in two step, fisrt step is to get the read count matrix from fastq files, the second step is to process to SpliceLauncher analysis

### get the read count from fastq files

This part of the pipeline is in the shell script **_pipelineRNAseq.sh_**  
To see the different options of this script `pipelineRNAseq.sh --help`

With the example data provided in this repository, single end RNAseq (2x75pb) on *BRCA1* and *BRCA2* transcripts:

    cd /path/to/SpliceLauncher
    bash ./pipelineRNAseq.sh \
     -F ./dataTest/fastq/ \
     -O ./testSpliceLauncher/ \
     -g ./genomeSTAR/ \
     --bedannot ./refData/refExons.bed \
     --star /path/to/STAR \
     --samtools /path/to/samtools \
     --bedtools /path/to/BEDtools/bin/

After running, two folders and one file are created in the directory output. The *Bam* folder contains the BAM files with their index and STAR log files. The folder *getClosestExons* contains the BED files and txt files that correspond to the junction coordinates and junction counts respectively. The file is the matrix count could be use by SpliceLauncher tool.

### SpliceLauncher analysis

---

To launch SpliceLauncher analysis, you need the matrix count and the transcript information file. An example of these two files is provide in *./dataTest/MatrixCountExample.txt* and *./refData/RefSpliceLauncher.txt* respectively.

An example of SpliceLauncher command:

    cd /path/to/SpliceLauncher/
    Rscript ./SpliceLauncher -I ./dataTest/MatrixCountExample.txt -R ./refData/RefSpliceLauncher.txt -O ./testSpliceLauncher

 All results (Excel file and pdf illustrations) are save in the folder *testSpliceLauncher_results*.


