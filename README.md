# SpliceLauncher

---

The RNAseq pipeline tool to study the alternative splicing

> **Cite as:** SpliceLauncher: a tool for detection, annotation and relative quantification of alternative junctions from RNAseq data. *Raphaël Leman, Grégoire Davy, Valentin Harter, Antoine Rousselin, Alexandre Atkinson, Laurent Castéra, Fréderic Lemoine, Pierre de la Grange, Dominique Vaur, Sophie Krieger*


**Table**

* [Repository contents](#1)
* [Install and configure SpliceLauncher pipeline](#2)
    * [Install STAR, samtools, BEDtools](#3)
        * [STAR](#4)
        * [Samtools](#5)
        * [BEDtools](#6)
    * [Install R libraries](#7)
    * [Download SpliceLauncher pipeline](#8)
    * [Make the reference files](#9)
        * [STAR genome](#10)
        * [Create exon BED annotations](#11)
        * [Create the transcripts information](#12)
* [Run SpliceLauncher Pipeline](#13)
    * [RNAseq pipeline, get the read count from fastq files](#14)
    * [SpliceLauncher analysis](#15)

## Repository contents<a id="1"></a>

---

* dataTest: example of input files
* refData: reference files that SpliceLauncher needs
* scripts: complementary scripts to run SpliceLauncher pipeline

## Install and configure SpliceLauncher pipeline<a id="2"></a>

---

The SpliceLauncher pipeline needs to start from fastq files
The tools:

* STAR (v2.6 or later)
* samtools (v1.3 or later)
* BEDtools (v2.17 or later)

The reference files:

* Human genome assembly (STAR)
* Exons annotations
* Transcripts information

An example of these two last files is provide in the [refData folder](https://github.com/raphaelleman/SpliceLauncher/tree/master/refData "tittle")

### Install STAR, samtools, BEDtools <a id="3"></a>

---

To use SpliceLauncher pipeline please install the three tools: STAR, samtools, BEDtools

#### STAR <a id="4"></a>

these following instruction were from the [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf "Title")

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

#### Samtools <a id="5"></a>

Download the samtools package at: http://www.htslib.org/download/

Configure samtools for linux:

    cd samtools-1.x    # and similarly for bcftools and htslib
    ./configure --prefix=/where/to/install
    make
    make install

For more information, please see the [samtools manual](http://www.htslib.org/doc/samtools.html "tittle")

#### BEDtools <a id="6"></a>

Installation of BEDtools for linux:

    wget https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz
    tar -zxvf bedtools-2.25.0.tar.gz
    cd bedtools2
    make

For more information, please see the [BEDtools tutorial](http://quinlanlab.org/tutorials/bedtools/bedtools.html "tittle")

### Install R libraries <a id="7"></a>

---

Open the R console:

    install.packages("WriteXLS")

    install.packages("Cairo")

**NB:** SpliceLauncher pipeline requires also the perl compiler but not particular perl libraries

### Download SpliceLauncher pipeline <a id="8"></a>

---

    git clone https://github.com/raphaelleman/SpliceLauncher
    cd ./SpliceLauncher

### Make the reference files <a id="9"></a>

---

#### STAR genome <a id="10"></a>
To create STAR genome you will need :

* genome fasta file
* RefSeq annot GTF file
* RefSeq annot BED file

##### Get genome fasta file

From [UCSC](http://hgdownload.soe.ucsc.edu/downloads.html#human "tittle"), with hg19 example:

    #the ftp URL depends on your assembly genome choice
    cd /path/to/SpliceLauncher/
    mkdir ./fastaGenome
    cd ./fastaGenome
    wget --timestamping ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz -O chromFa.tar.gz
    tar xvzf ./chromFa.tar.gz
    cd ..

##### Get RefSeq annot GTF file

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

##### Get RefSeq annot BED file

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

##### Create STAR genome

**First step, generation of the intron coordinates file from BED RefSeq file**
the generated file is an sjdb file that will use by STAR to get the junction reads

    cd /path/to/SpliceLauncher/
    Rscript ./scripts/generateRefSeqsjdb.r -i /path/to/RefSeqAnnot.bed -o ./RefSeqAnnot.sjdb

**Second step, generation of the genome STAR files**

    mkdir ./genomeSTAR
    /path/to/STAR \
     --runMode genomeGenerate \
     --runThreadN 5 #define here the number of thread to use \
     --genomeDir ./genomeSTAR \
     --genomeFastaFiles ./fastaGenome/*.fa \
     --sjdbFileChrStartEnd ./RefSeqAnnot.sjdb \
     --sjdbGTFfile /path/to/RefSeqAnnot.gtf \
     --sjdbOverhang 99

#### Create exon BED annotations <a id="11"></a>

An example is provide in this repository at [refExons.bed](https://github.com/raphaelleman/SpliceLauncher/tree/master/refData/refExons.bed "tittle")

The exon BED annatations uses the RefSeqAnnot.bed file generated in the **get RefSeq annot BED file** section.
The command is:

    Rscript ./scripts/generateExonBEDRef.r -i /path/to/RefSeqAnnot.bed -o ./refExons.bed
    sort -k1,1 -k2,2n ./refExons.bed > ./refExons.bed

#### Create the transcripts information <a id="12"></a>

An example is provide in this repository at [RefSpliceLauncher.txt](https://github.com/raphaelleman/SpliceLauncher/tree/master/refData/RefSpliceLauncher.txt "tittle")

**In first step**, you need to download the RefSeq annotation database from [UCSC table browser](https://genome.ucsc.edu/cgi-bin/hgTables "tittle"), with hg19 example:


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


## Run SpliceLauncher Pipeline <a id="13"></a>

---

The pipeline works in two step, fisrt step is to get the read count matrix from fastq files, the second step is to process to SpliceLauncher analysis

### RNAseq pipeline, get the read count from fastq files <a id="14"></a>

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

### SpliceLauncher analysis <a id="15"></a>

---

To launch SpliceLauncher analysis, you need the matrix count and the transcript information file. An example of these two files is provide in [MatrixCountExample.txt](https://github.com/raphaelleman/SpliceLauncher/tree/master/dataTest/MatrixCountExample.txt "tittle") and [RefSpliceLauncher.txt](https://github.com/raphaelleman/SpliceLauncher/tree/master/refData/RefSpliceLauncher.txt "tittle") respectively.

An example of SpliceLauncher command:

    cd /path/to/SpliceLauncher/
    Rscript ./SpliceLauncher.r -I ./dataTest/MatrixCountExample.txt -R ./refData/RefSpliceLauncher.txt -O ./

 The results are saved in the folder *MatrixCountExample_results*.

#### SpliceLauncher Options

**-I, --input**
The read count matrix used by SpliceLauncher
**-O, --output**
