# SpliceLauncher

---

The RNAseq pipeline tool to study the alternative splicing. The pipeline works in two step, fisrt step is to get the read count matrix from fastq files, the second step is to process to SpliceLauncher analysis.

![SpliceLauncher](https://github.com/raphaelleman/SpliceLauncher/blob/master/refData/Figure1.png)

**Table**

* [Repository contents](#1)
* [Getting started](#2)
    * [Prerequisites](#3)
        * [STAR](#4)
        * [Samtools](#5)
        * [BEDtools](#6)
        * [Install R libraries](#7)
    * [Installing SpliceLauncher](#8)
        * [Creating the reference files](#9)
* [Running the SpliceLauncher tests](#10)
    * [RNAseq pipeline, get the read count from fastq files](#11)
    * [SpliceLauncher analysis](#12)
* [RNAseq pipeline Options](#13)
* [SpliceLauncher Options](#14)
* [Authors](#15)
* [License](#16)

## Repository contents<a id="1"></a>

---

* dataTest: example of input files
* refData: reference files that SpliceLauncher needs, here in hg19 assambly
* scripts: complementary scripts to run SpliceLauncher pipeline

## Getting started<a id="2"></a>

---

### Prerequisites<a id="3"></a>
The SpliceLauncher pipeline needs to start from fastq files
The tools:

* STAR (v2.6 or later)
* samtools (v1.3 or later)
* BEDtools (v2.17 or later)
* R librairies *WriteXLS* and *Cairo*

**NB:** SpliceLauncher pipeline requires also the perl compiler but not particular perl libraries

#### STAR <a id="4"></a>

these following instruction were from the [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf "Title")

Get the g++ compiler for linux

```Bash
sudo apt-get update
sudo apt-get install g++
sudo apt-get install make
```

Download the [latest release](https://github.com/alexdobin/STAR/releases "Title") from and uncompress it

```Bash
# Get latest STAR source from releases
wget https://github.com/alexdobin/STAR/archive/2.7.0c.tar.gz
tar -xzf 2.7.0c.tar.gz
cd STAR-2.7.0c

# Alternatively, get STAR source using git
git clone https://github.com/alexdobin/STAR.git
```

Compile under Linux

```Bash
# Compile
cd STAR/source
make STAR
```

#### Samtools <a id="5"></a>

Download the samtools package at: http://www.htslib.org/download/

Configure samtools for linux:

```Bash
cd samtools-1.x    # and similarly for bcftools and htslib
./configure --prefix=/where/to/install
make
make install
```

For more information, please see the [samtools manual](http://www.htslib.org/doc/samtools.html "tittle")

#### BEDtools <a id="6"></a>

Installation of BEDtools for linux:

```Bash
wget https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz
tar -zxvf bedtools-2.25.0.tar.gz
cd bedtools2
make
```

For more information, please see the [BEDtools tutorial](http://quinlanlab.org/tutorials/bedtools/bedtools.html "tittle")

#### Install R libraries <a id="7"></a>

Open the R console:

```R
install.packages("WriteXLS")
install.packages("Cairo")
```

## Installing SpliceLauncher <a id="8"></a>

---

```Bash
git clone https://github.com/raphaelleman/SpliceLauncher
cd ./SpliceLauncher
```

### Creating the reference files <a id="9"></a>

SpliceLauncher needs the reference genome and its annotation file in GFF3 format.

1. Donwload genome fasta file
2. Donwload annot GFF file
3. Extract databases used by SpliceLauncher
4. Compile STAR genome

Steps:
1. Download Fasta genome: from [UCSC](http://hgdownload.soe.ucsc.edu/downloads.html#human "tittle"), with hg19 example:
    ```Bash
    #the ftp URL depends on your assembly genome choice
    cd /path/to/SpliceLauncher/
    mkdir ./fastaGenome
    cd ./fastaGenome
    wget --timestamping ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz -O chromFa.tar.gz
    tar xvzf ./chromFa.tar.gz
    cd ..
    ```

2. Donwload the GFF annotation file, either from [RefSeq ftp server](ftp://ftp.ncbi.nlm.nih.gov/refseq/ "tittle") or from [Gencode](https://www.gencodegenes.org/ "tittle").
For human hg19 annotation file from RefSeq: ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.gff.gz
```Bash
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.gff.gz
gunzip ./GRCh37_latest_genomic.gff.gz
head ./GRCh37_latest_genomic.gff
##gff-version 3
#!gff-spec-version 1.21
#!processor NCBI annotwriter
#!genome-build GRCh37.p13
#!genome-build-accession NCBI_Assembly:GCF_000001405.25
#!annotation-date
#!annotation-source
##sequence-region NC_000001.10 1 249250621
##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=9606
NC_000001.10	RefSeq	region	1	249250621	.	+	.	ID=id0;Dbxref=taxon:9606;Name=1;chromosome=1;gbkey=Src;genome=chromosome;mol_type=genomic DNA
NC_000001.10	BestRefSeq	gene	11874	14409	.	+	.	ID=gene0;Dbxref=GeneID:100287102,HGNC:HGNC:37102;Name=DDX11L1;description=DEAD/H-box helicase 11 like 1;gbkey=Gene;gene=DDX11L1;gene_biotype=misc_RNA;pseudo=true
NC_000001.10	BestRefSeq	transcript	11874	14409	.	+	.	ID=rna0;Parent=gene0;Dbxref=GeneID:100287102,Genbank:NR_046018.2,HGNC:HGNC:37102;Name=NR_046018.2;gbkey=misc_RNA;gene=DDX11L1;product=DEAD/H-box helicase 11 like 1;transcript_id=NR_046018.2
NC_000001.10	BestRefSeq	exon	11874	12227	.	+	.	ID=id1;Parent=rna0;Dbxref=GeneID:100287102,Genbank:NR_046018.2,HGNC:HGNC:37102;gbkey=misc_RNA;gene=DDX11L1;product=DEAD/H-box helicase 11 like 1;transcript_id=NR_046018.2
NC_000001.10	BestRefSeq	exon	12613	12721	.	+	.	ID=id2;Parent=rna0;Dbxref=GeneID:100287102,Genbank:NR_046018.2,HGNC:HGNC:37102;gbkey=misc_RNA;gene=DDX11L1;product=DEAD/H-box helicase 11 like 1;transcript_id=NR_046018.2
NC_000001.10	BestRefSeq	exon	13221	14409	.	+	.	ID=id3;Parent=rna0;Dbxref=GeneID:100287102,Genbank:NR_046018.2,HGNC:HGNC:37102;gbkey=misc_RNA;gene=DDX11L1;product=DEAD/H-box helicase 11 like 1;transcript_id=NR_046018.2
```

3. Extract databases used by SpliceLauncher

    ```Bash
    cd /path/to/SpliceLauncher/
    Rscript ./scripts/generateSpliceLauncherDB.r -i /path/to/GRCh37_latest_genomic.gff -o ./RefSeqHg19
    ```

    An example is provide in this repository at [refExons.bed](https://github.com/raphaelleman/SpliceLauncher/tree/master/refData/refExons.bed "tittle")

4. Compile STAR genome

    ```Bash
    mkdir ./genomeSTAR
    /path/to/STAR \
     --runMode genomeGenerate \
     --runThreadN 5 #define here the number of thread to use \
     --genomeDir ./genomeSTAR/ \
     --genomeFastaFiles ./fastaGenome/*.fa \
     --sjdbFileChrStartEnd ./RefSeqHg19/RefSeqAnnot.sjdb \
     --sjdbGTFfile /path/to/RefSeqAnnot.gtf \
     --sjdbOverhang 99
    ```

## Running the SpliceLauncher tests<a id="10"></a>

---

The example files are provided in [dataTest](https://github.com/raphaelleman/SpliceLauncher/tree/master/dataTest "tittle")

### RNAseq pipeline, get the read count from fastq files <a id="11"></a>

This part uses the shell script **_pipelineRNAseq.sh_**
To see the different options of this script `pipelineRNAseq.sh --help`

With the example data provided in this repository, single end RNAseq (2x75pb) on *BRCA1* and *BRCA2* transcripts:

```Bash
cd /path/to/SpliceLauncher
bash ./pipelineRNAseq.sh \
 -F ./dataTest/fastq/ \
 -O ./testSpliceLauncher/ \
 -g ./genomeSTAR/ \
 --bedannot ./refData/refExons.bed \
 --star /path/to/STAR \
 --samtools /path/to/samtools \
 --bedtools /path/to/BEDtools/bin/
 ```

After running, two folders and one file are created in the directory output. The *Bam* folder contains the BAM files with their index and STAR log files. The folder *getClosestExons* contains the BED files and txt files that correspond to the junction coordinates and junction counts respectively. The file is the matrix count could be use by SpliceLauncher tool.

### SpliceLauncher analysis <a id="12"></a>

To launch SpliceLauncher analysis, you need the matrix count and the transcript information file. An example of these two files is provide in [MatrixCountExample.txt](https://github.com/raphaelleman/SpliceLauncher/tree/master/dataTest/MatrixCountExample.txt "tittle") and [RefSpliceLauncher.txt](https://github.com/raphaelleman/SpliceLauncher/tree/master/refData/RefSpliceLauncher.txt "tittle") respectively.

SpliceLauncher command:

```shell
cd /path/to/SpliceLauncher/
Rscript ./SpliceLauncher.r -I ./dataTest/MatrixCountExample.txt -R ./refData/RefSpliceLauncher.txt -O ./
```

 The results are saved in the folder *MatrixCountExample_results*. The scheme of this report is:

| Column names | Example | Description |
|------------:|:--------:|:------------:|
| Conca | chr13_32915333_32920963 | The junction id (chr_start_end) |
| chr | chr13 | Chromosome number |
| start | 32915333 | Genomic coordinate of start junction<br/>End if on reverse strand |
| end | 32920963 | Genomic coordinate of end junction<br/>Start if on reverse strand |
| strand | + | Strand of the junction ('+': forward;<br/> '-':reverse) |
| brin | forward | Strand of transcript |
| NM | NM_000059 | The transcript id according RefSeq nomenclature |
| Gene | BRCA2 | Gene symbol |
| *Sample* | 2250 | Read count |
| *P_Sample* | 15.25659623 | % of relative expression |
| constitutive | NoPhysio | If junction supports the reference transcripts, Physio else NoPhysio |
| calcul | SkipEx | The nature of junction:<br/>Physio: Natural junction<br/>SkipEx: Exon skipping<br/>5AS: Donor splice site shift<br/>3AS: Acceptor splice site shift<br/>NoData: Unannotated juntion |
| AnnotJuncs | ∆12 | The junction names |
| cStart | c.6841 | Transcriptomic start coordinate of the junction |
| cEnd | c.6938 | Transcriptomic end coordinate of the junction |
| mean_pourcentage | 12.60242 | Average in % of relative expression across samples |
| read_mean | 2683.769231 | Average of read count across samples |
| nbSamp | 11 | Number of time that the junction has been seen in samples |
| DistribAjust | - | The Distribution of junction expression (Gamma/N.binom) |
| Significative | NO | If a sample shown an abnormal expression of the junction |


## RNAseq pipeline Options <a id="13"></a>

---

**-F, --fastq** /path/to/Fastq files

+ The directory of FastQ files in *.fastq.gz format. If paired end analysis please names the files: XXXXXXX1.fastq.gz for read 1 and XXXXXXX2.fastq.gz for read 2.

**-O, --output** /path/to/output

+ Directory of output files. If inexisting it creates the Directory.

**--star** /path/to/STAR

+ Path to the STAR executable

**-g, --genome** /path/to/STAR genome

+ Directory to the STAR genome (see [STAR genome](#10) section).

**--samtools** /path/to/samtools

+ Directory to the SAMtools executable.

**--bedtools** /path/to/bedtools/bin/

+ Directory of the folder with BEDtools executables (pipeline used bamToBed and closestBed).

**--bedannot** /path/to/BEDannot.bed

+ Exon coordinates of transcripts described in RefSeq database (see [Create exon BED annotations](#11) section).

**-p**

+ Process to a paired-end analysis

**-t, --threads** N

+ Number of CPUs used for the STAR alignment.

**--perlscript** /path/to/perl scripts

+ Path to the Perl script used by the pipeline, by default they are in [scripts](https://github.com/raphaelleman/SpliceLauncher/tree/master/scripts "tittle") folder.

## SpliceLauncher Options <a id="14"></a>

---

**-I, --input** /path/to/matrix of read count

+ The read count matrix file (see an example [MatrixCountExample.txt](https://github.com/raphaelleman/SpliceLauncher/tree/master/dataTest/MatrixCountExample.txt "tittle")) used by SpliceLauncher.

**-O, --output** /path/to/output/

+ The ouput directory. If inexisting SpliceLauncher create the Directory.

**-R, --RefSeqAnnot** /path/to/RefSeannot.txt

+ The RefSeq database withe the transcript details (used by default: [RefSpliceLauncher.txt](https://github.com/raphaelleman/SpliceLauncher/tree/master/refData/RefSpliceLauncher.txt "tittle")). To change this database please refer to [Create the transcripts information](#12) section.

**-S, --SampleNames** sample1|sample2|sample3|...

+ Permit to set the name of samples, used for the pdf files.

**-t, --TranscriptList** /path/to/transcriptList.txt

+ Permits to set the reference transcripts to study the alternative splicing. Does not support several transcripts by gene (see an example [transcriptsToSelect.txt](https://github.com/raphaelleman/SpliceLauncher/blob/master/dataTest/transcriptsToSelect.txt "tittle")). The option *-m, --MergeTranscrit* cann't affect the gene with the selected transcripts.

**--removeOther**

+ Remove genes that the transcripts are not in the list of reference transcripts.

**-b, --BEDannot**

+ Generate a BED file of alternative splicing events, permits to plot the junction on a genome browser such as [UCSC Genome browser](https://genome.ucsc.edu/ "tittle").

**-g, --Graphics**

+ Generate pdf file for each sample, plots the alternative splicing junctions on the reference transcripts. Warning, this option increases significantly the runtime.

**--threshold** N

+ Threshold (in percentage of relative expression) to display the alternative junctions in the pdf RefFiles

**-s, --Statistical**

+ Launch the statistical analysis to detect junctions abnomarly expressed acroos samples, needs at least 5 samples in all to work.

**-a, --Adjust** TRUE

+ Adjust the p-value by Bonferroni method (TRUE/FALSE).

**-n, --NbIntervals** N

+ Number of intervals used in estimation of Negative Binomial distribution

## Authors <a id="15"></a>

---

* Raphael Leman - [raphaelleman](https://github.com/raphaelleman/ "tittle")
    * You can contact me at: r.leman@baclesse.unicancer.fr or raphael.leman@orange.fr

> **Cite as:** SpliceLauncher: a tool for detection, annotation and relative quantification of alternative junctions from RNAseq data.
*Raphaël Leman, Grégoire Davy, Valentin Harter, Antoine Rousselin, Etienne Muller, Alexandre Atkinson, Laurent Castéra, Fréderic Lemoine, Pierre de la Grange, Marine Guillaud-Bataille, Dominique Vaur, Sophie Krieger*

## License <a id="16"></a>

---

This project is licensed under the MIT License - see the [LICENSE](https://github.com/raphaelleman/SpliceLauncher/blob/master/LICENSE "tittle") file for details
