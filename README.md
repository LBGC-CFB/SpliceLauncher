# SpliceLauncher

---

SpliceLauncher is a pipeline tool to study the alternative splicing. It works in three steps:
* Get a read count matrix from fastq files, by a dedicated RNAseq pipeline (A step in diagram below).
* Generate data files used hereafter (B step in diagram below)
* Run SpliceLauncher from a read count matrix (C step and furthermore in diagram below).

![SpliceLauncher](https://github.com/LBGC-CFB/SpliceLauncher/blob/master/doc/Figure1.png)

**Table of contents**

* [Repository contents](#1)
* [Prerequisites to install SpliceLauncher](#2)
    * [STAR](#3)
    * [Samtools](#4)
    * [BEDtools](#5)
    * [Install R libraries](#6)
* [Installing SpliceLauncher](#7)
    * [Singularity image](#8)
    * [Download the reference files](#9)
    * [Configure SpliceLauncher with INSTALL mode](#10)
* [Running the SpliceLauncher tests](#11)
    * [Output directory tree](#12)
* [SpliceLauncher options](#13)
    * [Option for INSTALL mode](#14)
    * [Option for Align mode](#15)
    * [Option for Count mode](#16)
    * [Option for SpliceLauncher mode](#17)
* [Authors](#18)
* [License](#19)

---

## Repository contents<a id="1"></a>


* dataTest: example of input files
* scripts: complementary scripts to run SpliceLauncher

## Prerequisites to install SpliceLauncher<a id="2"></a>


The SpliceLauncher pipeline needs to install the following tools and R librairies:

* STAR (v2.7 or later)
* samtools (v1.3 or later)
* BEDtools (v2.17 or later)
* R with *WriteXLS* and *Cairo* packages
* Perl

### STAR <a id="3"></a>

Following instruction were from the [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf "Title")

Get the g++ compiler for linux

```Bash
sudo apt-get update
sudo apt-get install g++
sudo apt-get install make
```

Download the [latest release](https://github.com/alexdobin/STAR/releases/latest "Title") and uncompress it

```Bash
# Get latest STAR source
version="2.7.0c"
wget https://github.com/alexdobin/STAR/archive/${version}.tar.gz
tar -xzf ${version}.tar.gz
cd STAR-${version}

# Alternatively, get STAR source using git
git clone https://github.com/alexdobin/STAR.git
```

Compile under Linux

```Bash
# Compile
cd STAR/source
make STAR
```

### Samtools <a id="4"></a>

Download the samtools package at: https://github.com/samtools/samtools/releases/latest

Configure samtools for linux:

```Bash
cd samtools-1.x
./configure --prefix=/where/to/install
make
make install
```

For more information, please see the [samtools manual](http://www.htslib.org/doc/samtools.html "Title")

### BEDtools <a id="5"></a>

Installation of BEDtools for linux:

```Bash
wget https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz
tar -zxvf bedtools-2.25.0.tar.gz
cd bedtools2
make
```

For more information, please see the [BEDtools tutorial](http://quinlanlab.org/tutorials/bedtools/bedtools.html "Title")

### Install R libraries <a id="6"></a>

The library WriteXLS allows to save result in xlsx format if you do not want to install it, use the *--txtOut* option. 
The library Cairo allows to print result in pdf format if you do not want to install it, do not add *--Graphics* option.
Open the R console:

```R
install.packages("WriteXLS")
install.packages("Cairo")
```

## Installing SpliceLauncher <a id="7"></a>


Download the latest release from of SpliceLauncher source using git

```Bash
git clone https://github.com/LBGC-CFB/SpliceLauncher
cd ./SpliceLauncher
```

### Singularity image <a id="8"></a>

As the Singularity image config file is not writable, you need to use a local version of the  [config file](https://github.com/LBGC-CFB/SpliceLauncher/blob/master/config.cfg "tittle")
SpliceLauncher and all its dependencies are also integrated in a Singularity image:

1. To build it:

```Bash
sudo singularity build /path/to/SpliceLauncher.simg /path/to/splicelauncher.recipe
```

2. To use it

```Bash
sudo singularity run /path/to/SpliceLauncher.simg --config /path/to/my_config.cfg --help
```

### Download the reference files <a id="9"></a>


The reference files are the genome (Fasta) and the corresponding annotation file (GFF3):

1. Reference genome in fasta format
2. The annotation file in GFF v3 format

Steps:
1. Download Fasta genome: from [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/) FTP server or from [Gencode](https://www.gencodegenes.org/ "tittle").

For example, human hg19 genome file from RefSeq:
```Bash
#the ftp URL depends on your assembly genome choice
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.fna.gz
gunzip ./GRCh37_latest_genomic.fna.gz
```

2. Download the GFF annotation file, either from [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/ "tittle") FTP server or from [Gencode](https://www.gencodegenes.org/ "tittle").

For example, human hg19 annotation file from RefSeq:
```Bash
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.gff.gz
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
```

3. [Optional] Convert contig of RefSeq to UCSC chromosome names, you will need to download the assembly report, an [example](https://github.com/LBGC-CFB/SpliceLauncher/blob/master/dataTest/Example_assembly_report.txt "tittle") of this report is provide in dataTest folder but it is a truncating example so do not use for your own genome. For url example of an assembly report of GRCh37 opf RefSeq can be dowload from https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_assembly_report.txt 
Use the **--assembly_report** of INSTALL mode to launch the convertion. This step usually takes one hour.

4. [Optionnal] To reduce needed memory, we can also restrict the analysis to the primary assembly, without unplaced contigs:
```
grep ">" GRCh37_latest_genomic.fna | grep -v "unplaced genomic contig"| grep -v "unlocalized genomic contig" | grep -v "genomic patch"| grep -v "alternate locus" | sed 's/^>//' > chr_names
seqtk subseq GRCh37_latest_genomic.fna chr_names > GRCh37_latest_genomic.sub.fna
cut -f 1 -d ' ' chr_names > chr_names_id
head -n 9 GRCh37_latest_genomic.gff > GRCh37_latest_genomic.sub.gff
grep -f chr_names_id GRCh37_latest_genomic.gff >> GRCh37_latest_genomic.sub.gff
```

### Configure SpliceLauncher with INSTALL mode <a id="10"></a>

SpliceLauncher comes with a ready to use config.cfg file. It contains the paths of software and files used by SpliceLauncher. The INSTALL mode of SpliceLauncher updates this config cfg file. If you define the path to GFF (v3) file and path to the FASTA genome, the INSTALL mode will extract all necessary information from this GFF and indexing the STAR genome. This informations are stored in a BED file that contains the exon coordinates, in a sjdb file that contains the intron coordinates and a text file that contains the details of transcript structures. You need to define where these files will saving by the `-O, --output` argument

Use INSTALL mode of SpliceLauncher:

```Bash
cd /path/to/SpliceLauncher/
mkdir ./refSpliceLauncher # Here this folder will contain the reference files used by SpliceLauncher
bash ./SpliceLauncher.sh --runMode INSTALL \
    -O ./refSpliceLauncher \
    --STAR /path/to/STAR \
    --samtools /path/to/samtools \
    --bedtools /path/to/bedtools \
    --gff /path/to/gff \
    --threads < number of thread > \
    --fasta /path/to/fasta
```

## Running the SpliceLauncher tests<a id="11"></a>


The example files are provided in [dataTest](https://github.com/LBGC-CFB/SpliceLauncher/tree/master/dataTest "tittle"), with the example data provided in single end RNAseq (1x75pb) on *BRCA1* and *BRCA2* transcripts:

```Bash
cd /path/to/SpliceLauncher
bash ./SpliceLauncher.sh --runMode Align,Count,SpliceLauncher -F ./dataTest/fastq/ -O ./testSpliceLauncher/ \
    # Optional \
    -t <number of thread> \
    -m <allowed memory in bits> \
    --Graphics \
    --tmpDir /path/to/tmpDir # path to save tmp file during alignment
 ```

### Output directory tree <a id="12"></a>
```
SpliceLauncher/outdir
├── Bam
│   ├── {sample}.Aligned.sortedByCoord.out.bam
|   ├── {sample}.Aligned.sortedByCoord.out.bam.csi
|   ├── {sample}.Aligned.sortedByCoord.out_juncs.bed
|   ├── {sample}.SJ.out.tab
│   └── ...
├── getClosestExons
│   ├── {sample}.Aligned.sortedByCoord.out.count
│   └── ...
├── {run name}_results
│   ├── {run name}_figures_output
│       ├── {run name}_{sample}.pdf
│       └── ...
│   ├── {run name}_outputSpliceLauncher.xlsx
│   └── {run name}.bed
├── {run name}.txt
├── {run name}_report_{run date}.txt
```

The results of SpliceLauncher analysis are in **{run name}_results**.

The final results are displayed in the file *{run name}_outputSpliceLauncher.xlsx*. The scheme of this file is:

| Column names | Example | Description |
|------------:|:--------:|:------------:|
| Conca | chr13_32915333_32920963 | The junction id (chr_start_end) |
| chr | chr13 | Chromosome number |
| start | 32915333 | Genomic coordinate of start junction<br/>End if on reverse strand |
| end | 32920963 | Genomic coordinate of end junction<br/>Start if on reverse strand |
| strand | + | Strand of the junction ('+': forward;<br/> '-':reverse) |
| Strand_transcript | forward | Strand of transcript |
| NM | NM_000059 | The transcript id according RefSeq nomenclature |
| Gene | BRCA2 | Gene symbol |
| *Sample* | 2250 | Read count |
| *P_Sample* | 15.25659623 | % of relative expression |
| event_type | SkipEx | The nature of junction:<br/>Physio: Natural junction<br/>SkipEx: Exon skipping<br/>5AS: Donor splice site shift<br/>3AS: Acceptor splice site shift<br/>NoData: Unannotated juntion |
| AnnotJuncs | ∆12 | The junction names |
| cStart | c.6841 | Transcriptomic start coordinate of the junction |
| cEnd | c.6938 | Transcriptomic end coordinate of the junction |
| mean_percent | 12.60242 | Average in % of relative expression across samples |
| read_mean | 2683.769231 | Average of read count across samples |
| nbSamp | 11 | Number of time that the junction has been seen in samples |
| DistribAjust | - | The Distribution of junction expression (Gamma/N.binom) |
| Significative | NO | If a sample shown an abnormal expression of the junction |
| filterInterpretation | Unique junction | If a sample had an abnormal expression: Aberrant junction<br/>For unmodelized junction, if max expression >1%: Unique junction  |

## SpliceLauncher options <a id="13"></a>


**--runMode** INSTALL,Align,Count,SpliceLauncher
* The runMode defines the steps of analysis with:
    * INSTALL: Updates the config.cfg file for SpliceLauncher pipeline
    * Align: Generates BAM files from the FASTQ files
    * Count: Generates the matrix read count from the BAM files
    * SpliceLauncher: Generates final output from the matrix read count

### Option for INSTALL mode <a id="14"></a>

**-C, --config** /path/to/configuration file/
* Path to the config.cfg file, **only** if you want to use your own config file

**-O, --output** /path/to/output/
* Directory to save the reference files (BED, sjdb, txt) and the indexed genome

**--STAR** /path/to/STAR
* Path to the STAR executable

**--samtools** /path/to/samtools
* Path to the samtools executable

**--bedtools** /path/to/bedtools
* Path to the bedtools executable

**--gff** /path/to/gff file
* Path to the GFF file (v3)

**--fasta** /path/to/fasta
* Path to the genome fasta file

**--mane** /path/to/MANElistFile.txt
* List of MANE transcripts, current version donwloaded from https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.1/MANE.GRCh38.v1.1.summary.txt.gz

**--assembly_report** /path/to/GRChXX_latest_assembly_report.txt
* Path to assembly report to convert contig in UCSC chromosome names

**-t, --threads** N
* Nb threads used to index the STAR genome

### Option for Align mode <a id="15"></a>

**-F, --fastq** /path/to/fastq/
* Repository of the FASTQ files

**-O, --output** /path/to/output/
* Repository of the output files

**-p**
* Processes to paired-end analysis

**-t, --threads** N
* Nb threads used for the alignment

**--tmpDir** /path/to/tmp directory of files from alignment
* path to directory where tempory files will saved during alignment, please point toward empty directory

**-g, --genome** /path/to/genome
* Path to the genome directory, **only** if you to use a genome directory different of the genome defined in config.cfg file

**--STAR** /path/to/STAR
* Path to the STAR executable, **only** if you to use a STAR software different of the STAR defined in config.cfg file

**--samtools** /path/to/samtools
* Path to the samtools executable, **only** if you to use a samtools software different of the samtools defined in config.cfg file

### Option for Count mode <a id="16"></a>

**-B, --bam** /path/to/BAM files
* Repository of the BAM folder

**-O, --output** /path/to/output/
* Repository of the output files

**--bedtools**\t/path/to/bedtools
* Path to the bedtools executable, **only** if you to use a bedtools software different of the bedtools defined in config.cfg file

**-b, --BEDannot** /path/to/your_annotation_file.bed
* Path to exon coordinates file (in BED format), **only** if you to use exon coordinates different of the coordinates defined in config.cfg file

**--fastCount**
* improve junction count runtime, in addtition increase the specificty but reduce the sensitivity

**-p**
* Processes to paired-end analysis, deprecated if use `--fastCount` option

**--samtools** /path/to/samtools
* Path to the samtools executable, **only** if you to use a samtools software different of the samtools defined in config.cfg file, deprecated if use `--fastCount` option

### Option for SpliceLauncher mode <a id="17"></a>

**-I, --input**/path/to/inputFile
* Read count matrix (.txt)

**-O, --output** /path/to/output/
* Directory to save the results

**--transcriptList** /path/to/transcriptList.txt
* Set the list of transcripts to use as reference

**--txtOut**
* Print main output in text instead of xls

**--bedOut**
* Get the output in BED format

**--Graphics**
* Display graphics of alternative junctions (Warnings: increase the runtime)

**-n, --NbIntervals** 10
* Nb interval of Neg Binom (Integer)

**--min_cov** 5
* Minimal number of read supporting a junction (Integer)

**--SampleNames** name1|name2|name3
* Sample names, '|'-separated, by default use the sample file names

If list of transcripts (--TranscriptList):
**--removeOther**
* Remove the genes with unselected transcripts to improve runtime

If graphics (-g, --Graphics):
**--threshold** 1
* Threshold to shown junctions (%)

**-R, --RefSeqAnnot** /path/to/RefSpliceLauncher.txt
* Transcript information file, **only** if you to use a transcript information file different of file defined in config.cfg file


## Authors <a id="18"></a>


* Raphael Leman - [raphaelleman](https://github.com/raphaelleman/ "tittle")
    * You can contact me at: r.leman@baclesse.unicancer.fr or raphael.leman@orange.fr

> **Cite as:** SpliceLauncher: a tool for detection, annotation and relative quantification of alternative junctions from RNAseq data
Raphaël Leman, Valentin Harter, Alexandre Atkinson, Grégoire Davy, Antoine Rousselin, Etienne Muller, Laurent Castéra, Fréderic Lemoine, Pierre de la Grange, Marine Guillaud-Bataille, Dominique Vaur, Sophie Krieger, [Bioinformatics](https://doi.org/10.1093/bioinformatics/btz784 "tittle")

## License <a id="19"></a>


This project is licensed under the MIT License - see the [LICENSE](https://github.com/LBGC-CFB/SpliceLauncher/blob/master/LICENSE "tittle") file for details
