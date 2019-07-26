#!/bin/bash

########## some useful functions
echo_on_stderr () {
    (>&2 echo "$*")
}

test_command_if_exist () {
    to_test=$*
    command -v $to_test >/dev/null 2>&1 || { echo_on_stderr "require $* but it's not installed. Will abort."; return 1; }
}

test_file_if_exist () {
    if [ ! -f $* ]; then
        echo -e "File $* not found! Will abort."
        return 1
    fi
}


workFolder=$(readlink -f $(dirname $0))
scriptPath="${workFolder}/scripts"
########## parsing config file
echo "Parsing configuration file..."
test_file_if_exist ${workFolder}/config.cfg
prev_state=$?
if [ $prev_state -eq 1 ]; then
        exit 1
    fi

source ${workFolder}/config.cfg

threads="1"
endType="single"
workFolder=$(readlink -f $(dirname $0))
ScriptPath=${workFolder}/"scripts"
BEDrefPath=${workFolder}/"refData/refExons.bed"
###runMode:
##INSTALL
##Align
##Count
##SpliceLauncher
install="FALSE"
align="FALSE"
count="FALSE"
spliceLauncher="FALSE"

createDB="FLASE"
createGenome="FLASE"

messageHelp="Usage: $0 [runMode] [options] <command>\n
    \n
    --runMode INSTALL,Align,Count,SpliceLauncher\n
    \tINSTALL\tConfigure the files for SpliceLauncher pipeline\n
    \tAlign\tGenerate BAM files from the FASTQ files\n
    \tCount\tGenerate BED files from the BAM files\n
    \tSpliceLauncher\tGenerate final output from the BED files\n\n
    Option for INSTALL mode\n
    \t-C, --config\t/path/to/configuration file/\t [default: ${workFolder}/config.cfg]\n
    \t-O, --output\t/path/to/output/\tdirectory of the output files\n
    \t--STAR\t/path/to/STAR executable \t[default: ${STAR}]\n
    \t--samtools\t/path/to/samtools executable \t[default: ${samtools}]\n
    \t--bedtools\t/path/to/bedtools/bin folder \t[default: ${bedtools}]\n
    \t--gff\t/path/to/gff file\n
    \t--fasta\t/path/to/fasta genome file\n\n
    Option for Align mode\n
    \t-F, --fastq /path/to/fastq/\n\t\trepository of the FASTQ files\n
    \t-O, --output /path/to/output/\n\t\trepository of the output files\n
    \t-g, --genome /path/to/genome\n\t\tpath to the STAR genome\t[default: ${genome}]\n
    \t--STAR /path/to/STAR\n\t\tpath to the STAR executable\t[default: ${STAR}]\n
    \t--samtools /path/to/samtools\n\t\tpath to samtools executable\t[default: ${samtools}]\n
    \t-p paired-end analysis\n\t\tprocesses to paired-end analysis\t[default: ${endType}]\n
    \t-t, --threads N\n\t\tNb threads used for the alignment\t[default: ${threads}]\n\n
    Option for Count mode\n
    \t-B, --bam /pat/to/BAM files\n
    \t-O, --output /path/to/output/\n\t\trepository of the output files\n
    \t--samtools /path/to/samtools\t\tpath to samtools executable\t[default: ${samtools}]\n
    \t--bedtools /path/to/bedtoolsFolder/\t\tpath to repository of bedtools executables\t[default: ${bedtools}]\n
    \t--bedannot /path/to/bedannot\t\tpath to exon coordinates (in BED format)\t[default: ${bed}]\n\n
    Option for SpliceLauncher mode\n
    \t-I, --input /path/to/inputFile\n\t\tRead count matrix (.txt)\n
    \t-O, --output /path/to/output/\n\t\tDirectory to save the results\n
    \t-R, --RefSeqAnnot /path/to/RefSpliceLauncher.txt\n\t\tRefSeq annotation file name \t[default: ${SL_DB}]\n
    \t--TranscriptList /path/to/transcriptList.txt\n\t\tSet the list of transcripts to use as reference\n
    \t--text\n\t\tPrint main output in txt instead of excel\n
    \t-b, --BEDannot\n\t\tGet the output in BED format\n
    \t--Graphics\n\t\tDisplay graphics of alternative junctions (Warnings: increase the runtime)\n
    \t-n, --NbIntervals 10\n\t\tNb interval of Neg Binom (Integer) [default= 10]\n
    \t--SampleNames name1|name2|name3\n\t\tSample names, '|'-separated, by default use the sample file names\n
    \tIf list of transcripts (--TranscriptList):\n
    \t\t--removeOther\n\t\tRemove the genes with unselected transcripts to improve runtime\n
    \tIf graphics (-g, --Graphics):\n
    \t\t--threshold 1\n\t\tThreshold to shown junctions (%) [default= 1]\n"


if [ $# -lt 1 ]; then
    echo -e $messageHelp
    exit
fi

while [[ $# -gt 0 ]]; do
   key=$1
   case $key in

       --runMode)
       runMode="`readlink -v -f $2`"
       shift 2 # shift past argument and past value
       ;;

       -C|--config)
       conf_file="$2"
       shift 2 # shift past argument and past value
       ;;

       -O|--output)
       out_path="`readlink -v -f $2`"
       shift 2 # shift past argument and past value
       ;;

       --STAR)
       STAR="`readlink -v -f $2`"
       shift 2 # shift past argument and past value
       ;;

       --samtools)
       samtools="`readlink -v -f $2`"
       shift 2 # shift past argument and past value
       ;;

       --bedtools)
       bedtools="`readlink -v -f $2`"
       shift 2 # shift past argument and past value
       ;;

       --gff)
       createDB="TRUE"
       gff_path="`readlink -v -f $2`"
       shift 2 # shift past argument and past value
       ;;

       --fasta)
       createGenome="TRUE"
       fasta_path="`readlink -v -f $2`"
       shift 2 # shift past argument and past value
       ;;

       -F|--fastq)
       fastq_path="`readlink -v -f $2`"
       shift 2 # shift past argument and past value
       ;;

       -g|--genome)
       fasta_genome="`readlink -v -f $2`"
       shift 2 # shift past argument and past value
       ;;

       -p)
       endType="paired"
       shift # shift past argument
       ;;

       -t|--threads)
       threads="$2"
       shift 2 # shift past argument and past value
       ;;

       -B|--bam)
       bam_path="`readlink -v -f $2`"
       shift 2 # shift past argument and past value
       ;;

       --bedannot)
       BEDrefPath="`readlink -v -f $2`"
       shift 2 # shift past argument and past value
       ;;

       -I|--input)
       input_path="`readlink -v -f $2`"
       shift 2 # shift past argument and past value
       ;;

       -R|--RefSeqAnnot)
       SL_DB="`readlink -v -f $2`"
       shift 2 # shift past argument and past value
       ;;

       --TranscriptList)
       TranscriptList="`readlink -v -f $2`"
       shift 2 # shift past argument and past value
       ;;

       --text)
       text="TRUE"
       shift 1 # shift past argument and past value
       ;;

       -b|--BEDannot)
       BEDannot="TRUE"
       shift 1 # shift past argument and past value
       ;;

       --Graphics)
       Graphics="TRUE"
       shift 1 # shift past argument and past value
       ;;

       -n|--NbIntervals)
       NbIntervals="`readlink -v -f $2`"
       shift 2 # shift past argument and past value
       ;;

       --SampleNames)
       SampleNames="`readlink -v -f $2`"
       shift 2 # shift past argument and past value
       ;;

       --removeOther)
       removeOther="TRUE"
       shift 1 # shift past argument and past value
       ;;

       --threshold)
       threshold="`readlink -v -f $2`"
       shift 2 # shift past argument and past value
       ;;

       *)  # unknown option
       POSITIONAL+=("$key") # save it in an array for later
       echo -e "    Unknown option ${key}"
       shift # shift past argument
      ;;
   esac
done

IFS=',' read -a array <<<"${runMode}" # ici, l'IFS n'est modifié que dans l'environnement du read.
for i in ${array[@]}; do
    if [ $i=="INSTALL" ]; then
        install="TRUE"
        echo "$i: ${install}"
    elif [ $i=="Align" ]; then
        align="TRUE"
        echo "$i: ${align}"
    elif [ $i=="Count" ]; then
        count="TRUE"
        echo "$i: ${count}"
    elif [ $i=="SpliceLauncher" ]; then
        spliceLauncher="TRUE"
        echo "$i: ${spliceLauncher}"
    else
        echo "Error in runMode selection!"
        exit
    fi
done

in_error=0 # variable initialization; will be 1 if a file or cmd not exist

# Test if cmd exist
for i in samtools bedtools STAR Rscript perl; do
    test_command_if_exist ${!i}
    prev_state=$?
    if [ $prev_state -eq 1 ]; then
        in_error=1
    else
        echo "${!i} OK."
    fi
done

#switch in INSTALL mode
if [ install=="TRUE" ]; then
    sed -i "s#^samtools=.*#samtools=\"${samtools}\"#" ${workFolder}/config.cfg
    sed -i "s#^bedtools=.*#bedtools=\"${bedtools}\"#" ${workFolder}/config.cfg
    sed -i "s#^STAR=.*#STAR=\"${STAR}\"#" ${workFolder}/config.cfg
fi
exit
# Test if files exist
for i in fasta_genome gff; do
    test_file_if_exist ${!i}
    prev_state=$?
    if [ $prev_state -eq 1 ]; then
        in_error=1
    else
        echo "${i} = ${!i}"
    fi
done

# exit if there is one error or more
if [ $in_error -eq 1 ]; then
    echo -e "=> Aborting."
    exit
fi

echo "Parsing OK."

########## lauch generateSpliceLauncherDB

echo "Will run generateSpliceLauncherDB."

cmd="${Rscript} ${scriptPath}/generateSpliceLauncherDB.r -i ${gff} -o ${out}"
echo -e "$cmd"
$cmd
