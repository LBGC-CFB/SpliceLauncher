#!/bin/bash

## initialize default value
threads="1"
endType=""
in_error=0 # will be 1 if a file or cmd not exist
workFolder=$(readlink -f $(dirname $0))
conf_file="${workFolder}/config.cfg"
scriptPath="${workFolder}/scripts"
BEDrefPath="${workFolder}/refData/refExons.bed"

########## some useful functions
echo_on_stderr () {
    (>&2 echo -e "$*")
}

test_command_if_exist () {
    command -v $* >/dev/null 2>&1 && echo 0 || { echo 1; }
}

# test_command_if_exist () {
    # if [ $(command -v $* >/dev/null 2>&1) -ne 0 ]; then
        # echo -e "not"
    # fi

test_file_if_exist () {
    if [ ! -f $* ]; then
        echo_on_stderr "File $* not found! Will abort."
        echo 1
    else
        echo 0
    fi
}

########## parsing config file
echo "Parsing configuration file..."
for (( i=1; i<=$#; i++)); do # for loop to find config path before reading all other arguments
    next=$((${i}+1))
    # echo "${i} ${!i} ${!next}"
    if [[ ${!i} = '-C' || ${!i} = '--config' ]]; then
        conf_file="${!next}"
        break
    fi
done

if [ $(test_file_if_exist "${conf_file}") -ne 0 ]; then
    in_error=1
else
    source "${conf_file}"
    echo -e "${conf_file} OK.\n"
fi

########## Help message
messageHelp="Usage: $0 [runMode] [options] <command>\n
    \n
    --runMode INSTALL,Align,Count,SpliceLauncher\n
    \tINSTALL\tConfigure the files for SpliceLauncher pipeline\n
    \tAlign\tGenerate BAM files from the FASTQ files\n
    \tCount\tGenerate BED files from the BAM files\n
    \tSpliceLauncher\tGenerate final output from the BED files\n
    \n
    Option for INSTALL mode\n
    \t-C, --config\t/path/to/configuration file/\t [default: ${conf_file}]\n
    \t-O, --output\t/path/to/output/\tdirectory of the output files\n
    \t--STAR\t/path/to/STAR executable \t[default: ${STAR}]\n
    \t--samtools\t/path/to/samtools executable \t[default: ${samtools}]\n
    \t--bedtools\t/path/to/bedtools/bin folder \t[default: ${bedtools}]\n
    \t--gff\t/path/to/gff file\n
    \t--fasta\t/path/to/fasta genome file\n
    \n
    Option for Align mode\n
    \t-F, --fastq /path/to/fastq/\n\t\trepository of the FASTQ files\n
    \t-O, --output /path/to/output/\n\t\trepository of the output files\n
    \t-g, --genome /path/to/genome\n\t\tpath to the STAR genome\t[default: ${genome}]\n
    \t--STAR /path/to/STAR\n\t\tpath to the STAR executable\t[default: ${STAR}]\n
    \t--samtools /path/to/samtools\n\t\tpath to samtools executable\t[default: ${samtools}]\n
    \t-p paired-end analysis\n\t\tprocesses to paired-end analysis\t[default: ${endType}]\n
    \t-t, --threads N\n\t\tNb threads used for the alignment\t[default: ${threads}]\n
    \n
    Option for Count mode\n
    \t-B, --bam /path/to/BAM files\n
    \t-O, --output /path/to/output/\n\t\tdirectory of the output files\n
    \t--samtools\t/path/to/samtools executable \t[default: ${samtools}]\n
    \t--bedtools\t/path/to/bedtools/bin folder \t[default: ${bedtools}]\n
    \t-b, --BEDannot /path/to/your_annotation_file.bed\n\t\tpath to exon coordinates file (in BED format)\t[default: ${bed}]\n
    \n
    Option for SpliceLauncher mode\n
    \t-I, --input /path/to/inputFile\n\t\tRead count matrix (.txt)\n
    \t-O, --output /path/to/output/\n\t\tDirectory to save the results\n
    \t-R, --RefSeqAnnot /path/to/RefSpliceLauncher.txt\n\t\tRefSeq annotation file name \t[default: ${SL_DB}]\n
    \t--TranscriptList /path/to/transcriptList.txt\n\t\tSet the list of transcripts to use as reference\n
    \t--txtOut\n\t\tPrint main output in text instead of xls\n
    \t--bedOut\n\t\tGet the output in BED format\n
    \t--Graphics\n\t\tDisplay graphics of alternative junctions (Warnings: increase the runtime)\n
    \t-n, --NbIntervals 10\n\t\tNb interval of Neg Binom (Integer) [default= 10]\n
    \t--SampleNames name1|name2|name3\n\t\tSample names, '|'-separated, by default use the sample file names\n
    \tIf list of transcripts (--TranscriptList):\n
    \t\t--removeOther\n\t\tRemove the genes with unselected transcripts to improve runtime\n
    \tIf graphics (-g, --Graphics):\n
    \t\t--threshold 1\n\t\tThreshold to shown junctions (%) [default= 1]\n"

## exit if not enough arguments
if [ $# -lt 1 ]; then
    echo -e $messageHelp
    exit
fi

########## runMode:
runMode=""
##INSTALL
install="FALSE"
##Align
align="FALSE"
##Count
count="FALSE"
##SpliceLauncher
spliceLauncher="FALSE"
createDB="FALSE"
createGenome="FALSE"

while [[ $# -gt 0 ]]; do
   key=$1
   case $key in

       --runMode)
       runMode="$2"
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

       -B|--BEDannot)
       BEDrefPath="`readlink -v -f $2`"
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

       -G|--genome)
       fasta_genome="`readlink -v -f $2`"
       shift 2 # shift past argument and past value
       ;;

       -p)
       endType="-p"
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

       -I|--input)
       input_path="`readlink -v -f $2`"
       shift 2 # shift past argument and past value
       ;;

       -R|--RefSeqAnnot)
       SL_DB="`readlink -v -f $2`"
       shift 2 # shift past argument and past value
       ;;

       --transcriptList)
       TranscriptList="`readlink -v -f $2`"
       shift 2 # shift past argument and past value
       ;;

       --text)
       text="TRUE"
       shift 1 # shift past argument and past value
       ;;

       --BEDout)
       BEDout="TRUE"
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

echo "runMode ${runMode}"

IFS=',' read -a array <<<"${runMode}" # here, IFS change is local
for i in ${array[@]}; do
    if [[ ${i} = "INSTALL" ]]; then
        install="TRUE"
        echo "$i: ${install}"
    elif [[ ${i} = "Align" ]]; then
        align="TRUE"
        echo "$i: ${align}"
    elif [[ ${i} = "Count" ]]; then
        count="TRUE"
        echo "$i: ${count}"
    elif [[ ${i} = "SpliceLauncher" ]]; then
        spliceLauncher="TRUE"
        echo "$i: ${spliceLauncher}"
    else
        echo -e "Error in runMode selection! ${i} unknown."
        in_error=1
        # exit
    fi
done

# Test if cmd exist
for i in samtools bedtools STAR Rscript perl; do
    if [[ -z ${!i} || $(test_command_if_exist ${!i}) -ne 0 ]]; then
        in_error=1
        echo_on_stderr "require ${i} but it's not installed. Will abort."
    else
        echo "${!i} OK."
    fi
done

if [[ -z ${out_path} ]]; then
    echo_on_stderr "require Output Path but it's not installed. Will abort."
    in_error=1
fi

echo "Parsing OK."

########## switch in INSTALL mode
if [[ ${install} = "TRUE" ]]; then
    sed -i "s#^samtools=.*#samtools=\"${samtools}\"#" ${conf_file}
    sed -i "s#^bedtools=.*#bedtools=\"${bedtools}\"#" ${conf_file}
    sed -i "s#^STAR=.*#STAR=\"${STAR}\"#" ${conf_file}

## launch generateSpliceLauncherDB

    # Test if files exist
    for i in fasta_path gff_path; do
        if [[ -z ${!i} || $(test_file_if_exist "${!i}") -ne 0 ]]; then
            in_error=1
            echo_on_stderr "${i} not found! Will abort."
        else
            echo "${i} = ${!i}"
        fi
    done

    # exit if there is one error or more
    if [ $in_error -eq 1 ]; then
        echo -e "=> Aborting."
        exit
    fi
    
    # run generateSpliceLauncherDB
    mkdir -p ${out_path}
    echo "Will run generateSpliceLauncherDB."
    cmd="${Rscript} ${scriptPath}/generateSpliceLauncherDB.r -i ${gff_path} -o ${out_path}"
    echo -e "$cmd"
    $cmd
    
    genome="${out_path}/STARgenome"
    sed -i "s#^genome=.*#genome=\"${genome}\"#" ${conf_file}
    mkdir -p ${genome}
    cmd="${STAR} \
    --runMode genomeGenerate \
    --runThreadN ${threads} \
    --genomeDir ${genome} \
    --genomeFastaFiles ${fasta_path} \
    --sjdbFileChrStartEnd ${out_path}/SJDBannotation.sjdb \
    --sjdbGTFfile ${gff_path} \
    --sjdbOverhang 99"
    echo -e "$cmd"
    $cmd
     
fi

########## lauch RNAseq
if [[ ${align} = "TRUE" ]]; then

## launch alignment

    # Test if files exist
    for i in fastq_path genome; do
        if [[ ! -d ${!i} ]]; then
            in_error=1
            echo_on_stderr "${i} not found! Will abort."
        else
            echo "${i} = ${!i}"
        fi
    done

    # exit if there is one error or more
    if [ $in_error -eq 1 ]; then
        echo -e "=> Aborting."
        exit
    fi
    
    # run alignment
    mkdir -p ${out_path}
    echo "Will run alignment."
    cmd="${scriptPath}/pipelineRNAseq.sh -F ${fastq_path} -O ${out_path} -g ${genome} --STAR ${STAR} --samtools ${samtools} -t ${threads} ${endType}"
    echo -e "$cmd"
    $cmd
fi
    
    
    
    
    
    
    
    
    
    
    