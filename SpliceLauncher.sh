#!/bin/bash
set -e

#########################
#RNAseq pipeline
#########################

#author Raphael Leman r.leman@baclesse.unicancer.fr, Center François Baclesse and Normandie University, Unicaen, Inserm U1245
#Copyright 2019 Center François Baclesse and Normandie University, Unicaen, Inserm U1245

#This software was developed from the work:
#SpliceLauncher: a tool for detection, annotation and relative quantification of alternative junctions from RNAseq data.
#Raphaël Leman, Grégoire Davy, Valentin Harter, Antoine Rousselin, Alexandre Atkinson, Laurent Castéra, Dominique Vaur,
#Fréderic Lemoine, Pierre De La Grange, Sophie Krieger

#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"),
#to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute,
#sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
#MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
#FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
#WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## initialize default value
threads="1"
memory=`free -t | grep "Mem:" | awk ' { print $2 } ' | awk '{printf "%.f", $1*0.5}'`
memory=`echo "${memory}*1000" | bc`
endType=""
min_cov="5"
in_error=0 # will be 1 if a file or cmd not exist
workFolder=$(readlink -f $(dirname $0))
conf_file="${workFolder}/config.cfg"
scriptPath="${workFolder}/scripts"
manePath="${scriptPath}/MANE.to.splicelauncher.txt"
BEDrefPath="${workFolder}/refData/refExons.bed"
fastCount=false
removeOther=""
tmpAlignDir=""
text=""
Graphics=""
NbIntervals=10
threshold=1
pathToConvert="-"

########## some useful functions
echo_on_stderr () {
    (>&2 echo -e "$*")
}

test_command_if_exist () {
    command -v $* >/dev/null 2>&1 && echo 0 || { echo 1; }
}

test_file_if_exist () {
    if [ ! -f $* ]; then
        echo_on_stderr "File $* not found! Will abort."
        echo 1
    else
        echo 0
    fi
}

########## parsing config file
echo -e "###############################################"
echo -e "####### Check SpliceLauncher environment"
echo -e "###############################################\n"

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
    \t-O, --output\t/path/to/output/\n\t\tDirectory of the output files\n
    \t--STAR\t/path/to/STAR executable \t[default: ${STAR}]\n
    \t--samtools\t/path/to/samtools executable \t[default: ${samtools}]\n
    \t--bedtools\t/path/to/bedtools/bin folder \t[default: ${bedtools}]\n
    \t--gff\t/path/to/gff file\n
    \t--fasta\t/path/to/fasta genome file\n
    \t-t, --threads N\n\t\tNb threads used to index genome\t[default: ${threads}]\n
    \t-m, --memory N\n\t\tMax Memory allowed to index genome, in bytes\t[default: ${memory}]\n
    \t--assembly_report /path/to/GRChXX_latest_assembly_report.txt\n\t\tPath to assembly report to convert contig in UCSC chromosome names\t[default: ${pathToConvert}]\n
    \t--mane /path/to/MANElistFile.txt\n\t\tList of MANE transcripts\t[default: ${manePath}]
    \n
    Option for Align mode\n
    \t-F, --fastq /path/to/fastq/\n\t\trepository of the FASTQ files\n
    \t-O, --output /path/to/output/\n\t\trepository of the output files\n
    \t-p paired-end analysis\n\t\tprocesses to paired-end analysis\t[default: ${endType}]\n
    \t-t, --threads N\n\t\tNb threads used for the alignment\t[default: ${threads}]\n
    \t-m, --memory\n\t\tMax Memory allowed for the alignment\t[default: ${memory}]\n
    \t-g, --genome /path/to/genome\n\t\tpath to the STAR genome\t[default: ${genome}]\n
    \t--STAR /path/to/STAR\n\t\tpath to the STAR executable\t[default: ${STAR}]\n
    \t--tmpDir /path/to/tmp directory of files from alignment\n
    \t--samtools /path/to/samtools\n\t\tpath to samtools executable\t[default: ${samtools}]\n
    \n
    Option for Count mode\n
    \t-B, --bam /path/to/BAM files\n
    \t-O, --output /path/to/output/\n\t\tdirectory of the output files\n
    \t--bedtools\t/path/to/bedtools/bin folder \t[default: ${bedtools}]\n
    \t-b, --BEDannot /path/to/your_annotation_file.bed\n\t\tpath to exon coordinates file (in BED format)\t[default: ${BEDrefPath}]\n
    \t--fastCount improve junction count runtime\n\t\tDecrease runtime but less sensitive in junction detection\t[default: ${fastCount}]\n
    \t-p paired-end analysis\n\t\tprocesses to paired-end analysis\t[default: ${endType}]\n
    \t--samtools /path/to/samtools\n\t\tpath to samtools executable\t[default: ${samtools}]\n
    \n
    Option for SpliceLauncher mode\n
    \t-I, --input /path/to/inputFile\n\t\tRead count matrix (.txt)\n
    \t-O, --output /path/to/output/\n\t\tDirectory to save the results\n
    \t-R, --RefSeqAnnot /path/to/RefSpliceLauncher.txt\n\t\tRefSeq annotation file name \t[default: ${spliceLaucherAnnot}]\n
    \t--transcriptList /path/to/transcriptList.txt\n\t\tSet the list of transcripts to use as reference\n
    \t--txtOut\n\t\tPrint main output in text instead of xlsx\n
    \t--bedOut\n\t\tGet the output in BED format\n
    \t--min_cov 5\n\t\tMinimal number of read supporting a junction [default= ${min_cov}]\n
    \t--Graphics\n\t\tDisplay graphics of alternative junctions (Warnings: increase the runtime)\n
    \t-n, --NbIntervals 10\n\t\tNb interval of Neg Binom (Integer) [default= ${NbIntervals}]\n
    \t--SampleNames name1|name2|name3\n\t\tSample names, '|'-separated, by default use the sample file names\n
    \tIf list of transcripts (--transcriptList):\n
    \t\t--removeOther\n\t\tRemove the genes with unselected transcripts to improve runtime\n
    \tIf graphics (--Graphics):\n
    \t\t--threshold 1\n\t\tThreshold to shown junctions (%) [default= ${threshold}]\n
    -h, --help\n\tPrint this help message and exit"

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
       STAR="$2"
       shift 2 # shift past argument and past value
       ;;

        --tmpDir)
        tmpAlignDir="--tmpDir $2"
        shift 2 # shift past argument and past value
        ;;

       --mane)
       manePath="$2"
       shift 2 # shift past argument and past value
       ;;

       --assembly_report)
       pathToConvert="`readlink -v -f $2`"
       shift 2 # shift past argument and past value
       ;;

       --samtools)
       samtools="$2"
       shift 2 # shift past argument and past value
       ;;

       --bedtools)
       bedtools="$2"
       shift 2 # shift past argument and past value
       ;;

       -b|--BEDannot)
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

       --fastCount)
       fastCount=true
       shift # shift past argument
       ;;

       -t|--threads)
       threads="$2"
       shift 2 # shift past argument and past value
       ;;

       -m|--memory)
       memory="$2"
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
       spliceLaucherAnnot="`readlink -v -f $2`"
       shift 2 # shift past argument and past value
       ;;

       --transcriptList)
       TranscriptList="`readlink -v -f $2`"
       shift 2 # shift past argument and past value
       ;;

       --txtOut)
       text="--text"
       shift 1 # shift past argument and past value
       ;;

       --bedOut)
       bedOut="--bedOut"
       shift 1 # shift past argument and past value
       ;;

       --Graphics)
       Graphics="--Graphics"
       shift 1 # shift past argument and past value
       ;;

       -n|--NbIntervals)
       NbIntervals="`readlink -v -f $2`"
       shift 2 # shift past argument and past value
       ;;

       --SampleNames)
       SampleNames="$2"
       shift 2 # shift past argument and past value
       ;;

       --removeOther)
       removeOther="--removeOther"
       shift 1 # shift past argument and past value
       ;;

       --threshold)
       threshold="$2"
       shift 2 # shift past argument and past value
       ;;

       --min_cov)
       min_cov="$2"
       shift 2 # shift past argument and past value
       ;;

       -h|--help)
       echo -e "${messageHelp}"
       exit # shift past argument and past value
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
    echo -e "###############################################"
    echo -e "####### Configure SpliceLauncher environment"
    echo -e "###############################################\n"

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

    # convert contig in UCSC chromosome names
    if [ -e ${pathToConvert} ]; then
        sed -i 's/\x0D//g' ${pathToConvert} # if file download from Windows
        echo -e "\nConvert contig in UCSC chromosome names using :\n ${pathToConvert}\n On the files:\n${gff_path}\n and\n${fasta_path}\n"
        grep -v "#" ${pathToConvert} | awk 'BEGIN{OFS="|"}{print $7,$11}' > tmp_convert.txt

        OLDIFS=$IFS
        IFS=$'\n'
        for line in $(cat tmp_convert.txt)
            do 
            IFS="|"
            read -a convert <<< $line
            IFS=$'\n'
            if [ `echo -e ${#convert[1]}` -gt 3 ]; then
                ucsc_name=${convert[1]}
                ucsc_name=${ucsc_name//\\n/""}
                echo -e "replace ${convert[0]} by ${ucsc_name}"
                sed -i "s/${convert[0]}/${ucsc_name}/g" ${gff_path}
                sed -i "s/${convert[0]}/${ucsc_name}/g" ${fasta_path}
            fi
        done
        IFS=$OLDIFS

        rm tmp_convert.txt

    else
        echo -e "No convert asking"
    fi

    # run generateSpliceLauncherDB
    mkdir -p ${out_path}
    echo "Will run generateSpliceLauncherDB."
    cmd="${Rscript} ${scriptPath}/generateSpliceLauncherDB.r -i ${gff_path} -o ${out_path} --mane ${manePath}"
    echo -e "$cmd"
    echo -e "${conf_file}"
    $cmd
    BEDrefPath=${out_path}/BEDannotation.bed
    spliceLaucherAnnot=${out_path}/SpliceLauncherAnnot.txt
    SJDBannot=${out_path}/SJDBannotation.sjdb

    echo -e "${conf_file}"

    # Test if output files exist
    for i in BEDrefPath spliceLaucherAnnot SJDBannot; do
        if [[ $(test_file_if_exist "${!i}") -ne 0 ]]; then
            in_error=1
            echo_on_stderr "${i} not found! Will abort."
        fi
    done
    if [ $in_error -eq 1 ]; then
        echo -e "=> Aborting."
        exit
    fi
    
    echo -e "${conf_file}"

    genome="${out_path}/STARgenome"
    sed -i "s#^genome=.*#genome=\"${genome}\"#" ${conf_file}
    sed -i "s#^BEDrefPath=.*#BEDrefPath=\"${BEDrefPath}\"#" ${conf_file}
    sed -i "s#^spliceLaucherAnnot=.*#spliceLaucherAnnot=\"${spliceLaucherAnnot}\"#" ${conf_file}
    sed -i "s#^SJDBannot=.*#SJDBannot=\"${SJDBannot}\"#" ${conf_file}

    genomeSAsparseD=`echo "31000000000 / ${memory}" | bc `
    genomeSAsparseD=`echo "${genomeSAsparseD} + 1" | bc `

    mkdir -p ${genome}
    cmd="${STAR} \
    --runMode genomeGenerate \
    --runThreadN ${threads} \
    --genomeSAsparseD ${genomeSAsparseD} \
    --limitGenomeGenerateRAM ${memory} \
    --genomeDir ${genome} \
    --genomeFastaFiles ${fasta_path} \
    --sjdbFileChrStartEnd ${SJDBannot} \
    --sjdbGTFfile ${gff_path} \
    --sjdbOverhang 99"
    echo -e "Running STAR = $cmd"
    $cmd
fi

########## launch RNAseq
if [[ ${align} = "TRUE" ]]; then
    echo -e "###############################################"
    echo -e "####### Launch aligment step"
    echo -e "###############################################\n"
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
    cmd="${scriptPath}/pipelineRNAseq.sh --runMode Align -F ${fastq_path} -O ${out_path} -g ${genome} --STAR ${STAR} --samtools ${samtools} -t ${threads} -m ${memory} ${endType} ${tmpAlignDir}"
    echo -e "cmd = $cmd"
    $cmd

    bam_path="${out_path}/Bam"
fi

########## count RNAseq
if [[ ${count} = "TRUE" ]]; then
    echo -e "###############################################"
    echo -e "####### Launch counting step"
    echo -e "###############################################\n"

## launch counting

    # Test if files exist
    for i in bam_path; do
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

    # run count
    mkdir -p ${out_path}
    echo "Will run counting."
    if ${fastCount}; then
        cmd="${scriptPath}/pipelineRNAseq.sh --runMode Count -B ${bam_path} -O ${out_path} --bedannot ${BEDrefPath} --bedtools ${bedtools} --fastCount --perlscript ${scriptPath} "
    else
        cmd="${scriptPath}/pipelineRNAseq.sh --runMode Count -B ${bam_path} -O ${out_path} --bedannot ${BEDrefPath} --bedtools ${bedtools} --perlscript ${scriptPath} --samtools ${samtools} ${endType}"
    fi
    echo -e "cmd = $cmd"
    $cmd

    input_path="${out_path}/$(basename ${out_path}).txt"
fi

if [[ ${spliceLauncher} = "TRUE" ]]; then
    echo -e "###############################################"
    echo -e "####### Launch SpliceLauncher step"
    echo -e "###############################################\n"


    mkdir -p ${out_path}
    echo "Will run SpliceLauncher"

    if [ -z ${TranscriptList+x} ]; then transcriptList_cmd=""; else transcriptList_cmd="--TranscriptList ${TranscriptList}"; fi
    if [ -z ${SampleNames+x} ]; then SampleNames_cmd=""; else SampleNames_cmd="--SampleNames ${SampleNames}"; fi


    cmd="${Rscript} ${scriptPath}/SpliceLauncherAnalyse.r --input ${input_path} -O ${out_path} --RefSeqAnnot ${spliceLaucherAnnot} -n ${NbIntervals} ${transcriptList_cmd} ${SampleNames_cmd} ${removeOther} ${text} ${bedOut} ${Graphics} --threshold ${threshold} --min_cov ${min_cov}"
    echo -e "cmd = $cmd"
    $cmd

fi
