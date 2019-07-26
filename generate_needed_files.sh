#!/bin/bash

messageHelp="Usage: $0\n
    \n
    [Mandatory]\n
        \t-O, --output\t/path/to/output/\n\t\t\tdirectory of the output files\n
    \n
    [Optional, in Generate mode]\n
        \t-C, --config\t/path/to/configuration file/\n\t\t\t [default: config.cfg]\n
    \n
    [Optional, in Run mode]\n
        \t-F, --fastq\t/path/to/fastq/\n\t\t\tdirectory of the FASTQ files\n
        \t-p\tpaired-end analysis\n\t\t\tprocesses to paired-end analysis [default: ${endType}]\n
    \n
    [General options]\n
        \t-t, --threads\tN\n\t\t\tNb threads used for the alignment[default: ${threads}]\n
        \t--bedannot\t/path/to/bedannot\n\t\t\tpath to exon coordinates (in BED format) [default: ${BEDrefPath}]\n
        \t--perlscript\t/path/to/perlscript/\n\t\t\tlocalization of perl scripts used by the pipeline [default: ${ScriptPath}]\n
    \n
   You could : $0 -F /patho/to/FASTQfolder/ -O /path/to/output/ -g /path/to/STARgenome/ --star /path/to/STAR
   --samtools /path/to/samtools --bedtools /path/to/BEDtools/bin/ --bedannot /path/to/BEDannot"

if [ $# -lt 1 ]; then
    echo -e $messageHelp
    exit
fi

workFolder=$(readlink -f $(dirname $0))
scriptPath="${workFolder}/scripts"

while [[ $# -gt 0 ]]; do
   key=$1
   case $key in

       -F|--fastq)
       pathToFastq="`readlink -v -f $2`"
       shift 2 # shift past argument and past value
       ;;

       -O|--output)
       out="`readlink -v -f $2`"
       shift 2 # shift past argument and past value
       ;;

        -g|--genome)
        genomeDirectory="`readlink -v -f $2`"
        shift 2 # shift past argument and past value
        ;;
        --star)
        STARPath="`readlink -v -f $2`"
        shift 2 # shift past argument and past value
        ;;

        --samtools)
        SamtoolsPath="`readlink -v -f $2`"
        shift 2 # shift past argument and past value
        ;;

        --bedtools)
        BEDtoolsPath="`readlink -v -f $2`"
        shift 2 # shift past argument and past value

        ;;
        --bedannot)
        BEDrefPath="`readlink -v -f $2`"
        shift 2 # shift past argument and past value
        ;;

        --perlscript)
        ScriptPath="`readlink -v -f $2`"
        shift 2 # shift past argument and past value
        ;;

       -t|--threads)
       threads="$2"
       shift 2 # shift past argument and past value
       ;;

       -p)
       endType="paired"
       shift # shift past argument
       ;;

       -C|--config)
       conf_file="$2"
       shift 2 # shift past argument and past value
       ;;
       
       *)  # unknown option
       POSITIONAL+=("$key") # save it in an array for later
       echo -e "    Unknown option ${key}"
       shift # shift past argument
      ;;
   esac
done

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

########## parsing config file

echo "Parsing configuration file..."
test_file_if_exist ${workFolder}/config.cfg
prev_state=$?
if [ $prev_state -eq 1 ]; then
        exit 1
    fi

source ${workFolder}/config.cfg

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






























