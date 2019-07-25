#!/bin/bash

messageHelp="Usage: $0\n

    [Mandatory] \n
        \t-F, --fastq /path/to/fastq/\n\t\trepository of the FASTQ files\n
        \t-O, --output /path/to/output/\n\t\trepository of the output files\n

    [Options]\n
        \t-p paired-end analysis\n\t\tprocesses to paired-end analysis [default: ${endType}]\n
        \t-t, --threads N\n\t\tNb threads used for the alignment[default: ${threads}]\n
        \t--bedannot /path/to/bedannot\n\t\tpath to exon coordinates (in BED format) [default: ${BEDrefPath}]\n
        \t--perlscript /path/to/perlscript/\n\t\trepository of perl scripts used by the pipeline [default: ${ScriptPath}]\n\n
   You could : $0 -F /patho/to/FASTQfolder/ -O /path/to/output/ -g /path/to/STARgenome/ --star /path/to/STAR
   --samtools /path/to/samtools --bedtools /path/to/BEDtools/bin/ --bedannot /path/to/BEDannot"

# if [ $# -lt 8 ]; then
    # echo -e $messageHelp
    # exit
# fi

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
source ./config.cfg

in_error=0 # variable initialization; will be 1 if a file or cmd not exist

# Test if cmd exist
for i in samtools bedtools STAR R perl; do
    test_command_if_exist ${!i}
    prev_state=$?
    if [ $prev_state -eq 1 ]; then
        in_error=1
    fi
done

# Test if files exist
for i in fasta_genome gff readme; do
    test_file_if_exist ${!i}
    prev_state=$?
    if [ $prev_state -eq 1 ]; then
        in_error=1
    fi
done

# exit if there is one error or more
if [ $in_error -eq 1 ]; then
    echo -e "=> Aborting."
    exit
fi



































