#!/usr/bin/bash

################
#Set parameters#
################
messageHelp="Usage: $0\n
    [Mandatory] \n
        \t[-F|--fastq /path/to/fastq/] \n
        \t[-O|--output /path/to/output/] \n
        \t[-g|--genome /path/to/genome] \n
        \t[--star /path/to/STAR] \n
        \t[--samtools /path/to/samtools] \n
        \t[--bedtools /path/to/bedtoolsFolder/] \n
        \t[--bedannot /path/to/bedannot] \n
    [Options] \n
        \t[-p paired-end analysis] \n
        \t[-t|--threads <int>] \n
        \t[--perlscript /path/to/perlscript/] \n
   You could : $0 -F /patho/to/FASTQfolder/ -O /path/to/output/ -g /path/to/STARgenome/ --star /path/to/STAR
   --samtools /path/to/samtools --bedtools /path/to/BEDtools/bin/ --bedannot /path/to/BEDannot"

if [ $# -lt 8 ]; then
    echo -e $messageHelp
    exit
fi

threads="1"
endType="single"
workFolder=$(readlink -f $(dirname $0))
ScriptPath=${workFolder}/"scripts"

while [[ $# -gt 0 ]]; do
   key=$1
   case $key in

       -F|--fastq)
       pathToFastq="`readlink -v -f $2`"
       shift 2 # shift past argument and past value
       ;;

       -O|--output)
       pathToRun="`readlink -v -f $2`"
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

        --perlscirpt)
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

       *)  # unknown option
       POSITIONAL+=("$key") # save it in an array for later
       echo -e "    Unknown option ${key}"
       shift # shift past argument
#        exit # end of prog... here, bash
      ;;

   esac
done

RunPath="${pathToRun}/"
runName=`basename $pathToRun`

mkdir -p ${pathToRun}

echo "####Your options:"
echo -e "    path to fastq = ${pathToFastq}"
echo -e "    path to output = ${pathToRun}"
echo -e "    path to STAR genome files = ${genomeDirectory}"
echo -e "    path to STAR = ${STARPath}"
echo -e "    path to samtools = ${SamtoolsPath}"
echo -e "    path to bedtools = ${BEDtoolsPath}"
echo -e "    path to the bedtools annot file = ${BEDrefPath}"
echo -e "    path to perl scripts = ${ScriptPath}"
echo -e "    name of run = ${runName}"
echo -e "    threads = ${threads}"
echo -e "    end type analysis = ${endType}"

chech_error=0 #stop program if chech_error = 1

#check depository
echo "####Check your option"
for depository in ${pathToFastq} ${genomeDirectory}; do
    if [ ! -d $depository ]
        then
        echo -e "****You must define the depository : ${depository}****"
        ${chech_error}=1
    else
        echo "    The depository: ${depository}... OK"
    fi
done

#check BED exon annotation
if [ ! -e $BEDrefPath ]
then
   echo -e "****We cannot found the bedtools annot file at ${BEDrefPath}****"
   ${chech_error}=1
else
    echo "    The bedtools annot file: ${BEDrefPath}... OK"
fi

#check the softs
for soft in ${SamtoolsPath} \
 ${STARPath} \
 "${ScriptPath}/bedBlocks2IntronsCoords.pl" \
 "${ScriptPath}/joinJuncFiles.pl" \
 "${BEDtoolsPath}/bamToBed" \
 "${BEDtoolsPath}/closestBed" ; do
    command -v ${soft} >/dev/null 2>&1 || { ${chech_error}=1; echo >&2 "I require ${soft} but it's not installed. Aborting.";}
done

if [ ${chech_error} -eq 1 ]; then
    echo -e $messageHelp
    exit 1;
else
    echo "    The Samtools soft: ${SamtoolsPath}... OK"
    echo "    The STAR soft: ${STARPath}... OK"
    echo "    The perl Scripts: ${ScriptPath}... OK"
    echo "    The BEDtools soft: ${BEDtoolsPath}... OK"
fi

###########################
# Alignement STAR v.2.6.0 #
###########################

##Parametres choisis##
# outSAMstrandField : intronMotif
# outFilterMismatchNmax : Nombre de mismatch autoris�s "2"
# outFilterMultimapNmax : 10
# runThreadN : Nombre de coeurs � allouer "2"
# outSAMunmapped : les reads non mapp�s ne sont pas "transform�s" en fichier sam "Within"
# outStd : Fichier de sortie ".sam" (NB : dans la nouvelle version de STAR, la sortie peut �tre en bam)
# genomeLoad : Gestion de la m�moire lors du chargement du g�nome "LoadAndKeep"

echo "###### Aligment process ######"

mkdir -p $RunPath

# Procedure
BamPath="${RunPath}Bam"
mkdir -p ${BamPath}

cd ${pathToFastq}

ls $pathToFastq

# unload genome reference of any last run
${STARPath} --genomeDir ${genomeDirectory} --genomeLoad Remove

if [ $endType = "paired" ]
then
    #paired-end
    for file in *1.fastq.gz; do
        if [ -e $file ]
        then
            echo "traitement de $file..."
            ${STARPath} \
                --outSAMstrandField intronMotif \
                --outFilterMismatchNmax 2 \
                --outFilterMultimapNmax 10 \
                --genomeDir ${genomeDirectory} \
                --readFilesIn $file ${file%1.fastq.gz}2.fastq.gz\
                --readFilesCommand zcat \
                --runThreadN ${threads} \
                --outSAMunmapped Within \
                --outSAMtype BAM SortedByCoordinate \
                --limitBAMsortRAM 15000000000 \
                --outSAMheaderHD \@HD VN:1.4 SO:SortedByCoordinate \
                --outFileNamePrefix ${BamPath}/${file%1.fastq.gz}. \
                --genomeLoad LoadAndKeep > ${BamPath}/${file%1.fastq.gz}.log 2>&1
        else
            echo -e "****the file $file dosn't exist, Please rename your file in XXXXX1(2).fastq.gz format****\n****XXXXX1.fastq.gz for R1 read and XXXXX2.fastq.gz for R2 read****"
            echo -e $messageHelp
            exit
        fi
    done
else
    #single-end
    for file in *.fastq.gz; do
        if [ -e $file ]
        then
            echo "traitement de $file..."
            ${STARPath} \
                --outSAMstrandField intronMotif \
                --outFilterMismatchNmax 2 \
                --outFilterMultimapNmax 10 \
                --genomeDir ${genomeDirectory} \
                --readFilesIn $file \
                --readFilesCommand zcat \
                --runThreadN ${threads} \
                --outSAMunmapped Within \
                --outSAMtype BAM SortedByCoordinate \
                --limitBAMsortRAM 15000000000 \
                --outSAMheaderHD \@HD VN:1.4 SO:SortedByCoordinate \
                --outFileNamePrefix ${BamPath}/${file%.fastq.gz}. \
                --genomeLoad LoadAndKeep > ${BamPath}/${file%.fastq.gz}.log 2>&1
        else
            echo "****the file $file don't exist, Please rename your file in XXXXX.fastq.gz****"
            echo -e $messageHelp
            exit
        fi
    done
fi
# unload genome reference from shared memory
${STARPath} --genomeDir ${genomeDirectory} --genomeLoad Remove

## Transformation des fichers sam en bam, trie, cr�ation de l'index
echo "###### Make BAM index ######"

cd $BamPath
for file in *sortedByCoord.out.bam; do
    echo "traitement de $file..."
    ${SamtoolsPath} index $file
done

#######################################################
#Creation fichier BAM sens et anti-sens et fichier BED#
#######################################################
echo "###### Make BED file and junction count ######"

ClosestExPath="${RunPath}getClosestExons"
mkdir -p ${ClosestExPath}

cd $BamPath

if [ $endType = "paired" ]
then
    for file in *.Aligned.sortedByCoord.out.bam; do
        echo "traitement de $file for forward read..."
        ${SamtoolsPath} view -b -f 0x40 $file | \
            ${BEDtoolsPath}/bamToBed -bed12 -i stdin | \
            awk '{if($10>1){print $0}}' | \
            perl ${ScriptPath}/bedBlocks2IntronsCoords.pl y - | \
            awk '{if($5==255){print $0}}'  > ${file%.Aligned.sortedByCoord.out.bam}_juncs.bed
        echo "traitement de $file for reverse read..."
        ${SamtoolsPath} view -b -f 0x80 $file | \
            ${BEDtoolsPath}/bamToBed -bed12 -i stdin | \
            awk '{if($10>1){print $0}}' | \
            perl ${ScriptPath}/bedBlocks2IntronsCoords.pl n - | \
            awk '{if($5==255){print $0}}'  >> ${file%.Aligned.sortedByCoord.out.bam}_juncs.bed
    done

    echo "###### count junctions ######"
    for file in *_juncs.bed; do
        echo "traitement de $file..."
        sort -k1,1 -k2,2n $file | \
            uniq -c | \
            awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5 "_" $1,$6,$7}' | \
            ${BEDtoolsPath}/closestBed -d -t first -a stdin -b ${BEDrefPath} | \
            awk 'BEGIN{OFS="\t"}{if($13<200){print $0}else{print $1,$2,$3,$4,$5,$6,".","-1","-1",".","-1",".","-1" }}' | \
            awk 'BEGIN{OFS="\t"}{if($8>0){split($4,counts,"_"); split($10,nm,"_");print $1,$2,$3,$6,nm[1] "_" nm[2],counts[1],counts[2],counts[3],counts[4]}}' \
            >  ${ClosestExPath}/${file%_juncs.bed}.txt
    done
else
    for file in *.Aligned.sortedByCoord.out.bam; do
        echo "traitement de $file..."
        ${BEDtoolsPath}/bamToBed -bed12 -i $file | \
            awk '{if($10>1){print $0}}' | \
            perl ${ScriptPath}/bedBlocks2IntronsCoords.pl y - | \
            awk '{if($5==255){print $0}}' | \
            sort -k1,1 -k2,2n | \
            uniq -c | \
            awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5 "_" $1,$6,$7}' | \
            tee ${ClosestExPath}/${file%.Aligned.sortedByCoord.out.bam}.sort_count.bed | \
            ${BEDtoolsPath}/closestBed -d -t first -a stdin -b ${BEDrefPath} | \
            awk 'BEGIN{OFS="\t"}{if($13<200){print $0}else{print $1,$2,$3,$4,$5,$6,".","-1","-1",".","-1",".","-1" }}' | \
            awk 'BEGIN{OFS="\t"}{if($8>0){split($4,counts,"_"); split($10,nm,"_");print $1,$2,$3,$6,nm[1] "_" nm[2],counts[1],counts[2],counts[3],counts[4]}}' \
            >  ${ClosestExPath}/${file%.Aligned.sortedByCoord.out.bam}.txt
    done
fi

#######################
#fin comptage jonction#
#######################
echo "###### Merge the junction counts ######"

cd $ClosestExPath

perl ${ScriptPath}/joinJuncFiles.pl *.txt > ${RunPath}/${runName}.txt

echo "###### Finished! ######"
