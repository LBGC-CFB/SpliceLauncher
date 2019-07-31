#!/usr/bin/env Rscript
#########################
#SpliceLauncher
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

options(stringsAsFactors=FALSE)

tryCatch({
library(WriteXLS)
},
	error=function(cond) {
		message("Here's the original error message:")
		message(cond)
		message("*****You need to install \'WriteXLS\' library")
})

tryCatch({
library(Cairo)
},
	error=function(cond) {
		message("Here's the original error message:")
		message(cond)
		message("*****You need to install \'Cairo\' library")
})
message("SpliceLauncherV1.3")

#####################
#import of arguments#
#####################
#5AS => Donor shift
#3AS => Acceptor shift
#SkipEx => exon skipping

#default Parameters
EchName=NULL
pathTranscript=NULL
checkBEDannot="NO"
DisplayGraph='NO'
removeOther='NO'
thr=1
StatAnalysis='YES'
negbinom.n=10
printText=FALSE

argsFull <- commandArgs()

scriptPath=dirname(normalizePath(sub("--file=","",argsFull[substr(argsFull,1,7)=="--file="])))
RefFile= paste(scriptPath,"refData/RefSpliceLauncher.txt",sep="/")
if(!file.exists(RefFile)){RefFile="***WARNING: REFSEQ FILE NOT FIND***"}

helpMessage=paste("Usage: SpliceLauncher.r\n
    [Mandatory] \n
        -I, --input /path/to/inputFile\n\t\tRead count matrix (.txt)
        -O, --output /path/to/output/\n\t\tDirectory to save the results\n
    [Options] \n
        -R, --RefSeqAnnot /path/to/RefSpliceLauncher.txt\n\t\tRefSeq annotation file name [default=",RefFile,"]
        --TranscriptList /path/to/transcriptList.txt\n\t\tSet the list of transcripts to use as reference
        --text\n\t\tPrint main output in txt instead of excel
        -b, --BEDannot\n\t\tGet the output in BED format
        --Graphics\n\t\tDisplay graphics of alternative junctions (Warnings: increase the runtime)
        -n, --NbIntervals 10\n\t\tNb interval of Neg Binom (Integer) [default=",negbinom.n,"]
        --SampleNames name1|name2|name3\n\t\tSample names, '|'-separated, by default use the sample file names\n
        If list of transcripts (--TranscriptList):
            --removeOther\n\t\tRemove the genes with unselected transcripts to improve runtime
        If graphics (-g, --Graphics):
            --threshold 1\n\t\tThreshold to shown junctions (%) [default=",thr,"]
    h, --help\n\t\tprint this help message and exit\n
   You could : Rscript SpliceLauncher.r -I ./dataTest/MatrixCountExample.txt -O ./outputTest/")

#get script argument
if (length(which(argsFull=="--args"))==0){message(helpMessage);q(save = "no")}

args = argsFull[(which(argsFull=="--args")+1):length(argsFull)]
if (length(args)<3){message(helpMessage);q(save = "no")}

i=1
while (i <= length(args)){
    if(args[i]=="-I"|args[i]=="--input"){
        inputFile=normalizePath(path=args[i+1]);i = i+2
        run=sub(".txt","",basename(inputFile))
    }else if(args[i]=="-R"|args[i]=="--RefSeqAnnot"){
        RefFile=normalizePath(args[i+1]);i = i+2
    }else if(args[i]=="--TranscriptList"){
        pathTranscript=normalizePath(args[i+1]);i = i+2
    }else if(args[i]=="-O"|args[i]=="--output"){
        outputDir=args[i+1];i = i+2
    }else if(args[i]=="--SampleNames"){
        EchName=args[i+1];i = i+2
    }else if(args[i]=="-b"|args[i]=="--BEDannot"){
        checkBEDannot="YES";i = i+1
    }else if(args[i]=="--Graphics"){
        DisplayGraph="YES";i = i+1
    }else if(args[i]=="--text"){
        printText=TRUE;i = i+1
    }else if(args[i]=="--threshold"){
        thr=args[i+1];i = i+2
    }else if(args[i]=="--removeOther"){
        removeOther="YES";i = i+1
    }else if(args[i]=="-n"|args[i]=="--NbIntervals"){
        negbinom.n=args[i+1];i = i+2
    }else if(args[i]=="-h"|args[i]=="--help"){
        message(helpMessage);q(save="no")
    }else{
        message(helpMessage);stop(paste("********Unknown option:",args[i],"\n"))
    }
}

message("######################")
message('#Parameters of SpliceLauncher')
message("######################")

if (!file.exists(outputDir)){dir.create(outputDir, recursive = TRUE)}
outputDir = normalizePath(outputDir)

if(!is.null(pathTranscript) & removeOther=="NO"){
    cat("You defined a transcript list\n Do you want to remove other transcript ? [y/n] ");
    removeOther <- readLines("stdin",n=1);
    removeOther = if(removeOther=="y"|removeOther=="Y"){"YES"}
}

message(paste("Matrix count:", inputFile))
message(paste("Output directory:",outputDir))
message(paste("Run name:",run))
message(paste("Sample Name:",EchName))
message(paste("RefSeq Annot:",RefFile))
message(paste("Output format:",if(printText){"text"}else{"excel"}))
message(paste("List of transcripts:",pathTranscript))
message(paste("Remove the other genes:",removeOther))
message(paste("BED annotation:",checkBEDannot))
message(paste("Display graphics:",DisplayGraph))
message(paste("Threshold:",thr))
message(paste("nb intervals for Neg Binom:",negbinom.n))

#Creating folders to save the output

chem_results = paste(outputDir,"/",run,"_results/",sep="")
dir.create(chem_results,showWarnings=FALSE)

#Path to report

chem_rapport=paste(outputDir,"/",sep="")
date_analyse=format(Sys.time(), "%a %b %d %Y %X")

message("######################")
message("#Data formating...")
message("######################")

#import of data

T1<-Sys.time()
message("   Import matrix count...")
tmp=read.table(file=inputFile, sep="\t", header=T)

SampleInput=names(tmp)[c(7:ncol(tmp))]

if(!is.null(EchName)){
    EchName = unlist(strsplit(EchName,"|",fixed=TRUE))
	if(length(SampleInput)!=length(EchName)){
		stop("***** All sample were not named")
	}else{
        SampleInput = EchName
        names(tmp)[c(7:ncol(tmp))] = SampleInput
    }
}else{
	EchName = SampleInput
}

tmp$Conca=paste(tmp$chr,tmp$start,tmp$end,sep="_")

message("   Download RefSeq File...")

tableConvert=read.table(file = RefFile, sep="\t", header=T)

message("   Get gene name...")

tableTomerge = tableConvert[,c("Gene","transcrit","Strand")]
tableTomerge = tableTomerge[-which(duplicated(tableTomerge$transcrit)),]

tmp <- merge(tmp, tableTomerge, by.x="NM", by.y="transcrit", all.x=TRUE)

missingNM = tmp[is.na(tmp$Gene),"NM"]
missingNM = missingNM[-which(duplicated(missingNM))]
message(paste("   I don't find",length(missingNM),"transcrit(s) (",paste(missingNM,collapse = ", "),")"))

tmp = tmp[!is.na(tmp$Gene),]

tmp$Strand_transcript[tmp$Strand=="+"] = "forward"
tmp$Strand_transcript[tmp$Strand=="-"] = "reverse"

getTransMaj <- function (x){
	trans = names(which(x==max(x)))
	return(trans[1])
}

if(!is.null(pathTranscript)){
    message('   Apply transcripts list...')
    SelectTranscrit = readLines(pathTranscript)
    tableSelect = tableTomerge[which(tableTomerge$transcrit%in%SelectTranscrit),c("Gene","transcrit")]
    if(max(as.numeric(table(tableSelect[,"Gene"])))>1){
        stop(paste("***** Multiple reference transcripts found for the",
                names(table(tableSelect[,"Gene"]))[as.numeric(table(tableSelect[,"Gene"]))>1]))
    }
    if(removeOther=="YES"){
        tmp = tmp[which(tmp$Gene%in%tableSelect[,"Gene"]),]
    }
    for(gene in tableSelect[,"Gene"]){
        tmp$NM[tmp$Gene==gene] <- tableSelect[tableSelect$Gene==gene,"transcrit"]
    }
}

message('   Merge transcrit...')
message(length(unique(tmp$NM)))
MatGeneTrans = table(tmp$NM,tmp$Gene)
TransMaj = apply(MatGeneTrans,2,getTransMaj)
convertTransMaj = data.frame(Gene = names(TransMaj),NMadjust = TransMaj,row.names=1:length(TransMaj))
tmp <- merge(tmp, convertTransMaj, by="Gene")
tmp$NM = tmp$NMadjust
message(length(unique(tmp$NM)))

message("   Generate matrix for SpliceLauncher calculation...")
data_junction = tmp[,c("Conca","chr","start","end","strand","Strand_transcript","NM","Gene",SampleInput)]

message("   Remove Opposite junctions...")

idDoublon = c(which(duplicated(data_junction$Conca,fromLast=T)),which(duplicated(data_junction$Conca)))
data_junctionUnique = data_junction[-idDoublon,]
data_junctionDoublon = data_junction[idDoublon,]

JuncNoFilt = nrow(data_junctionDoublon)
data_junctionDoublon = data_junctionDoublon[(data_junctionDoublon$Strand_transcript=="forward" & data_junctionDoublon$strand=="+")|
						(data_junctionDoublon$Strand_transcript=="reverse" & data_junctionDoublon$strand=="-"),]
JuncFilt = nrow(data_junctionDoublon)
nbJuncOpp = JuncNoFilt-JuncFilt

data_junction = rbind(data_junctionUnique,data_junctionDoublon)

message(paste("   Nb junctions removed:",nbJuncOpp))

#index of column containing read counts
input= 9:dim(data_junction)[2]

n_ech=length(input)

if(n_ech<5){
	message('***Not enough samples to perform statistical analysis\n    Param for statistical analysis is fixed to \'NO\'')
	StatAnalysis='NO'
}else if(n_ech>=5 & n_ech<10){
	message('***WARNINGS: low number of samples statistical analysis may be not relevant')
}

#index of column to the output
SampleOutput = paste("P",SampleInput,sep="_")
eval(parse(text=paste("data_junction$",SampleOutput," <- 0",sep="",collapse=" ; ")))
output = c((dim(data_junction)[2]-(n_ech-1)):dim(data_junction)[2])

jid.column<-"Conca"  # column with junction id
sid<-names(data_junction)[output]  # column with individual data to model

#################
#index of genes#
#################
message("   Index of gene...")

data_junction=data_junction[order(data_junction$Gene),]

data_junction$Gene = data_junction$Gene
data_junction$Gene = as.factor(data_junction$Gene)

data_junction$ID_gene=as.numeric(data_junction$Gene)

nb_gene=max(data_junction$ID_gene)
nb_gene
#####################
#junction annotation#
#####################

BlackInvTriangle = intToUtf8(0x25BC)
DeltaSymb = intToUtf8(0x2206)

annotjunc = ""
constitutive = ""
calcul = ""

#Name Junctions

getAnnotJuncs <- function(Start, End){
	constitutive = "NoPhysio"
	calcul = "NoData"
	nameJunc = "no annotate"

	if(nrow(tableConvertAnnot)>=1){
	strand = tableConvertAnnot$Strand[1]
	# Nearest exon to the start position
		NearestStart = tableConvertAnnot$gStart[which(abs(tableConvertAnnot$gStart - Start) == min(abs(tableConvertAnnot$gStart-  Start)))]
		NearestEnd = tableConvertAnnot$gEnd[which(abs(tableConvertAnnot$gEnd - Start) == min(abs(tableConvertAnnot$gEnd - Start)))]
		NearestStart = NearestStart[1]
		NearestEnd = NearestEnd[1]
		if(abs(NearestStart-Start)<abs(NearestEnd-Start)){
			NearestEnd = tableConvertAnnot$gEnd[tableConvertAnnot$gStart==NearestStart]
		}else if(abs(NearestStart-Start)>abs(NearestEnd-Start)){
			NearestStart = tableConvertAnnot$gStart[tableConvertAnnot$gEnd==NearestEnd]
		}
		NearestExonStart = tableConvertAnnot$idEx[tableConvertAnnot$gStart==NearestStart]

	# Nearest Exon to the end position
		NearestStart = tableConvertAnnot$gStart[which(abs(tableConvertAnnot$gStart - End) == min(abs(tableConvertAnnot$gStart-  End)))]
		NearestEnd = tableConvertAnnot$gEnd[which(abs(tableConvertAnnot$gEnd - End) == min(abs(tableConvertAnnot$gEnd - End)))]
		NearestStart = NearestStart[1]
		NearestEnd = NearestEnd[1]
		if(abs(NearestStart-End)<abs(NearestEnd-End)){
			NearestEnd = tableConvertAnnot$gEnd[tableConvertAnnot$gStart==NearestStart]
		}else if(abs(NearestStart-End)>abs(NearestEnd-End)){
			NearestStart = tableConvertAnnot$gStart[tableConvertAnnot$gEnd==NearestEnd]
		}
		NearestExonEnd = tableConvertAnnot$idEx[tableConvertAnnot$gStart==NearestStart]

		if(strand=="+"){
			transcritStart = min(tableConvertAnnot$gStart)
			transcritEnd = max(tableConvertAnnot$gEnd)
			if(Start < transcritStart | transcritEnd < End){
				nameJunc = "Outside Transcript"
			}else{

				if(tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonStart]==Start){
					if(tableConvertAnnot$gStart[tableConvertAnnot$idEx==NearestExonEnd]==End){
						if((NearestExonEnd-NearestExonStart)==1){
							nameJunc = paste(NearestExonStart,NearestExonEnd,sep="_")
							constitutive="Physio"
							calcul="Physio"
						}else if((NearestExonEnd-NearestExonStart)==2){
							nameJunc = paste(DeltaSymb, (NearestExonStart + 1), sep="")
							calcul="SkipEx"
						}else if((NearestExonEnd-NearestExonStart)>2){
							nameJunc = paste(DeltaSymb, (NearestExonStart + 1),"_", (NearestExonEnd - 1), sep="")
							calcul="SkipEx"
						}else{
							nameJunc = "exon excision"
						}
					}else{
						if((NearestExonEnd-NearestExonStart)==1){
							if(End < tableConvertAnnot$gStart[tableConvertAnnot$idEx==NearestExonEnd]){
								nameJunc = paste(BlackInvTriangle,NearestExonStart,"p","(",
												abs(End - tableConvertAnnot$gStart[tableConvertAnnot$idEx==NearestExonEnd]),")",sep="")
							}else if(End < tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonEnd]){
								nameJunc = paste(DeltaSymb,NearestExonEnd,"p","(",
												abs(End - tableConvertAnnot$gStart[tableConvertAnnot$idEx==NearestExonEnd]),")",sep="")
							}else{
								nameJunc = paste(BlackInvTriangle,NearestExonEnd,"p","(",abs(End - tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonEnd]),
												"),",DeltaSymb,NearestExonEnd,abs(End - tableConvertAnnot$gStart[tableConvertAnnot$idEx==NearestExonEnd]),sep="")
							}
						}else if((NearestExonEnd-NearestExonStart)==2){
							if(End < tableConvertAnnot$gStart[tableConvertAnnot$idEx==NearestExonEnd]){
								nameJunc = paste(BlackInvTriangle,NearestExonStart,"p","(",
												abs(End - tableConvertAnnot$gStart[tableConvertAnnot$idEx==NearestExonEnd]),")",",",
												DeltaSymb ,NearestExonStart+1,sep="")
							}else if(End < tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonEnd]){
								nameJunc = paste(DeltaSymb,NearestExonStart+1,",",NearestExonEnd,"p","(",
												abs(End - tableConvertAnnot$gStart[tableConvertAnnot$idEx==NearestExonEnd]),")",sep="")
							}else{
								nameJunc = paste(BlackInvTriangle,NearestExonEnd,"p","(",abs(End - tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonEnd]),
												"),",DeltaSymb,NearestExonStart+1, "_",NearestExonEnd-1,sep="")
							}
						}else if((NearestExonEnd-NearestExonStart) > 2){
							if(End < tableConvertAnnot$gStart[tableConvertAnnot$idEx==NearestExonEnd]){
								nameJunc = paste(BlackInvTriangle,NearestExonEnd-1,"p","(",
												abs(End - tableConvertAnnot$gStart[tableConvertAnnot$idEx==NearestExonEnd]),")",",",
												DeltaSymb ,NearestExonStart+1,"_",NearestExonEnd-1,sep="")
							}else if(End < tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonEnd]){
								nameJunc = paste(DeltaSymb,NearestExonStart+1,"_",NearestExonEnd-1,",",NearestExonEnd,"p","(",
												abs(End - tableConvertAnnot$gStart[tableConvertAnnot$idEx==NearestExonEnd]),")",sep="")
							}else{
								nameJunc = paste(BlackInvTriangle,NearestExonEnd,"p","(",abs(End - tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonEnd]),
												"),",DeltaSymb,NearestExonStart+1, "_",NearestExonEnd-1,sep="")
							}
						}else if((NearestExonEnd-NearestExonStart) == 0){
							nameJunc = paste(BlackInvTriangle,NearestExonStart,"p","(",
											abs(End - tableConvertAnnot$gStart[tableConvertAnnot$idEx==NearestExonStart+1]),")",sep="")
						}
						calcul="3AS"
					}
				}else if(tableConvertAnnot$gStart[tableConvertAnnot$idEx==NearestExonEnd]==End){
					if((NearestExonEnd-NearestExonStart)==1){
						if(Start < tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonStart]){
							nameJunc = paste(DeltaSymb,NearestExonStart,"q","(",
											abs(Start - tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonStart]),")",sep="")
						}else if(Start < tableConvertAnnot$gStart[tableConvertAnnot$idEx==NearestExonStart]){
							nameJunc = paste(BlackInvTriangle,NearestExonStart-1,"q","(",
											abs(Start - tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonStart-1]),
											")",DeltaSymb,NearestExonStart,sep="")
						}else if(Start > tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonStart]){
							nameJunc = paste(BlackInvTriangle,NearestExonStart,"q","(",
											abs(Start - tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonStart]),")",sep="")
						}
					}else if((NearestExonEnd-NearestExonStart)==2){
						if(Start < tableConvertAnnot$gStart[tableConvertAnnot$idEx==NearestExonStart]){
							nameJunc = paste(BlackInvTriangle,NearestExonStart-1,"q","(",
											abs(Start - tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonStart-1]),
											")",DeltaSymb,NearestExonStart,"_",NearestExonEnd-1,sep="")
						}else if(Start < tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonStart]){
							nameJunc = paste(DeltaSymb,NearestExonStart,"q","(",
											abs(Start - tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonStart])
											,"),",DeltaSymb,NearestExonStart+1,sep="")
						}else if(Start > tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonStart]){
							nameJunc = paste(BlackInvTriangle,NearestExonStart,"q","(",
											abs(Start - tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonStart]),
											"),",DeltaSymb,NearestExonStart+1,sep="")
						}
					}else if ((NearestExonEnd-NearestExonStart)>2){
						if(Start < tableConvertAnnot$gStart[tableConvertAnnot$idEx==NearestExonStart]){
							nameJunc = paste(BlackInvTriangle,NearestExonStart-1,"q","(",
											abs(Start - tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonStart-1]),
											")",DeltaSymb,NearestExonStart,"_",NearestExonEnd-1,sep="")
						}else if(Start < tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonStart]){
							nameJunc = paste(DeltaSymb,NearestExonStart,"q","(",
											abs(Start - tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonStart])
											,"),",DeltaSymb,NearestExonStart+1,"_",NearestExonEnd-1,sep="")
						}else if(Start > tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonStart]){
							nameJunc = paste(BlackInvTriangle,NearestExonStart,"q","(",
											abs(Start - tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonStart]),
											"),",DeltaSymb,NearestExonStart+1,"_",NearestExonEnd-1,sep="")
						}
					}else if((NearestExonEnd-NearestExonStart) == 0){
							nameJunc = paste(BlackInvTriangle,NearestExonStart-1,"q","(",
											abs(Start - tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonStart-1]),")",sep="")
					}
					calcul="5AS"
				}else{
					nameJunc = "Event too complex"
				}
			}
		}else if (strand == "-"){
			transcritStart = max(tableConvertAnnot$gStart)
			transcritEnd = min(tableConvertAnnot$gEnd)

			if(Start > transcritStart | transcritEnd > End){
				nameJunc = "Outside Transcript"
			}else{
				if(tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonStart]==Start){
					if(tableConvertAnnot$gStart[tableConvertAnnot$idEx==NearestExonEnd]==End){
						if((NearestExonEnd-NearestExonStart)==1){
							nameJunc = paste(NearestExonStart,NearestExonEnd,sep="_")
							constitutive="Physio"
							calcul="Physio"
						}else if((NearestExonEnd-NearestExonStart)==2){
							nameJunc = paste(DeltaSymb, (NearestExonStart + 1), sep="")
							calcul="SkipEx"
						}else if((NearestExonEnd-NearestExonStart)>2){
							nameJunc = paste(DeltaSymb, (NearestExonStart + 1),"_", (NearestExonEnd - 1), sep="")
							calcul="SkipEx"
						}else{
							nameJunc = "exon excision"
						}
					}else{
						if((NearestExonEnd-NearestExonStart)==1){
							if(End > tableConvertAnnot$gStart[tableConvertAnnot$idEx==NearestExonEnd]){
								nameJunc = paste(BlackInvTriangle,NearestExonStart,"p","(",
												abs(End - tableConvertAnnot$gStart[tableConvertAnnot$idEx==NearestExonEnd]),")",sep="")
							}else if(End > tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonEnd]){
								nameJunc = paste(DeltaSymb,NearestExonEnd,"p","(",
												abs(End - tableConvertAnnot$gStart[tableConvertAnnot$idEx==NearestExonEnd]),")",sep="")
							}else{
								nameJunc = paste(BlackInvTriangle,NearestExonEnd,"p","(",abs(End - tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonEnd]),
												"),",DeltaSymb,NearestExonEnd,abs(End - tableConvertAnnot$gStart[tableConvertAnnot$idEx==NearestExonEnd]),sep="")
							}
						}else if((NearestExonEnd-NearestExonStart)==2){
							if(End > tableConvertAnnot$gStart[tableConvertAnnot$idEx==NearestExonEnd]){
								nameJunc = paste(BlackInvTriangle,NearestExonStart,"p","(",
												abs(End - tableConvertAnnot$gStart[tableConvertAnnot$idEx==NearestExonEnd]),"),",
												DeltaSymb ,NearestExonStart+1,sep="")
							}else if(End > tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonEnd]){
								nameJunc = paste(DeltaSymb,NearestExonStart+1,",",NearestExonEnd,"p","(",
												abs(End - tableConvertAnnot$gStart[tableConvertAnnot$idEx==NearestExonEnd]),")",sep="")
							}else{
								nameJunc = paste(BlackInvTriangle,NearestExonEnd,"p","(",abs(End - tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonEnd]),
												"),",DeltaSymb,NearestExonStart+1, "_",NearestExonEnd-1,sep="")
							}
						}else if((NearestExonEnd-NearestExonStart) > 2){
							if(End > tableConvertAnnot$gStart[tableConvertAnnot$idEx==NearestExonEnd]){
								nameJunc = paste(BlackInvTriangle,NearestExonEnd-1,"p","(",
												abs(End - tableConvertAnnot$gStart[tableConvertAnnot$idEx==NearestExonEnd]),")",",",
												DeltaSymb ,NearestExonStart+1,"_",NearestExonEnd-1,sep="")
							}else if(End > tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonEnd]){
								nameJunc = paste(DeltaSymb,NearestExonStart+1,"_",NearestExonEnd-1,",",NearestExonEnd,"p","(",
												abs(End - tableConvertAnnot$gStart[tableConvertAnnot$idEx==NearestExonEnd]),")",sep="")
							}else{
								nameJunc = paste(BlackInvTriangle,NearestExonEnd,"p","(",abs(End - tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonEnd]),
												"),",DeltaSymb,NearestExonStart+1, "_",NearestExonEnd-1,sep="")
							}
						}else if((NearestExonEnd-NearestExonStart) == 0){
							nameJunc = paste(BlackInvTriangle,NearestExonStart,"p","(",
											abs(End - tableConvertAnnot$gStart[tableConvertAnnot$idEx==NearestExonStart+1]),")",sep="")
						}
						calcul="3AS"
					}
				}else if(tableConvertAnnot$gStart[tableConvertAnnot$idEx==NearestExonEnd]==End){
					if((NearestExonEnd-NearestExonStart)==1){
						if(Start > tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonStart]){
							nameJunc = paste(DeltaSymb,NearestExonStart,"q","(",
											abs(Start - tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonStart]),")",sep="")
						}else if(Start > tableConvertAnnot$gStart[tableConvertAnnot$idEx==NearestExonStart]){
							nameJunc = paste(BlackInvTriangle,NearestExonStart-1,"q","(",
											abs(Start - tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonStart-1]),
											")",DeltaSymb,NearestExonStart,sep="")
						}else if(Start < tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonStart]){
							nameJunc = paste(BlackInvTriangle,NearestExonStart,"q","(",
											abs(Start - tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonStart]),")",sep="")
						}
					}else if((NearestExonEnd-NearestExonStart)==2){
						if(Start > tableConvertAnnot$gStart[tableConvertAnnot$idEx==NearestExonStart]){
							nameJunc = paste(BlackInvTriangle,NearestExonStart-1,"q","(",
											abs(Start - tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonStart-1]),
											")",DeltaSymb,NearestExonStart,"_",NearestExonEnd-1,sep="")
						}else if(Start > tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonStart]){
							nameJunc = paste(DeltaSymb,NearestExonStart,"q","(",
											abs(Start - tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonStart])
											,"),",DeltaSymb,NearestExonStart+1,sep="")
						}else if(Start < tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonStart]){
							nameJunc = paste(BlackInvTriangle,NearestExonStart,"q","(",
											abs(Start - tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonStart]),
											"),",DeltaSymb,NearestExonStart+1,sep="")
						}
					}else if ((NearestExonEnd-NearestExonStart)>2){
						if(Start > tableConvertAnnot$gStart[tableConvertAnnot$idEx==NearestExonStart]){
							nameJunc = paste(BlackInvTriangle,NearestExonStart-1,"q","(",
											abs(Start - tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonStart-1]),
											")",DeltaSymb,NearestExonStart,"_",NearestExonEnd-1,sep="")
						}else if(Start > tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonStart]){
							nameJunc = paste(DeltaSymb,NearestExonStart,"q","(",
											abs(Start - tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonStart])
											,"),",DeltaSymb,NearestExonStart+1,"_",NearestExonEnd-1,sep="")
						}else if(Start < tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonStart]){
							nameJunc = paste(BlackInvTriangle,NearestExonStart,"q","(",
											abs(Start - tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonStart]),
											"),",DeltaSymb,NearestExonStart+1,"_",NearestExonEnd-1,sep="")
						}
					}else if((NearestExonEnd-NearestExonStart) == 0){
							nameJunc = paste(BlackInvTriangle,NearestExonStart-1,"q","(",
											abs(Start - tableConvertAnnot$gEnd[tableConvertAnnot$idEx==NearestExonStart-1]),")",sep="")
					}
					calcul="5AS"
				}else{
					nameJunc = "Event too complex"
				}
			}
		}
	}else{
		nameJunc = "Transcript Not found"
	}
	annot <<-c(constitutive,calcul,nameJunc)
}

#Conversion gNomen in cNomen

convertTocNomenFor <- function(gPos){

	if(nrow(tableConvertcNomen)>=1){
		NearestStart = tableConvertcNomen$gStart[which(abs(tableConvertcNomen$gStart - gPos) == min(abs(tableConvertcNomen$gStart-  gPos)))]
		NearestEnd = tableConvertcNomen$gEnd[which(abs(tableConvertcNomen$gEnd - gPos) == min(abs(tableConvertcNomen$gEnd - gPos)))]
		NearestStart = NearestStart[1]
		NearestEnd = NearestEnd[1]

		if(abs(NearestStart-gPos)<abs(NearestEnd-gPos)){
			NearestEnd = tableConvertcNomen$gEnd[tableConvertcNomen$gStart==NearestStart]
		}else if(abs(NearestStart-gPos)>abs(NearestEnd-gPos)){
			NearestStart = tableConvertcNomen$gStart[tableConvertcNomen$gEnd==NearestEnd]
		}

		startTranscrit = tableConvertcNomen$gStart[1]
		endTranscrit = tableConvertcNomen$gEnd[nrow(tableConvertcNomen)]
		gCDSstart = tableConvertcNomen$gCDSstart[1]
		gCDSend = tableConvertcNomen$gCDSend[1]

		if(gPos < startTranscrit){
			delta = startTranscrit - gPos
			cPos <<- paste("c.",tableConvertcNomen$cStart[1]-delta,sep="")
		}else if(gPos > endTranscrit){
			delta = gPos - endTranscrit
			cPos <<- paste("c.*",tableConvertcNomen$cEnd[nrow(tableConvertcNomen)]+delta-1,sep="")
		}else if(nrow(tableConvertcNomen[tableConvertcNomen$gStart==gPos,])==1){
			cPos <<- paste("c.",tableConvertcNomen$cStart[tableConvertcNomen$gStart==gPos],sep="")
		}else if(nrow(tableConvertcNomen[tableConvertcNomen$gEnd==gPos,])==1){
			cPos <<- paste("c.",tableConvertcNomen$cEnd[tableConvertcNomen$gEnd==gPos],sep="")
		}else if (gPos < gCDSstart){
			if(gPos > NearestStart & gPos < NearestEnd){
				delta = gPos - NearestStart
				cPos <<- paste("c.",tableConvertcNomen$cStart[tableConvertcNomen$gStart==NearestStart]+delta,sep="")
			}else if(gPos < NearestStart){
				delta = NearestStart - gPos
				cPos <<- paste("c.",tableConvertcNomen$cStart[tableConvertcNomen$gStart==NearestStart],"-",delta,sep="")
			}else if(gPos > NearestEnd){
				delta = gPos - NearestEnd
				cPos <<- paste("c.",tableConvertcNomen$cEnd[tableConvertcNomen$gEnd==NearestEnd],"+",delta,sep="")
			}
		}else if (gPos > gCDSend){
			if(gPos > NearestStart & gPos < NearestEnd){
				delta = NearestEnd - gPos
				cPos <<- paste("c.*",tableConvertcNomen$cEnd[tableConvertcNomen$gEnd==NearestEnd]-delta-1,sep="")
			}else if(gPos < NearestStart){
				delta = NearestStart - gPos
				cPos <<- paste("c.*",tableConvertcNomen$cStart[tableConvertcNomen$gStart==NearestStart]-1,"-",delta,sep="")
			}else if(gPos > NearestEnd){
				delta = gPos - NearestEnd
				cPos <<- paste("c.*",tableConvertcNomen$cEnd[tableConvertcNomen$gEnd==NearestEnd]-1,"+",delta,sep="")
			}
		}else if (gPos > NearestStart & gPos < NearestEnd){
			delta = gPos - NearestStart
			cPos <<- paste("c.",tableConvertcNomen$cStart[tableConvertcNomen$gStart==NearestStart]+delta,sep="")
		}else if (gPos < NearestStart){
			delta = NearestStart - gPos
			cPos <<- paste("c.",tableConvertcNomen$cStart[tableConvertcNomen$gStart==NearestStart],"-",delta,sep="")
		}else if(gPos > NearestEnd){
			delta = gPos - NearestEnd
			cPos <<- paste("c.",tableConvertcNomen$cEnd[tableConvertcNomen$gEnd==NearestEnd],"+",delta,sep="")
		}
	}else{
	cPos <<- 0
	}
}

convertTocNomenRev <- function(gPos){

	if(nrow(tableConvertcNomen)>=1){
		NearestStart = tableConvertcNomen$gStart[which(abs(tableConvertcNomen$gStart - gPos) == min(abs(tableConvertcNomen$gStart-  gPos)))]
		NearestEnd = tableConvertcNomen$gEnd[which(abs(tableConvertcNomen$gEnd - gPos) == min(abs(tableConvertcNomen$gEnd - gPos)))]
		NearestStart = NearestStart[1]
		NearestEnd = NearestEnd[1]

		if(abs(NearestStart-gPos)<abs(NearestEnd-gPos)){
			NearestEnd = tableConvertcNomen$gEnd[tableConvertcNomen$gStart==NearestStart]
		}else if(abs(NearestStart-gPos)>abs(NearestEnd-gPos)){
			NearestStart = tableConvertcNomen$gStart[tableConvertcNomen$gEnd==NearestEnd]
		}

		startTranscrit = tableConvertcNomen$gStart[nrow(tableConvertcNomen)]
		endTranscrit = tableConvertcNomen$gEnd[1]
		gCDSstart = tableConvertcNomen$gCDSstart[1]
		gCDSend = tableConvertcNomen$gCDSend[1]

		if(gPos > startTranscrit){
			delta = gPos - startTranscrit
			cPos <<- paste("c.",tableConvertcNomen$cStart[nrow(tableConvertcNomen)]-delta,sep="")
		}else if(gPos < endTranscrit){
			delta = endTranscrit - gPos
			cPos <<- paste("c.*",tableConvertcNomen$cEnd[1]+delta-1,sep="")
		}else if(nrow(tableConvertcNomen[tableConvertcNomen$gStart==gPos,])==1){
			cPos <<- paste("c.",tableConvertcNomen$cStart[tableConvertcNomen$gStart==gPos],sep="")
		}else if(nrow(tableConvertcNomen[tableConvertcNomen$gEnd==gPos,])==1){
			cPos <<- paste("c.",tableConvertcNomen$cEnd[tableConvertcNomen$gEnd==gPos],sep="")
		}else if (gPos > gCDSstart){
			if(gPos < NearestStart & gPos > NearestEnd){
				delta = NearestStart - gPos
				cPos <<- paste("c.",tableConvertcNomen$cStart[tableConvertcNomen$gStart==NearestStart]+delta,sep="")
			}else if(gPos > NearestStart){
				delta = gPos - NearestStart
				cPos <<- paste("c.",tableConvertcNomen$cStart[tableConvertcNomen$gStart==NearestStart],"-",delta,sep="")
			}else if(gPos < NearestEnd){
				delta = NearestEnd - gPos
				cPos <<- paste("c.",tableConvertcNomen$cEnd[tableConvertcNomen$gEnd==NearestEnd],"+",delta,sep="")
			}
		}else if (gPos < gCDSend){
			if(gPos < NearestStart & gPos > NearestEnd){
				delta = gPos - NearestEnd
				cPos <<- paste("c.*",tableConvertcNomen$cEnd[tableConvertcNomen$gEnd==NearestEnd]-delta-1,sep="")
			}else if(gPos > NearestStart){
				delta = gPos - NearestStart
				cPos <<- paste("c.*",tableConvertcNomen$cStart[tableConvertcNomen$gStart==NearestStart]-1,"-",delta,sep="")
			}else if(gPos < NearestEnd){
				delta = NearestEnd - gPos
				cPos <<- paste("c.*",tableConvertcNomen$cEnd[tableConvertcNomen$gEnd==NearestEnd]-1,"+",delta,sep="")
			}
		}else if (gPos < NearestStart & gPos > NearestEnd){
			delta = NearestStart - gPos
			cPos <<- paste("c.",tableConvertcNomen$cStart[tableConvertcNomen$gStart==NearestStart]+delta,sep="")
		}else if (gPos > NearestStart){
			delta = gPos - NearestStart
			cPos <<- paste("c.",tableConvertcNomen$cStart[tableConvertcNomen$gStart==NearestStart],"-",delta,sep="")
		}else if(gPos < NearestEnd){
			delta = NearestEnd - gPos
			cPos <<- paste("c.",tableConvertcNomen$cEnd[tableConvertcNomen$gEnd==NearestEnd],"+",delta,sep="")
		}
	}else{
	cPos <<- 0
	}
}

#Annotation functions

FuncAnnotJuncsFor <- function(Start, End){

	annotjunc = getAnnotJuncs(Start, End)
	constitutive = annotjunc[1]
	calcul = annotjunc[2]
	nameJunc = annotjunc[3]
	cStart = convertTocNomenFor(Start)
	cEnd =convertTocNomenFor(End+1)
	result <<-list(constitutive,calcul,nameJunc,cStart,cEnd)
}

FuncAnnotJuncsRev <- function(Start, End){

	annotjunc = getAnnotJuncs(End, Start)
	constitutive = annotjunc[1]
	calcul = annotjunc[2]
	nameJunc = annotjunc[3]
	cEnd = convertTocNomenRev(Start)
	cStart = convertTocNomenRev(End+1)
	result <<-list(constitutive,calcul,nameJunc,cStart,cEnd)
}

message("######################")
message("#Annot junctions...")
message("######################")

prog=0
nbTrans = length(unique(data_junction$NM))
IdTrans = NULL
constitutive = NULL
calcul = NULL
AnnotJuncs = NULL
cStart = NULL
cEnd = NULL
data_junction$constitutive = ""
data_junction$event_type = ""
data_junction$AnnotJuncs = ""
data_junction$cStart = ""
data_junction$cEnd = ""
pb = txtProgressBar(min = 0, max = nbTrans, initial = prog, char = "=", style = 3)

for (i in unique(data_junction$NM)){

	prog = prog+1
	setTxtProgressBar(pb,prog)

	strand = data_junction$Strand_transcript[data_junction$NM==i][1]
	tableConvertAnnot = tableConvert[tableConvert$transcrit==i,]
	tmpIdTrans = which(data_junction$NM==i)
	IdTrans = c(IdTrans,tmpIdTrans)

	if(strand=="forward"){
		tableConvertcNomen = tableConvertAnnot
		tableConvertcNomen$gStart = tableConvertcNomen$gStart + 1

		tmpAnnot <- mapply(FuncAnnotJuncsFor,data_junction$start[tmpIdTrans],data_junction$end[tmpIdTrans])

	}else if(strand=="reverse"){
		tableConvertcNomen = tableConvertAnnot
		tableConvertcNomen$gEnd= tableConvertcNomen$gEnd + 1

		tmpAnnot <- mapply(FuncAnnotJuncsRev,data_junction$start[tmpIdTrans],data_junction$end[tmpIdTrans])

	}

	constitutive = c(constitutive,unlist(tmpAnnot[1,]))
	calcul = c(calcul,unlist(tmpAnnot[2,]))
	AnnotJuncs = c(AnnotJuncs,unlist(tmpAnnot[3,]))
	cStart = c(cStart,unlist(tmpAnnot[4,]))
	cEnd = c(cEnd,unlist(tmpAnnot[5,]))

}

data_junction$constitutive[IdTrans] = constitutive
data_junction$event_type[IdTrans] = calcul
data_junction$AnnotJuncs[IdTrans] = AnnotJuncs
data_junction$cStart[IdTrans] = cStart
data_junction$cEnd[IdTrans] = cEnd

#######################################
#perrcentage calculation of junctions #
#######################################
message("\n######################")
message("#Calcul Junction...")
message("######################")


getJunctionPhysioDonFor <- function(PosCryptic){

	tableTranscrit=tableTranscrit[order(tableTranscrit$gStart),]

	if(PosCryptic<=tableTranscrit$gEnd[1]){
		JunctPhysio = paste(tableTranscrit$Chr[1],tableTranscrit$gEnd[1],tableTranscrit$gStart[2],sep="_")
	}else if (PosCryptic>tableTranscrit$gEnd[(nrow(tableTranscrit)-1)]){
		JunctPhysio = paste(tableTranscrit$Chr[1],tableTranscrit$gEnd[(nrow(tableTranscrit)-1)],tableTranscrit$gStart[nrow(tableTranscrit)],sep="_")
	}else{
		IdSupp = which(tableTranscrit$gStart>PosCryptic)[1]
		IdInf = which(tableTranscrit$gStart<PosCryptic)[length(which(tableTranscrit$gStart<PosCryptic))]
		JunctPhysio = paste(tableTranscrit$Chr[1],tableTranscrit$gEnd[IdInf],tableTranscrit$gStart[IdSupp],sep="_")
	}
	return(JunctPhysio)
}

getJunctionPhysioDonRev <- function(PosCryptic){

	tableTranscrit=tableTranscrit[order(tableTranscrit$gStart,decreasing=TRUE),]

	if(PosCryptic>=tableTranscrit$gEnd[1]){
		JunctPhysio = paste(tableTranscrit$Chr[1],tableTranscrit$gStart[2],tableTranscrit$gEnd[1],sep="_")
	}else if (PosCryptic<=tableTranscrit$gEnd[(nrow(tableTranscrit)-1)]){
		JunctPhysio = paste(tableTranscrit$Chr[1],tableTranscrit$gStart[nrow(tableTranscrit)],tableTranscrit$gEnd[(nrow(tableTranscrit)-1)],sep="_")
	}else{
		IdSupp = which(tableTranscrit$gStart<PosCryptic)[1]
		IdInf = which(tableTranscrit$gStart>PosCryptic)[length(which(tableTranscrit$gStart>PosCryptic))]
		JunctPhysio = paste(tableTranscrit$Chr[1],tableTranscrit$gStart[IdSupp],tableTranscrit$gEnd[IdInf],sep="_")
	}
	return(JunctPhysio)
}

getJunctionPhysioAccFor <- function(PosCryptic){

	tableTranscrit=tableTranscrit[order(tableTranscrit$gStart),]

	if(PosCryptic<=tableTranscrit$gStart[2]){
		JunctPhysio = paste(tableTranscrit$Chr[1],tableTranscrit$gEnd[1],tableTranscrit$gStart[2],sep="_")
	}else if (PosCryptic>=tableTranscrit$gStart[nrow(tableTranscrit)]){
		JunctPhysio = paste(tableTranscrit$Chr[1],tableTranscrit$gEnd[(nrow(tableTranscrit)-1)],tableTranscrit$gStart[nrow(tableTranscrit)],sep="_")
	}else{
		IdSupp = which(tableTranscrit$gEnd>PosCryptic)[1]
		IdInf = which(tableTranscrit$gEnd<PosCryptic)[length(which(tableTranscrit$gEnd<PosCryptic))]
		JunctPhysio = paste(tableTranscrit$Chr[1],tableTranscrit$gEnd[IdInf],tableTranscrit$gStart[IdSupp],sep="_")
	}
	return(JunctPhysio)
}

getJunctionPhysioAccRev <- function(PosCryptic){

	tableTranscrit=tableTranscrit[order(tableTranscrit$gStart,decreasing=TRUE),]

		if(PosCryptic>=tableTranscrit$gStart[2]){
			JunctPhysio = paste(tableTranscrit$Chr[1],tableTranscrit$gStart[2],tableTranscrit$gEnd[1],sep="_")
		}else if (PosCryptic<=tableTranscrit$gStart[nrow(tableTranscrit)]){
			JunctPhysio = paste(tableTranscrit$Chr[1],tableTranscrit$gStart[nrow(tableTranscrit)],tableTranscrit$gEnd[(nrow(tableTranscrit)-1)],sep="_")
		}else{
			IdSupp = which(tableTranscrit$gEnd<PosCryptic)[1]
			IdInf = which(tableTranscrit$gEnd>PosCryptic)[length(which(tableTranscrit$gEnd>PosCryptic))]
			JunctPhysio = paste(tableTranscrit$Chr[1],tableTranscrit$gStart[IdSupp],tableTranscrit$gEnd[IdInf],sep="_")
		}
	return(JunctPhysio)
}


calcJuncTo1sampleFor <- function(Conca,calcul, start, end){
	P_junc=rep(0,n_ech)

	if (calcul=="5AS"){
		JunctPhysio <- getJunctionPhysioDonFor(start)

		read_physio_don = mean(dataTmpPhysio[JunctPhysio,input])

		if (length(read_physio_don)==0) {
			calcul="No Ref physio"
		}else{
			read_alt_don=as.numeric(dataTmp[Conca,input])
			P_junc =as.numeric((read_alt_don/read_physio_don)*100)
		}
	}else if (calcul=="3AS"){
		JunctPhysio <- getJunctionPhysioAccFor(end)

		read_physio_acc = mean(dataTmpPhysio[JunctPhysio,input])

		if (length(read_physio_acc)==0) {
			calcul="No Ref physio"
		}else{
			read_alt_acc=as.numeric(dataTmp[Conca,input])
			P_junc =as.numeric((read_alt_acc/read_physio_acc)*100)
		}
	}else if (calcul=="SkipEx"){

		start_saut=start
		end_saut=end
		read_physio_saut=dataTmpPhysio[dataTmpPhysio$start>=start_saut & dataTmpPhysio$end<=end_saut,input]

		if (length(read_physio_saut)==0) {
			calcul="No Ref physio"
		}else{
			read_physio_saut=mean(as.numeric(read_physio_saut))
			read_saut=as.numeric(na.omit(dataTmp[Conca, input]))
			P_junc=as.numeric((read_saut/read_physio_saut)*100)
		}
	}
	result <<- list(calcul, P_junc)
}

calcJuncTo1sampleRev <- function(Conca,calcul, start, end){
	P_junc=rep(0,n_ech)
	if (calcul=="5AS"){
		JunctPhysio <- getJunctionPhysioDonRev(end)

		read_physio_don = mean(dataTmpPhysio[JunctPhysio ,input])

		if (length(read_physio_don)==0) {
			calcul="No Ref physio"
		}else{
			read_alt_don=as.numeric(dataTmp[Conca,input])
			P_junc =as.numeric((read_alt_don/read_physio_don)*100)
		}
	}else if (calcul=="3AS"){

		JunctPhysio <- getJunctionPhysioAccRev(start)

		read_physio_acc = mean(dataTmpPhysio[JunctPhysio,input])

		if (length(read_physio_acc)==0) {
			calcul="No Ref physio"
		}else{
			read_alt_acc=as.numeric(dataTmp[Conca,input])
			P_junc =as.numeric((read_alt_acc/read_physio_acc)*100)
		}
	}else if (calcul=="SkipEx"){

		start_saut=start
		end_saut=end
		read_physio_saut=dataTmpPhysio[dataTmpPhysio$start>=start_saut & dataTmpPhysio$end<=end_saut,input]

		if (length(read_physio_saut)==0) {
			calcul="No Ref physio"
		}else{
			read_physio_saut=mean(as.numeric(read_physio_saut))
			read_saut=as.numeric(na.omit(dataTmp[Conca, input]))
			P_junc=as.numeric((read_saut/read_physio_saut)*100)
		}
	}
	result <<- list(calcul, P_junc)
}

calcJuncToMultiSampleFor <- function(Conca,calcul, start, end){
	P_junc=rep(0,n_ech)

	if (calcul=="5AS"){

		JunctPhysio <- getJunctionPhysioDonFor(start)
		read_physio_don = colMeans(dataTmpPhysio[JunctPhysio,input])

		if (is.null(read_physio_don)){
			calcul="No Ref physio"
		}else{
			read_alt=as.numeric(dataTmp[Conca,input])
			P_junc =as.numeric((read_alt/read_physio_don)*100)
		}
	}else if (calcul=="3AS"){

		JunctPhysio <- getJunctionPhysioAccFor(end)
		read_physio_acc = colMeans(dataTmpPhysio[JunctPhysio,input])

		if (is.null(read_physio_acc)){
			calcul="No Ref physio"
		}else{
			read_alt_acc=as.numeric(dataTmp[Conca,input])
			P_junc =as.numeric((read_alt_acc/read_physio_acc)*100)
		}
	}else if (calcul=="SkipEx"){

		start_saut=start
		end_saut=end

		ref_physio=dataTmpPhysio[dataTmpPhysio$start>=start_saut & dataTmpPhysio$end<=end_saut,input]
		if (is.null(ref_physio)){
			calcul="No Ref physio"
		}else{
			ref_calcul = as.numeric(colMeans(ref_physio))
			read_saut=as.numeric(na.omit(dataTmp[Conca, input]))
			P_junc=as.numeric((read_saut/ref_calcul)*100)
		}
	}
	result <<- list(calcul, P_junc)
}

calcJuncToMultiSampleRev <- function(Conca,calcul, start, end){
	P_junc=rep(0,n_ech)

	if (calcul=="5AS"){

		JunctPhysio <- getJunctionPhysioDonRev(end)
		read_physio_don = colMeans(dataTmpPhysio[JunctPhysio,input])

		if (is.null(read_physio_don)){
			calcul="No Ref physio"
		}else{
			read_alt=as.numeric(dataTmp[Conca,input])
			P_junc =as.numeric((read_alt/read_physio_don)*100)
		}
	}else if (calcul=="3AS"){

		JunctPhysio <- getJunctionPhysioAccRev(start)
		read_physio_acc = colMeans(dataTmpPhysio[JunctPhysio,input])

		if (is.null(read_physio_acc)){
			calcul="No Ref physio"
		}else{
			read_alt_acc=as.numeric(dataTmp[Conca,input])
			P_junc =as.numeric((read_alt_acc/read_physio_acc)*100)
		}

	}else if (calcul=="SkipEx"){

		start_saut=start
		end_saut=end

		ref_physio=dataTmpPhysio[dataTmpPhysio$start>=start_saut & dataTmpPhysio$end<=end_saut,input]
		if (is.null(ref_physio)){
			calcul="No Ref physio"
		}else{
			ref_calcul = as.numeric(colMeans(ref_physio))
			read_saut=as.numeric(na.omit(dataTmp[Conca, input]))
			P_junc=as.numeric((read_saut/ref_calcul)*100)
		}
	}
	result <<- list(calcul, P_junc)
}

prog=0
nbTrans = length(unique(data_junction$NM))
pb1 = txtProgressBar(min = 0, max = nbTrans, initial = prog, char = "=", style = 3)

if(n_ech==1){
	for( i in unique(data_junction$NM)){
		prog = prog+1
        setTxtProgressBar(pb1,prog)

		tmpIdTrans = which(data_junction$NM==i)
		strand = data_junction$Strand_transcript[tmpIdTrans][1]
		dataTmp = data_junction[tmpIdTrans,]
		row.names(dataTmp) = dataTmp$Conca
		dataTmpPhysio = dataTmp[dataTmp$constitutive=="Physio",]
		tableTranscrit = tableConvert[tableConvert$transcrit==i,]

		if(strand=="forward"){
			tmp = mapply(calcJuncTo1sampleFor,Conca = dataTmp$Conca,calcul = dataTmp$event_type, start = dataTmp$start, end = dataTmp$end)
		}else{
			tmp = mapply(calcJuncTo1sampleRev,Conca = dataTmp$Conca,calcul = dataTmp$event_type, start = dataTmp$start, end = dataTmp$start)
		}
		dataTmp[,"event_type"]  = unlist(tmp[1,])
		dataTmp[,SampleOutput] = unlist(tmp[2,])
		data_junction[tmpIdTrans,c(SampleOutput,"event_type")] = dataTmp[,c(SampleOutput,"event_type")]
	}
	data_junction$mean_percent = data_junction[,output]
	data_junction$read_mean = data_junction[,input]
}else{
	for( i in unique(data_junction$NM)){
		prog = prog+1
        setTxtProgressBar(pb1,prog)

		tmpIdTrans = which(data_junction$NM==i)

		strand = data_junction$Strand_transcript[tmpIdTrans][1]
		dataTmp = data_junction[tmpIdTrans,]
		row.names(dataTmp) = dataTmp$Conca
		dataTmpPhysio = dataTmp[dataTmp$constitutive=="Physio",]
		tableTranscrit = tableConvert[tableConvert$transcrit==i,]

		if(strand=="forward"){
			tmp = mapply(calcJuncToMultiSampleFor,Conca = dataTmp$Conca,calcul = dataTmp$event_type, start = dataTmp$start, end = dataTmp$end)
		}else{
			tmp = mapply(calcJuncToMultiSampleRev,Conca = dataTmp$Conca,calcul = dataTmp$event_type, start = dataTmp$start, end = dataTmp$end)
		}
		dataTmp[,"event_type"]  = unlist(tmp[1,])
		dataTmp[,SampleOutput] = matrix(unlist(tmp[2,]),ncol=n_ech,byrow=TRUE)
		data_junction[tmpIdTrans,c(SampleOutput,"event_type")] = dataTmp[,c(SampleOutput,"event_type")]
		T2bis<-as.numeric(format(Sys.time(), "%s"))
	}
	data_junction$mean_percent = as.numeric(rowMeans(data_junction[,output]))
	data_junction$read_mean = as.numeric(rowMeans(data_junction[,input]))
}

#numerique convertion of outputs

eval(parse(text=paste("data_junction[,",output,"] <- as.numeric(data_junction[,",output,"])",sep="",collapse=" ; ")))

#infinite value are fixed to 100 %

data_junction$mean_percent[data_junction$mean_percent=='Inf']=100

data_junction=data_junction[order(data_junction$Gene),]

#comptage nb times of junctions
nbTimes <- function(x){length(x[x!=0])}
data_junction$nbSamp = apply(data_junction[,input],1,nbTimes)

name_output=paste(run,"_outputSpliceLauncher.xlsx",sep="")

#get BED annot Junctions
getBEDannot<-function(chr,start,end,name,strand,cStart,cEnd,meanP,rowSamples,calcul){
    #get RGB
    meanP[is.na(meanP)]=0
    ind=1:length(chr)
    ind=ind[order(meanP,decreasing = TRUE)]
    sortCalc = calcul[order(meanP,decreasing = TRUE)]
    nJunc=table(calcul)
    intCol = seq(from=240,to=0,by=-10)
    sizeSE = ceiling(as.numeric(nJunc["SkipEx"])/length(intCol))
    sizeSA = ceiling(as.numeric(nJunc["3AS"])/length(intCol))
    sizeSD = ceiling(as.numeric(nJunc["5AS"])/length(intCol))
    rangSE = rep(intCol,each=sizeSE)
    rangSA = rep(intCol,each=sizeSA)
    rangSD = rep(intCol,each=sizeSD)
    rgbSE = paste(max(intCol)-rangSE,max(intCol)-rangSE,255,sep=",")[1:nJunc["SkipEx"]]
    rgbSA = paste(max(intCol)-rangSA,255,max(intCol)-rangSA,sep=",")[1:nJunc["3AS"]]
    rgbSD = paste(255,max(intCol)-rangSD,max(intCol)-rangSD,sep=",")[1:nJunc["5AS"]]
    rgb = rep("0,0,0",length(chr))
    rgb[ind[sortCalc=="SkipEx"]]=rgbSE
    rgb[ind[sortCalc=="3AS"]]=rgbSA
    rgb[ind[sortCalc=="5AS"]]=rgbSD
    #get annot
    annot = apply(rowSamples,1,function(x){paste(paste(SampleOutput,x,sep=": "),collapse="|")})
    #get BED
    name = sub(BlackInvTriangle,"ins_",name)
    name = sub(DeltaSymb,"del_",name)
    bedAnnot = c(paste("track name=",run,
                    " type=bedDetail description=\"Alternative junction from ",run,
                    " RNAseq\" db=hg19 itemRgb=\"On\" visibility=\"full\"",sep=""),
                paste(as.character(chr), start, end,name,meanP,as.character(strand),start,end,rgb,1,
                    paste(cStart,cEnd,sep="_"),annot,sep="\t"))
    return(bedAnnot)
}

if(checkBEDannot=="YES"){
    message("\nGet BED annot...")
    bedAnnot = getBEDannot(data_junction$chr[calcul!="Physio"&calcul!="NoData"],data_junction$start[calcul!="Physio"&calcul!="NoData"],
                            data_junction$end[calcul!="Physio"&calcul!="NoData"],data_junction$AnnotJuncs[calcul!="Physio"&calcul!="NoData"],
                            data_junction$strand[calcul!="Physio"&calcul!="NoData"],data_junction$cStart[calcul!="Physio"&calcul!="NoData"],
                            data_junction$cEnd[calcul!="Physio"&calcul!="NoData"],data_junction$mean_percent[calcul!="Physio"&calcul!="NoData"],
                            data_junction[calcul!="Physio"&calcul!="NoData",SampleOutput],data_junction$event_type[calcul!="Physio"&calcul!="NoData"])
    cat(bedAnnot,file=paste(chem_results,run,".bed",sep=""),sep="\n")
}

checkSizeOutput <- function(data){
	if (nrow(data)>1000000){
		message("***More than 1,000,000 junctions to save")
		message(paste("   remove unique junctions:",nrow(data[data$nbSamp==1,])))
		data = data[data$nbSamp>1,]
	}
	return(data)
}

printInText <- function(data,path,name){
    data$AnnotJuncs = sub(BlackInvTriangle,"ins_",data$AnnotJuncs)
    data$AnnotJuncs = sub(DeltaSymb,"del_",data$AnnotJuncs)
    write.table(data,file=paste(path,name,sep=""),row.names=F,dec=".",sep="\t",quote=FALSE)
}

dataToPrint = subset( data_junction, select = -c(ID_gene, constitutive))
dataToPrint = checkSizeOutput(dataToPrint)

if(StatAnalysis=="NO"){
	message("\n   Data are saving...")
    if(!printText){
        tryCatch({
         WriteXLS(dataToPrint, ExcelFileName =paste(chem_results,name_output,sep=""),row.names=F,SheetNames=run,Encoding = "UTF-8",na="")
        },
         error=function(cond) {
             message("Here's the original error message:")
             message(cond)
             name_output=paste(run,"_outputSpliceLauncher.csv",sep="")
             write.table(dataToPrint,file=paste(chem_results,name_output,sep=""),row.names=F,dec=".",sep=";")
         },
         warning=function(cond) {
             message("Here's the original warning message:")
             message(cond)
             name_output=paste(run,"_outputSpliceLauncher.csv",sep="")
             write.table(dataToPrint,file=paste(chem_results,name_output,sep=""),row.names=F,dec=".",sep=";")
         })
    }else if(printText){
        printInText(dataToPrint,chem_results,paste(run,"_outputSpliceLauncher.txt",sep=""))
    }

}else{
message("\n######################")
message("#Statistical analysis...")
message("######################")

#####################################################
#Code R Valentin to calculate junctions distribution#
#####################################################
setwd(paste(chem_results,"/",sep=""))

a<-data_junction

#Functions used by the Valentin's script#

count<-function(x, values){
  res<-rep(NA,length(values))
  for(i in 1:length(values)) res[i]<-sum(x==values[i])
  return(res)
}
msd<-function(x){
  m<-mean(x,na.rm=T) ; sd<-sd(x,na.rm=T)
  if(m<1e-3){ m<-"<1e-3" }else{ m<-format(m,digits=3,nsmall=3) }
  if(sd<1e-3){ sd<-"<1e-3" }else{ sd<-format(sd,digits=3,nsmall=3) }
  return(paste(m,sd,sep=";"))
}
format.p<-function(p){
  if(p<1e-16) return("<1e-16")
  if(p<1e-6) return("<1e-6")
  if(p<1e-3) return("<1e-3")
  return(format(p,digits=3,nsmall=3))
}

cuter<-function (x, cuts = NULL, g = 2, mod = "[[", minmax = FALSE)
{
  mod <- match.arg(mod, c("[[", "]]"))
  if (is.null(cuts)) {
    cuts <- quantile(x, seq(0, 1, 1/g), na.rm = T)
    cuts <- cuts[-c(1, length(cuts))]
  }
  res <- x
  ncut <- length(cuts)
  levels <- rep(NA, ncut + 1)
  if (mod == "[[") {
    res[x < cuts[1]] <- 1
    res[x >= cuts[ncut]] <- ncut + 1
    levels[1] <- paste("<", cuts[1], sep = "")
    levels[ncut + 1] <- paste(">=", cuts[ncut], sep = "")
    if (minmax) {
      m1 <- min(x)
      m2 <- max(x)
      if (m1 < cuts[1]) {
        levels[1] <- paste("[", format.p(m1), ",", cuts[1],
                           "[", sep = "")
      }
      if (m2 >= cuts[ncut]) {
        levels[ncut + 1] <- paste("[", cuts[ncut], ",",
                                  format.p(m2), "]", sep = "")
      }
    }
    if (ncut > 1) {
      for (i in 1:(ncut - 1)) {
        res[x >= cuts[i] & x < cuts[i + 1]] <- i + 1
        levels[i + 1] <- paste("[", cuts[i], ",", cuts[i +
                                                         1], "[", sep = "")
      }
    }
  }
  else if (mod == "]]") {
    res[x <= cuts[1]] <- 1
    res[x > cuts[ncut]] <- ncut + 1
    levels[1] <- paste("<=", cuts[1], sep = "")
    levels[ncut + 1] <- paste(">", cuts[ncut], sep = "")
    if (minmax) {
      m1 <- min(x)
      m2 <- max(x)
      if (m1 <= cuts[1]) {
        levels[1] <- paste("[", format.p(m1), ",", cuts[1],
                           "]", sep = "")
      }
      if (m2 > cuts[ncut]) {
        levels[ncut + 1] <- paste("]", cuts[ncut], ",",
                                  format.p(m2), "]", sep = "")
      }
    }
    if (ncut > 1) {
      for (i in 1:(ncut - 1)) {
        res[x > cuts[i] & x <= cuts[i + 1]] <- i + 1
        levels[i + 1] <- paste("]", cuts[i], ",", cuts[i +
                                                         1], "]", sep = "")
      }
    }
  }
  res <- factor(res, levels = as.character(1:(ncut + 1)))
  levels(res) <- levels
  return(res)
}

vh.like <-
  function(pattern, x){
    w<-regexpr(pattern,x)
    w<-which(na(w==1) & na(nchar(x)==attr(w,"match.length")))
    res<-x ; res[]<-NA ; class(res)<-"logical" ; res[!is.na(x)]<-FALSE ; res[w]<-TRUE
    return(res)
  }

na<-function(x, to = FALSE){
    x[is.na(x)]<-to
    x
}

#Function to transform a table jid*sid in sid*jid
# a : initial table
# sid : samples ID (names of columns from a)
# jid.column : name of column with junctions ID

js2sj<-function(a, sid, jid.column){

  if(!all(c(jid.column,sid)%in%names(a))){
    stop("***** The 'jid' and 'sid' don't find in the columns of table")
  }else if(length(unique(a[,jid.column]))!=nrow(a)){
    stop("***** The junctions were not identify in unique way")
  }
  a = a[a$constitutive!='Physio',]
  a = a[a$event_type!='NoData',]
  a = a[a$nbSamp>=5,]

  for(i in 1:n_ech){
  a = a[a[,sid[i]]<=1000,]
  }
  data = as.data.frame (t(a[,sid]))
  names(data) = a[,jid.column]
  row.names(data) = sid
  data = data[,apply(data,2,function(data) !any(is.na(data)))]
  return(data)
}

# Function to create Gamma / NegBinomial models from dataset
# data : table of data sid*jid
# jid : Junctions ID to modelisate
# negbinom.n : the min number of classes to create for negative binomial ajustment )


fit.gamma.negbinomial<-function(data, jid, negbinom.n=10){

  model<-list()

  for(j in jid){

    model[[j]]<-list()
    xp<-x<-data[,j]

    x<-x[!is.na(x)]
    bx<-boxplot.stats(x)
    tot<-bx$stats[4]+1.5*(bx$stats[4]-bx$stats[2])
    model[[j]][["outliers.ub"]]<-tot
    nto<-0
    if(length(bx$out>0)){
      nto<-sum(x%in%bx$out)
      x<-x[!x%in%bx$out]
    }
    model[[j]][["outliers.nb"]]<-nto

    # Modelisation
    m<-mean(x, na.rm=T)
    s<-sd(x, na.rm=T)

    model[[j]][["mean"]]<-m
    model[[j]][["sd"]]<-s

    if(bx$stats[4]==0){  # 3rd quartile to 0
      model[[j]][["model"]]<-"inexistant event"
      next
    }

    # Assay gamma distribution
    scale<-s**2/m
    shape<-m/scale
	tryCatch({
		test<-ks.test(x, pgamma, shape=shape, scale=scale)$p.value
	},
	error=function(cond) {
		message("Here's the original error message:")
		message(cond)
		test = 0
	})
	if(is.na(test)) test=0
    if(test>=0.01){ # acceptable adequation to the gamma distribution
      model[[j]][["model"]]<-"Gamma"
      model[[j]][["model.p"]]<-test
      model[[j]][["scale"]]<-scale
      model[[j]][["shape"]]<-shape
      next
    }

    # Assay negative binomiale distribution with intervalls of negbinom.n of same size
    step<-max(x)/negbinom.n
    y<-as.numeric(cuter(x, cuts=c(seq(step, max(x), step)), mod="]]"))-1

    my<-mean(y, na.rm=T)
    vy<-var(y, na.rm=T)
    prob<-my/vy  ;  size<-my*prob/(1-prob)

    obs<-count(y, 0:max(y))
    th<-dnbinom(0:(max(y)-1), size=size, prob=prob) ; th<-c(th,1-sum(th))
	tryCatch({
		test<-chisq.test(obs, p=th)$p.value
	},
	error=function(cond) {
		message("Here's the original error message:")
		message(cond)
		test = 0
	})
	if(is.na(test)) test=0
    if(test>=0.01){ # acceptable adequation to the negative binomiale distribution
      model[[j]][["model"]]<-paste("Negative binomial ",negbinom.n,"c",sep="")
      model[[j]][["model.p"]]<-test
      model[[j]][["step"]]<-step
      model[[j]][["prob"]]<-prob
      model[[j]][["size"]]<-size
      next
    }

    model[[j]][["model"]]<-"Aucun"
  }

  return(model)
}

#Function to apply the previously model on a serie of values
# model : is the model got by 'fit.gamma.negbinomial'
# x : a serie of values
# jid : junction ID designating the modele to apply

predict.gamma.negbinomial<-function(model, x, jid){

  m<-model[[jid]]

  if(m$model=="Gamma"){
    val<-pgamma(x, shape=m$shape, scale=m$scale, lower.tail=F)
    val[na(val==0)]<-1e-324
  }else if(substr(m$model,1,17)=="Negative binomial"){
    y<-as.numeric(cuter(x, cuts=c(seq(m$step, max(x,na.rm=T), m$step)), mod="]]"))-1
    val<-pnbinom(y, size=m$size, prob=m$prob, lower.tail=F)
    val[na(val==0)]<-1e-324
  }else{
    val<-rep(NA,length(x))
  }
  return(val)
}

###########################
#Aplication on the dataset#
###########################

#########################
# 2 - Fromating to an input ind*evenement : data
message("   Creating matrix...")
# Create table of data analysis
data<-js2sj(a, sid, jid.column)
jid<-names(data)[-1]

#########################
# 3 - Modelization : model
message("   Fitting model...")

model<-fit.gamma.negbinomial(data, jid, negbinom.n=negbinom.n)

#########################
# 4 - Application of model on the dataset
message("   Model apliccation...")

newdata<-data # forme sid*jid

# Count of the number of tests carried out
nb.tests<-0
for(j in names(model)){
    if(j%in%names(newdata) & vh.like("(Gamma)|(Negative binomial.*)",model[[j]][["model"]])){
        nb.tests<-nb.tests+sum(!is.na(newdata[,j]))
    }
}

message("   tmp file with p-values...")
output1<-"output with adjustments.csv"
write(paste("Jonction;Nb Tuckey's outliers;Ajustement;Moyenne;Ecart-type;p-value Goodness of fit;",paste(EchName,collapse=";"), sep=""), output1)

err = NULL
for(j in names(model)){

  m<-model[[j]]

  if(j%in%names(newdata) & vh.like("(Gamma)|(Negative binomial.*)",m$model)){

    val<-predict.gamma.negbinomial(model, x=newdata[,j], jid=j)
    val<-p.adjust(val, "holm", n=nb.tests)
    write(paste(j, m$outliers.nb, m$model, m$mean, m$sd, m$model.p, paste(val,collapse=";"),sep=";"), output1, append=T)

  }
  # Evenements non modelises
  jid.other<-setdiff(names(data), c("sid",names(model)))
  if(length(jid.other)){
	err = c(err,j)
  }
}
message(paste("   Following modeled junctions are not in the models:",(length(err)/nrow(data_junction))*100,'%'))

data_junction_pvalue = read.table("output with adjustments.csv",header=T,dec=".",sep=";")
data_junction_pvalue = data_junction_pvalue[order(data_junction_pvalue$Jonction),]
data_junction = data_junction[order(data_junction$Conca),]
getSign = function(v,n){
    id = which(v<0.05)
    if(length(id)==0){sign="No"}else{sign = paste("Yes",paste(n[id],v[id],collapse="; ",sep=", p-value = "),sep=": ")}
    return(sign)
}

data_junction$DistribAjust = NA
data_junction$DistribAjust[data_junction$Conca%in%data_junction_pvalue$Jonction] = data_junction_pvalue$Ajustement
significative = apply(data_junction_pvalue[,7:ncol(data_junction_pvalue)],1,getSign,names(data_junction)[input])
data_junction$Significative[data_junction$Conca%in%data_junction_pvalue$Jonction] = significative

dataToPrint = subset( data_junction, select = -c(ID_gene, constitutive))
dataToPrint=dataToPrint[order(dataToPrint$Gene),]
dataToPrint = checkSizeOutput(dataToPrint)
message("   Data are saving...")

file.remove("output with adjustments.csv")
file.remove("output without adjustment.csv")

if(!printText){
    tryCatch({
     WriteXLS(dataToPrint, ExcelFileName =paste(chem_results,name_output,sep=""),row.names=F,SheetNames=run,Encoding = "UTF-8",na="")
    },
     error=function(cond) {
         message("Here's the original error message:")
         message(cond)
         name_output=paste(run,"_outputSpliceLauncher.csv",sep="")
         write.table(dataToPrint,file=paste(chem_results,name_output,sep=""),row.names=F,dec=".",sep=";")
     },
     warning=function(cond) {
         message("Here's the original warning message:")
         message(cond)
         name_output=paste(run,"_outputSpliceLauncher.csv",sep="")
         write.table(dataToPrint,file=paste(chem_results,name_output,sep=""),row.names=F,dec=".",sep=";")
     })
}else if(printText){
    printInText(dataToPrint,chem_results,paste(run,"_outputSpliceLauncher.txt",sep=""))
}
}

if(DisplayGraph=="YES"){

	chem_dessin=paste(chem_results ,run,"_figures_output/",sep="")
	dir.create(chem_dessin,showWarnings=FALSE)
	#########################
	#graphic representation #
	#########################
	message("######################")
	message("#Graphic representation...")
	message("######################")

	getSigniDon <- function(data){
		if(StatAnalysis=="YES"){
			data_pvalue = data_junction_pvalue[,c(1,6+j)]
			data_pvalue[,1] = as.character(data_pvalue[,1])
			data= data[order(data$Conca),]

			Junc_don=as.character(na.omit(data$Conca))
			Signi_don = rep(' ',length(Junc_don))
			Pvalue_don = rep(1,length(Junc_don))
			Pvalue_don[Junc_don%in%data_pvalue$Jonction] = data_pvalue[data_pvalue$Jonction%in%Junc_don,2]
			Signi_don[Pvalue_don<0.05& Pvalue_don>0.01] = '*'
			Signi_don[Pvalue_don<0.01& Pvalue_don>0.001] = '**'
			Signi_don[Pvalue_don<0.001] = '***'

		}else{
			Junc_don=as.character(na.omit(data$Conca))
			Signi_don = rep(' ',length(Junc_don))
		}
		return(Signi_don)
	}

	getSigniAcc <- function(data){
		if(StatAnalysis=="YES"){
			data_pvalue = data_junction_pvalue[,c(1,6+j)]
			data_pvalue[,1] = as.character(data_pvalue[,1])
			data= data[order(data$Conca),]

			Junc_acc=as.character(na.omit(data$Conca))
			Signi_acc = rep(' ',length(Annot_acc))
			Pvalue_acc = rep(1,length(Junc_acc))
			Pvalue_acc[Junc_acc%in%data_pvalue$Jonction] = data_pvalue[data_pvalue$Jonction%in%Junc_acc,2]
			Signi_acc[Pvalue_acc<0.05& Pvalue_acc>0.01] = '*'
			Signi_acc[Pvalue_acc<0.01& Pvalue_acc>0.001] = '**'
			Signi_acc[Pvalue_acc<0.001] = '***'
		}else{
			Junc_acc=as.character(na.omit(data$Conca))
			Signi_acc = rep(' ',length(Junc_acc))
		}
		return(Signi_acc)
	}

	getSigniSaut <- function(data){
		if(StatAnalysis=="YES"){
			data_pvalue = data_junction_pvalue[,c(1,6+j)]
			data_pvalue[,1] = as.character(data_pvalue[,1])
			data= data[order(data$Conca),]

			Junc_saut=as.character(na.omit(data$Conca))
			Signi_saut = rep(' ',length(Junc_saut))
			Pvalue_saut = rep(1,length(Junc_saut))
			Pvalue_saut[Junc_saut%in%data_pvalue$Jonction] = data_pvalue[data_pvalue$Jonction%in%Junc_saut,2]
			Signi_saut[Pvalue_saut<0.05& Pvalue_saut>0.01] = '*'
			Signi_saut[Pvalue_saut<0.01& Pvalue_saut>0.001] = '**'
			Signi_saut[Pvalue_saut<0.001] = '***'
		}else{
			Junc_saut=as.character(na.omit(data$Conca))
			Signi_saut = rep(' ',length(Junc_saut))
		}
		return(Signi_saut)
	}

	getFont <- function(Signi){
		Font = rep(1,length(Signi))
		Font[Signi!=" "]=4
		return(Font)
	}

    maxReadCount = apply(data_junction[data_junction$constitutive=="Physio",input],1,max,na.rm=TRUE)
    plottedGene = as.character(data_junction$Gene[data_junction$constitutive=="Physio"])
    plottedGene = unique(plottedGene[maxReadCount>100])
    countPhysio = table(data_junction$Gene[data_junction$constitutive=="Physio"])
    plottedGene = c(plottedGene,names(countPhysio)[as.numeric(countPhysio)>2])
    plottedGene = plottedGene[duplicated(plottedGene)]
    message(paste("  ",length(plottedGene),"plotted gene(s)"))
    data_junction = data_junction[which(data_junction$Gene %in% plottedGene),]
	data_junction=data_junction[order(data_junction$start),]
	lty='dotted'
    id_gene = unique(data_junction$ID_gene)
    id_gene = id_gene[order(id_gene)]

	for(j in 1:length(EchName)){
		message(paste("   Graphics for:",EchName[j]))

		CairoPDF(file=paste(chem_dessin,run,"_",EchName[j],".pdf",sep=""),width =60, height = 12)

		data_junction$mean_percent = data_junction[,SampleOutput[j]]
		data_junction$mean_percent[data_junction$mean_percent=='Inf']=100
        GeneMis = 0
		for(i in id_gene){
			data_gene = data_junction[data_junction$ID_gene==i & data_junction$event_type!="NoData",]

			if(max(data_gene[data_gene$constitutive=="Physio",input[j]])<100){
				GeneMis = GeneMis+1
            }else{

                if(data_gene$Strand_transcript[1]=="forward"){

					start_physio=data_gene$start[data_gene$constitutive=="Physio"]
					lab_physio=c(1:(length(start_physio)+1))
					end_physio=data_gene$end[data_gene$constitutive=="Physio"]
					mean_read_physio=data_gene[data_gene$constitutive=="Physio",SampleInput[j]]
					h_physio=log(mean_read_physio)/20

					start_don=as.numeric(na.omit(data_gene$start[data_gene$event_type=="5AS" & data_gene$mean_percent>=thr]))
					end_don=as.numeric(na.omit(data_gene$end[data_gene$event_type=="5AS" & data_gene$mean_percent>=thr]))
					P_don=as.numeric(na.omit(data_gene$mean_percent[data_gene$event_type=="5AS" & data_gene$mean_percent>=thr]))
					Annot_don=as.character(na.omit(data_gene$AnnotJunc[data_gene$event_type=="5AS" & data_gene$mean_percent>=thr]))
					Signi_don = getSigniDon(data_gene[data_gene$event_type=="5AS" & data_gene$mean_percent>=thr,])
					font_Don = getFont(Signi_don)

					start_acc=as.numeric(na.omit(data_gene$start[data_gene$event_type=="3AS" & data_gene$mean_percent>=thr]))
					end_acc=as.numeric(na.omit(data_gene$end[data_gene$event_type=="3AS" & data_gene$mean_percent>=thr]))
					P_acc=as.numeric(na.omit(data_gene$mean_percent[data_gene$event_type=="3AS" & data_gene$mean_percent>=thr]))
					Annot_acc=as.character(na.omit(data_gene$AnnotJunc[data_gene$event_type=="3AS" & data_gene$mean_percent>=thr]))
					Signi_acc = getSigniAcc(data_gene[data_gene$event_type=="3AS" & data_gene$mean_percent>=thr,])
					font_Acc = getFont(Signi_acc)

					start_saut=as.numeric(na.omit(data_gene$start[data_gene$event_type=="SkipEx" & data_gene$mean_percent>=thr]))
					P_saut=as.numeric(na.omit(data_gene$mean_percent[data_gene$event_type=="SkipEx" & data_gene$mean_percent>=thr]))
					Annot_saut=as.character(na.omit(data_gene$AnnotJunc[data_gene$event_type=="SkipEx" & data_gene$mean_percent>=thr]))
					end_saut=as.numeric(na.omit(data_gene$end[data_gene$event_type=="SkipEx" & data_gene$mean_percent>=thr]))
					Signi_saut = getSigniSaut(data_gene[data_gene$event_type=="SkipEx" & data_gene$mean_percent>=thr,])
					font_Saut = getFont(Signi_saut)

					min_x=min(c(start_don,start_physio,start_saut))
					max_x=max(c(end_don,end_physio,end_saut,end_acc))

					plot(x=0,xlim=c(min_x,max_x),ylim=c(0.8,2.5),xlab=data_gene$Strand_transcript[1],
						main=paste(data_gene$Gene[1],paste("(",unique(data_gene$NM),")",sep=""),
						">",thr,"%",run),cex.main=4,cex.lab=3, yaxt="n", ylab="")
					segments(x0=min(start_physio),y0=1.3,y1=1.3,x1=max(end_physio),col='black',lwd=6)
					rect(xleft=start_physio[-1],xright=end_physio[-length(end_physio)],ytop=rep(1.40,length(start_physio[-1])),
						ybottom=rep(1.20,length(end_physio)),col='blue')
					rect(xleft=min(c(start_don,start_physio,start_saut)),xright=start_physio[1],ytop=1.45,ybottom=1.15,col='blue')
					abline(h=seq(from=2.1,to=2.3,by=0.2),col='gray85',lty=2,lwd=1)
					abline(h=seq(from=2,to=2.5,by=0.2),col='gray70',lty=2,lwd=1)
					rect(xleft=start_physio + (end_physio - start_physio)/2 -((max_x-min_x)/800),
						xright=start_physio + (end_physio - start_physio)/2+((max_x-min_x)/800),
						ytop=h_physio+1.7,ybottom=2,col='blue')
					rect(xleft=end_physio[length(end_physio)],xright=max(c(end_physio,end_acc,end_saut)),ytop=1.45,ybottom=1.15,col='blue')
					segments(x0=start_don,y0=rep(1.25,length(start_don)),y1=rep(1.35,length(start_don)),x1=start_don,col='red',lwd=3)
					segments(x0=end_acc,y0=rep(1.25,length(end_acc)),y1=rep(1.35,length(end_acc)),x1=end_acc,col='green',lwd=3)
					segments(x0=start_don,y0=rep(1.3,length(start_don)),y1=rep(1.0,length(start_don)),x1=start_don+(end_don-start_don)/2,col='red',
						lty=lty,lwd=3)
					segments(x0=start_don+(end_don-start_don)/2,y0=rep(1.0,length(start_don)),y1=rep(1.3,length(start_don)),x1=end_don,col='red',
						lty=lty,lwd=3)
					segments(x0=start_acc,y0=rep(1.3,length(start_acc)),y1=rep(1.0,length(start_acc)),x1=start_acc+(end_acc-start_acc)/2,col='green',
						lty=lty,lwd=3)
					segments(x0=start_acc+(end_acc-start_acc)/2,y0=rep(1.0,length(start_acc)),y1=rep(1.3,length(start_acc)),x1=end_acc,col='green',
						lty=lty,lwd=3)
					segments(x0=start_saut,y0=rep(1.4,length(start_saut)),y1=rep(1.6,length(start_saut)),x1=start_saut+(end_saut-start_saut)/2,col='blue',
						lty=lty,lwd=3)
					segments(x0=start_saut+(end_saut-start_saut)/2,y0=rep(1.6,length(start_saut)),y1=rep(1.4,length(start_saut)),x1=end_saut,col='blue',
						lty=lty,lwd=3)
					text(y=rep(1.1,length(c(start_physio,max(end_physio)))),cex=2,x=c(min(start_physio),
						start_physio[-1]+( end_physio[-length(end_physio)]- start_physio[-1])/2,max(end_physio)),
						labels=c(lab_physio))
					text(y=2.5,x=min_x+(max_x-min_x)*0.05,labels="Physiological expression:",bg='white',cex=2.5,font=4)
					text(y=1.9,x=min_x+(max_x-min_x)*0.05,labels="Alternative splicing:",bg='white',cex=2.5,font=4)
					text(y=2.42,x=min_x+(max_x-min_x)*0.05,labels=paste("Global expression =",round(mean(mean_read_physio),0)),bg='white',cex=2)
					text(y=h_physio+1.75,x=start_physio + (end_physio - start_physio)/2 -((max_x-min_x)/800),
						labels=round(mean_read_physio,0),bg='white',cex=1)
					text(y=c(1.7,1.8),x=start_saut+(end_saut-start_saut)/2,labels=paste(Annot_saut,":",round(P_saut,2),"%",Signi_saut),
						bg='white',cex=1.5,col='blue',font = font_Saut)
					text(y=c(0.95,0.9),x=start_don+(end_don-start_don)/2,labels=paste(Annot_don,":",round(P_don,2),"%",Signi_don),
						bg='white',cex=1.5,col='red',font = font_Don)
					text(y=c(0.85,0.8),x=start_acc+(end_acc-start_acc)/2,labels=paste(Annot_acc,":",round(P_acc,2),"%",Signi_acc),
						bg='white',cex=1.5,col='green',font = font_Acc)

				}else if (data_gene$Strand_transcript[1]=="reverse"){

					data_gene=data_gene[order(data_gene$start),]
					end_physio=data_gene$end[data_gene$constitutive=="Physio"]
					start_physio=data_gene$start[data_gene$constitutive=="Physio"]
					lab_physio=c((length(end_physio)+1):1)
					mean_read_physio=data_gene[data_gene$constitutive=="Physio",SampleInput[j]]
					h_physio=log(mean_read_physio)/20

					start_don=as.numeric(na.omit(data_gene$end[data_gene$event_type=="5AS" & data_gene$mean_percent>=thr]))
					end_don=as.numeric(na.omit(data_gene$start[data_gene$event_type=="5AS" & data_gene$mean_percent>=thr]))
					P_don=as.numeric(na.omit(data_gene$mean_percent[data_gene$event_type=="5AS" & data_gene$mean_percent>=thr]))
					Annot_don=as.character(na.omit(data_gene$AnnotJunc[data_gene$event_type=="5AS" & data_gene$mean_percent>=thr]))
					Signi_don = getSigniDon(data_gene[data_gene$event_type=="5AS" & data_gene$mean_percent>=thr,])
					font_Don = getFont(Signi_don)

					start_acc=as.numeric(na.omit(data_gene$end[data_gene$event_type=="3AS" & data_gene$mean_percent>=thr]))
					end_acc=as.numeric(na.omit(data_gene$start[data_gene$event_type=="3AS" & data_gene$mean_percent>=thr]))
					P_acc=as.numeric(na.omit(data_gene$mean_percent[data_gene$event_type=="3AS" & data_gene$mean_percent>=thr]))
					Annot_acc=as.character(na.omit(data_gene$AnnotJunc[data_gene$event_type=="3AS" & data_gene$mean_percent>=thr]))
					Signi_acc = getSigniAcc(data_gene[data_gene$event_type=="3AS" & data_gene$mean_percent>=thr,])
					font_Acc = getFont(Signi_acc)

					start_saut=as.numeric(na.omit(data_gene$end[data_gene$event_type=="SkipEx" & data_gene$mean_percent>=thr]))
					end_saut=as.numeric(na.omit(data_gene$start[data_gene$event_type=="SkipEx" & data_gene$mean_percent>=thr]))
					P_saut=as.numeric(na.omit(data_gene$mean_percent[data_gene$event_type=="SkipEx" & data_gene$mean_percent>=thr]))
					Annot_saut=as.character(na.omit(data_gene$AnnotJunc[data_gene$event_type=="SkipEx" & data_gene$mean_percent>=thr]))
					Signi_saut = getSigniSaut(data_gene[data_gene$event_type=="SkipEx" & data_gene$mean_percent>=thr,])
					font_Saut = getFont(Signi_saut)

					min_x=min(c(end_physio,end_acc,end_saut,start_don,start_physio,start_saut))
					max_x=max(c(start_don,start_physio,start_saut,end_physio,end_acc,end_saut))

					plot(x=0,xlim=c(min_x,max_x),ylim=c(0.8,2.5),cex.main=4,cex.lab=3,xlab=data_gene$Strand_transcript[1],
						main=paste(data_gene$Gene[1],paste("(",unique(data_gene$NM),")",sep=""),
						">",thr,"%",run), yaxt="n", ylab="")
					segments(x0=min(start_physio),y0=1.3,y1=1.3,x1=max(end_physio),col='black',lwd=6)
					rect(xleft=start_physio[-1],xright=end_physio[-length(end_physio)],ytop=rep(1.40,length(start_physio[-1])),
						ybottom=rep(1.20,length(end_physio)),col='blue')
					abline(h=seq(from=2.1,to=2.3,by=0.2),col='gray85',lty=2,lwd=1)
					abline(h=seq(from=2,to=2.4,by=0.2),col='gray70',lty=2,lwd=1)
					rect(xleft=end_physio + (start_physio-end_physio)/2+((max_x-min_x)/800),
						xright=end_physio + (start_physio-end_physio)/2-((max_x-min_x)/800),
						ytop=h_physio+1.7,ybottom=2,col='blue')
					rect(xleft=min(c(start_don,start_physio,start_saut)),xright=start_physio[1],ytop=1.45,ybottom=1.15,col='blue')
					rect(xleft=end_physio[length(end_physio)],xright=max(c(end_physio,end_acc,end_saut)),ytop=1.45,ybottom=1.15,col='blue')
					segments(x0=start_don,y0=rep(1.25,length(start_don)),y1=rep(1.35,length(start_don)),x1=start_don,col='red',lwd=3)
					segments(x0=end_acc,y0=rep(1.25,length(end_acc)),y1=rep(1.35,length(end_acc)),x1=end_acc,col='green',lwd=3)
					segments(x0=start_saut,y0=rep(1.4,length(start_saut)),y1=rep(1.6,length(start_saut)),x1=start_saut+(end_saut-start_saut)/2,col='blue',
						lty=lty,lwd=3)
					segments(x0=start_saut+(end_saut-start_saut)/2,y0=rep(1.6,length(start_saut)),y1=rep(1.4,length(start_saut)),x1=end_saut,col='blue',
						lty=lty,lwd=3)
					segments(x0=start_don,y0=rep(1.3,length(start_don)),y1=rep(1,length(start_don)),x1=start_don+(end_don-start_don)/2,col='red',
						lty=lty,lwd=3)
					segments(x0=start_don+(end_don-start_don)/2,y0=rep(1,length(start_don)),y1=rep(1.3,length(start_don)),x1=end_don,col='red',
						lty=lty,lwd=3)
					segments(x0=start_acc,y0=rep(1.3,length(start_acc)),y1=rep(1,length(start_acc)),x1=start_acc+(end_acc-start_acc)/2,col='green',
						lty=lty,lwd=3)
					segments(x0=start_acc+(end_acc-start_acc)/2,y0=rep(1,length(start_acc)),y1=rep(1.3,length(start_acc)),x1=end_acc,col='green',
						lty=lty,lwd=3)
					text(y=rep(1.1,length(c(end_physio,max(start_physio)))),cex=2,x=c(min(start_physio),
						end_physio[-length(end_physio)]+( start_physio[-1]- end_physio[-length(end_physio)])/2,max(end_physio)),
						labels=c(lab_physio))
					text(y=2.5,x=min_x+(max_x-min_x)*0.05,labels="Physiological expression:",bg='white',cex=2.5,font=4)
					text(y=1.9,x=min_x+(max_x-min_x)*0.05,labels="Alternative splicing:",bg='white',cex=2.5,font=4)
					text(y=2.42,x=min_x+(max_x-min_x)*0.05,labels=paste("Global expression =",round(mean(mean_read_physio),0)),bg='white',cex=2)
					text(y=h_physio+1.75,x=end_physio + (start_physio-end_physio)/2+((max_x-min_x)/800),
						labels=round(mean_read_physio,0),bg='white',cex=1)
					text(y=c(1.7,1.8),x=end_saut+(start_saut-end_saut)/2,labels=paste(Annot_saut,":",round(P_saut,2),"%",Signi_saut),
						bg='white',cex=1.5,col='blue',font = font_Saut)
					text(y=c(0.95,0.9),x=end_don+(start_don-end_don)/2,labels=paste(Annot_don,":",round(P_don,2),"%",Signi_don),
						bg='white',cex=1.5,col='red',font = font_Don)
					text(y=c(0.85,0.8),x=end_acc+(start_acc-end_acc)/2,labels=paste(Annot_acc,":",round(P_acc,2),"%",Signi_acc),
						bg='white',cex=1.5,col='green',font = font_Acc)
				}
			}
		}
        message(paste("    Not enough read on physiological junctions for",if(GeneMis!=0){GeneMis}else{"NONE"},"gene(s)"))
		dev.off()
	}
}

#########################
#edit report of analysis#
#########################
T2<-Sys.time()

Tdiff= difftime(T2, T1, units ="mins")
if(DisplayGraph=="NO"){
chem_dessin = ""
}

chapitre=c("date of analysis","name of input file", "run id", "Nb of samples", rep("sample",n_ech),
"path to output", "name of output file", "path to graphics output",
"Junctions out of strand","average depth")
name_input=inputFile
donnee=c(paste(date_analyse,"; Time to analysis:",Tdiff,"min"),name_input,run,n_ech,names(data_junction)[input],
			outputDir,name_output,chem_dessin,nbJuncOpp,mean(data_junction$read_mean))

rapport=data.frame(chapitre,donnee)
write.table(rapport,file=paste(chem_rapport,run,"_report_",format(Sys.time(), "%m-%d-%Y"),".txt",sep=""),dec=".",row.names=F,
col.names=F,sep="\t",quote=F)
