#!/usr/bin/Rscript
options(scipen=50)

#############################################################
#Convert RefSeq BED file into BED file for BEDTools analysis#
#############################################################
helpMessage="Usage: generateExonBEDRef.r\n
    [Mandatory] \n
        \t[-i|--input /path/to/inputFile]
        \t\tRefSeq BED database, downloadable at UCSC: https://genome.ucsc.edu/cgi-bin/hgTables\n
        \t[-o|--output /path/to/outputFile]
        \t\tOutput to exon BED file\n
    [Options] \n
        \t[-t|--transcript]
        \t\tList of target transcripts in txt file\n
        \t[-h|--help]
        \t\tprint this help message and exit\n
   You could : Rscript generateExonBEDRef.r -i ./RefSeqAnnot.bed -o ./RefExons.bed"

#get script argument
transcript=NULL
args <- commandArgs(trailingOnly = TRUE)

if (length(args)<2){message(helpMessage);stop()}

i=1
while (i <= length(args)){
    if(args[i]=="-i"|args[i]=="--input"){
        inputFile=normalizePath(path=args[i+1]);i = i+2
    }else if(args[i]=="-o"|args[i]=="--output"){
        outputFile=args[i+1];i = i+2
    }else if(args[i]=="-t"|args[i]=="--MergeTranscrit"){
        transcript=normalizePath(args[i+1]);i = i+2
    }else if(args[i]=="-h"|args[i]=="--help"){
        message(helpMessage);stop()
    }else{
        message(paste("********Unknown option:",args[i],"\n"));message(helpMessage);stop()
    }
}

message(paste("Input:",inputFile))
message(paste("Output:",outputFile))

message("import RefSeq BED file...")
dataRefSeq = read.table(inputFile,header=F, sep="\t")
dataRefSeq[,4] = as.character(dataRefSeq[,4])
a = unlist(strsplit(dataRefSeq[,4],".",fixed=T))
dataRefSeq[,4] = a[grep("_",a)]

if(!is.null(transcript)){
	message("select the targeted transcripts...")
	SelectTranscrit = readLines(transcript)
	dataRefSeq = dataRefSeq[which(dataRefSeq$V4%in%SelectTranscrit),]
}

Chr=NULL
strand=NULL
nbEx=NULL
Transcrit=NULL
posEX=NULL
StartEx = NULL
EndEx = NULL

convertBEDfile <- function(chr,posStart,transcrit,sens,taille,tailleCum){

	if(sens=="+"){
		tailleCum = strsplit(as.character(tailleCum),split=",")
		tailleCum = as.numeric(unlist(tailleCum))
		taille = strsplit(as.character(taille),split=",")
		taille = as.numeric(unlist(taille))
		posAcc = posStart+tailleCum
		posDon = posStart+tailleCum+taille

		Chr = c(Chr,rep(chr,length(posDon)))
		strand = c(strand,rep(sens,length(posDon)))
		Transcrit = c(Transcrit,rep(transcrit,length(posDon)))
		posEX = c(posEX,posAcc+1)
		StartEx = c(StartEx,posAcc)
		EndEx = c(EndEx,posDon)

	}else if(sens=="-"){
		posEnd=posStart
		tailleCum=strsplit(as.character(tailleCum),split=",")
		tailleCum=as.numeric(unlist(tailleCum))
		taille=strsplit(as.character(taille),split=",")
		taille=as.numeric(unlist(taille))
		posDon=posEnd+tailleCum
		posAcc=posEnd+tailleCum+taille

		Chr = c(Chr,rep(chr,length(posAcc)))
		strand = c(strand,rep(sens,length(posAcc)))
		Transcrit = c(Transcrit,rep(transcrit,length(posAcc)))
		posEX = c(posEX,posDon+1)
		StartEx = c(StartEx,posDon)
		EndEx = c(EndEx,posAcc)
	}else{
		print("error")
	}
result <<-list(Chr,StartEx,EndEx,Transcrit,posEX,strand)
}

message("Convert in Exon...")
tmp = mapply(convertBEDfile, chr = as.character(dataRefSeq[,1]), posStart = dataRefSeq[,2], transcrit = as.character(dataRefSeq[,4]), sens = as.character(dataRefSeq[,6]), taille = dataRefSeq[,11], tailleCum = dataRefSeq[,12])

dataOutput=data.frame(chr = as.character(unlist(tmp[1,])),
					StartEx = unlist(tmp[2,]),
					EndEx = unlist(tmp[3,]),
					Annot = paste(unlist(tmp[4,]),"exon",0,0,unlist(tmp[1,]),unlist(tmp[5,]), unlist(tmp[6,]),sep="_"),
					Score = rep(0,length(unlist(tmp[1,]))),
					strand = unlist(tmp[6,])
)

dataOutput[,4] = sub(pattern = "_-", replacement="_r", dataOutput[,4],fixed = TRUE)
dataOutput[,4] = sub(pattern = "_+", replacement="_f", dataOutput[,4],fixed = TRUE)

message("File is saving...")
write.table(dataOutput,outputFile,quote=F, row.names=F, col.names=F,sep="\t")

message("Finish !!!")
