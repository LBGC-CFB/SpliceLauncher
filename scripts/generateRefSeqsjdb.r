#!/usr/bin/Rscript
options(scipen=50)

###################################################################
#Convert RefSeq BED file into sjdb file to craete genome assembly #
###################################################################
helpMessage="Usage: generateRefSeqsjdb.r\n
    [Mandatory] \n
        \t[-i|--input /path/to/inputFile]
        \t\tRefSeq BED database, downloadable at UCSC: https://genome.ucsc.edu/cgi-bin/hgTables\n
        \t[-o|--output /path/to/outputFile]
        \t\tOutput to the sjdb file\n
        \t[-h|--help]
        \t\tprint this help message and exit\n
   You could : Rscript generateRefSeqsjdb.r -i ./RefSeqAnnot.bed -o ./RefSeqAnnot.sjdb"

#get script argument
args <- commandArgs(trailingOnly = TRUE)

if (length(args)<2){message(helpMessage);stop()}

i=1
while (i <= length(args)){
    if(args[i]=="-i"|args[i]=="--input"){
        inputFile=normalizePath(path=args[i+1]);i = i+2
    }else if(args[i]=="-o"|args[i]=="--output"){
        outputFile=args[i+1];i = i+2
    }else if(args[i]=="-h"|args[i]=="--help"){
        message(helpMessage);stop()
    }else{
        message(paste("********Unknown option:",args[i],"\n"));message(helpMessage);stop()
    }
}

message(paste("Input:",inputFile))
message(paste("Output:",outputFile))

Chr=NULL
strand=NULL
StartInt = NULL
EndInt = NULL

message("import RefSeq BED file...")

dataRefSeq = read.table(inputFile,header=F, sep="\t")
dataRefSeq[,4] = as.character(dataRefSeq[,4])
a = unlist(strsplit(dataRefSeq[,4],".",fixed=T))
dataRefSeq[,4] = a[grep("_",a)]

convertBEDfile <- function(chr,posStart,sens,taille,tailleCum){

	if(sens=="+"){
		tailleCum = strsplit(as.character(tailleCum),split=",")
		tailleCum = as.numeric(unlist(tailleCum))
		taille = strsplit(as.character(taille),split=",")
		taille = as.numeric(unlist(taille))
		posAcc = posStart+tailleCum
		posDon = posStart+tailleCum+taille
		posDon = posDon[-length(posDon)]
		posAcc = posAcc[-1]

		Chr = c(Chr,rep(chr,length(posDon)))
		strand = c(strand,rep(sens,length(posDon)))
		StartInt = c(StartInt,posDon+1)
		EndInt = c(EndInt,posAcc)

	}else if(sens=="-"){
		posEnd=posStart
		tailleCum=strsplit(as.character(tailleCum),split=",")
		tailleCum=as.numeric(unlist(tailleCum))
		taille=strsplit(as.character(taille),split=",")
		taille=as.numeric(unlist(taille))
		posDon=posEnd+tailleCum
		posAcc=posEnd+tailleCum+taille
		posAcc=posAcc[-length(posAcc)]
		posDon = posDon[-1]

		Chr = c(Chr,rep(chr,length(posAcc)))
		strand = c(strand,rep(sens,length(posAcc)))
		StartInt = c(StartInt,posAcc+1)
		EndInt = c(EndInt,posDon)
	}else{
		print("error")
	}
result <<-list(Chr,StartInt,EndInt,strand)
}

message("Convert in Intron coordinates...")
tmp = mapply(convertBEDfile, chr = as.character(dataRefSeq[,1]), posStart = dataRefSeq[,2], sens = as.character(dataRefSeq[,6]), taille = dataRefSeq[,11], tailleCum = dataRefSeq[,12])

dataOutput=data.frame(chr = as.character(unlist(tmp[1,])),
					StartInt = unlist(tmp[2,]),
					EndInt = unlist(tmp[3,]),
					strand = unlist(tmp[4,])
)

message("File is saving...")
write.table(dataOutput,outputFile,quote=F, row.names=F, col.names=F,sep="\t")
message("Finish !!!")
