#!/usr/bin/Rscript
options(scipen=50)

tryCatch({
library("optparse")
},
	error=function(cond) {
		message("Here's the original error message:")
		message(cond)
		message("*****You need to install \'optparse\' library")
})

#############################################################
#Convert RefSeq BED file into BED file for BEDTools analysis#
#############################################################

option_list = list(make_option(opt_str = c("-i","--input"), action="store", type="character", default=NULL,
						help="RefSeq BED database", metavar="character"),
					make_option(c("-o",'--output'), type="character", default=NULL,
						help="output to exon BED file", metavar="character"),
					 make_option(c("-t",'--transcript'), type="character", default="no file",
						help="target transcripts in txt file  [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if(is.null(opt$input)|is.null(opt$output)){
	print_help(opt_parser)
	stop()
}

inputFile <- opt$input
outputFile <- opt$output

message("import RefSeq BED file...")
dataRefSeq = read.table(inputFile,header=F, sep="\t")
dataRefSeq[,4] = as.character(dataRefSeq[,4])
a = unlist(strsplit(dataRefSeq[,4],".",fixed=T))
dataRefSeq[,4] = a[grep("_",a)]

if(opt$transcript!="no file"){
	message("select the targeted transcripts...")
	SelectTranscrit = readLines(opt$transcript)
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
