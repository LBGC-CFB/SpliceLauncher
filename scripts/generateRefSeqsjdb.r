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

###################################################################
#Convert RefSeq BED file into sjdb file to craete genome assembly #
###################################################################

option_list = list(make_option(opt_str = c("-i","--input"), action="store", type="character", default=NULL,
						help="RefSeq BED database", metavar="character"),
					make_option(c("-o",'--output'), type="character", default=NULL,
						help="output to sjdb file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if(is.null(opt$input)|is.null(opt$output)){
	print_help(opt_parser)
	stop()
}

inputFile <- opt$input
outputFile <- opt$output

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
