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

##########################################################
#Convert RefSeq database into Ref file for SpliceLauncher#
##########################################################

option_list = list(make_option(opt_str = c("-i","--input"), action="store", type="character", default=NULL,
						help="RefSeq txt database", metavar="character"),
					make_option(c("-o",'--output'), type="character", default=NULL,
						help="output to SpliceLauncher file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if(is.null(opt$input)|is.null(opt$output)){
	print_help(opt_parser)
	stop()
}
inputFile <- opt$input
outputFile <- opt$output

message("Import RefSeq data...")
dataRefSeq = read.table(inputFile,header=F, sep="\t")
names(dataRefSeq) = c("bin","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds","score","name2","cdsStartStat","cdsEndStat","exonFrames")

dataRefSeq = na.omit(dataRefSeq)
dataConvertPool=NULL
a = unlist(strsplit(as.character(dataRefSeq$name),".",fixed=T))
dataRefSeq$name = a[grep("_",a)]

Gene_pool = NULL
Strand_pool = NULL
gCDSstart_pool = NULL
gCDSend_pool = NULL
transcrit_pool = NULL
Chr_pool = NULL
idEx_pool = NULL
lenEx_pool = NULL
gStart_pool = NULL
gEnd_pool = NULL
cStart_pool = NULL
cEnd_pool = NULL

convertcNomenIngNomen <- function(transcrit = dataRefSeq[,"name"], sens = dataRefSeq[,"strand"],
									gene = dataRefSeq[,"name2"], chr = dataRefSeq[,"chrom"],
									exonStarts = dataRefSeq[,"exonStarts"], exonEnds = dataRefSeq[,"exonEnds"],
									gCDSstart = dataRefSeq[,"cdsStart"], gCDSend = dataRefSeq[,"cdsEnd"]){

	tailleExon = abs(as.numeric(unlist(strsplit(as.character(exonStarts),",")))-
						as.numeric(unlist(strsplit(as.character(exonEnds),","))))

	if(sens=="+"){

		posAcc = as.numeric(unlist(strsplit(as.character(exonStarts),",")))
		posDon = as.numeric(unlist(strsplit(as.character(exonEnds),",")))

		dataConvert=data.frame(Gene=rep(gene,length(tailleExon)),Strand=rep(sens,length(tailleExon)),gCDSstart=rep(gCDSstart,length(tailleExon)),
							gCDSend=rep(gCDSend,length(tailleExon)),transcrit=rep(transcrit,length(tailleExon)),Chr = rep(chr,length(tailleExon)),
							idEx=c(1:length(tailleExon )),lenEx=tailleExon,gStart=posAcc,gEnd=posDon,cStart=0,cEnd=0 )

		ExCDSstart=dataConvert$idEx[dataConvert$gStart<=gCDSstart & dataConvert$gEnd>=gCDSstart]
		ExCDSend=dataConvert$idEx[dataConvert$gStart<=gCDSend & dataConvert$gEnd>=gCDSend ]

		dataConvert$cStart[dataConvert$idEx==ExCDSstart]=dataConvert$gStart[dataConvert$idEx==ExCDSstart]-gCDSstart
		dataConvert$cEnd[dataConvert$idEx==ExCDSstart]=dataConvert$gEnd[dataConvert$idEx==ExCDSstart]-gCDSstart
		if(ExCDSstart<ExCDSend){
			if(ExCDSstart < max(dataConvert$idEx)){
				for (i in seq(from=ExCDSstart+1,to=ExCDSend,by=1)){
					dataConvert$cStart[dataConvert$idEx==i]=dataConvert$cEnd[dataConvert$idEx==(i-1)]+1
					dataConvert$cEnd[dataConvert$idEx==i]=dataConvert$cStart[dataConvert$idEx==i]+(dataConvert$lenEx[dataConvert$idEx==i]-1)
				}
			}
		}
		dataConvert$cEnd[dataConvert$idEx==ExCDSend]=dataConvert$gEnd[dataConvert$idEx==ExCDSend]-gCDSend +1

		if(ExCDSstart > 1){
			for (i in seq(from=ExCDSstart-1,to=1,by=-1)){
				dataConvert$cEnd[dataConvert$idEx==i]=dataConvert$cStart[dataConvert$idEx==(i+1)]-1
				dataConvert$cStart[dataConvert$idEx==i]=dataConvert$cEnd[dataConvert$idEx==i]-(dataConvert$lenEx[dataConvert$idEx==i]-1)
			}
		}
		if(ExCDSend < max(dataConvert$idEx)){
			for (i in seq(from=ExCDSend+1,to=max(dataConvert$idEx),by=1)){
				dataConvert$cStart[dataConvert$idEx==i]=dataConvert$cEnd[dataConvert$idEx==i-1]+1
				dataConvert$cEnd[dataConvert$idEx==i]=dataConvert$cStart[dataConvert$idEx==i]+(dataConvert$lenEx[dataConvert$idEx==i]-1)
			}
		}
	}else if(sens=="-"){

		gCDSstart2 = gCDSend
		gCDSend2 = gCDSstart
		posDon = as.numeric(unlist(strsplit(as.character(exonStarts),",")))
		posAcc = as.numeric(unlist(strsplit(as.character(exonEnds),",")))

		dataConvert=data.frame(Gene=rep(gene,length(tailleExon)),Strand=rep(sens,length(tailleExon)),gCDSstart2=rep(gCDSstart2,length(tailleExon)),
							gCDSend2=rep(gCDSend2,length(tailleExon)),transcrit=rep(transcrit,length(tailleExon)),Chr = rep(chr,length(tailleExon)),
							idEx=c(length(tailleExon):1),lenEx=tailleExon,gStart=posAcc,gEnd=posDon,cStart=0,cEnd=0)

		ExCDSstart=dataConvert$idEx[dataConvert$gStart>=gCDSstart2 & dataConvert$gEnd<=gCDSstart2]
		ExCDSend=dataConvert$idEx[dataConvert$gStart>=gCDSend2 & dataConvert$gEnd<=gCDSend2 ]

		dataConvert$cStart[dataConvert$idEx==ExCDSstart]=gCDSstart2-dataConvert$gStart[dataConvert$idEx==ExCDSstart]
		dataConvert$cEnd[dataConvert$idEx==ExCDSstart]=gCDSstart2-dataConvert$gEnd[dataConvert$idEx==ExCDSstart]
		if(ExCDSstart<ExCDSend){
			if(ExCDSstart < max(dataConvert$idEx)){
				for (i in seq(from=ExCDSstart+1,to=ExCDSend,by=1)){
					dataConvert$cStart[dataConvert$idEx==i]=dataConvert$cEnd[dataConvert$idEx==(i-1)]+1
					dataConvert$cEnd[dataConvert$idEx==i]=dataConvert$cStart[dataConvert$idEx==i]+(dataConvert$lenEx[dataConvert$idEx==i]-1)
				}
			}
		}
		dataConvert$cEnd[dataConvert$idEx==ExCDSend]=gCDSend2-dataConvert$gEnd[dataConvert$idEx==ExCDSend]+1

		if(ExCDSstart > 1){
			for (i in seq(from=ExCDSstart-1,to=1,by=-1)){
				dataConvert$cEnd[dataConvert$idEx==i]=dataConvert$cStart[dataConvert$idEx==(i+1)]-1
				dataConvert$cStart[dataConvert$idEx==i]=dataConvert$cEnd[dataConvert$idEx==i]-(dataConvert$lenEx[dataConvert$idEx==i]-1)
			}
		}
		if(ExCDSend < max(dataConvert$idEx)){
			for (i in seq(from=ExCDSend+1,to=max(dataConvert$idEx),by=1)){
				dataConvert$cStart[dataConvert$idEx==i]=dataConvert$cEnd[dataConvert$idEx==i-1]+1
				dataConvert$cEnd[dataConvert$idEx==i]=dataConvert$cStart[dataConvert$idEx==i]+(dataConvert$lenEx[dataConvert$idEx==i]-1)
			}
		}
		dataConvert$gCDSstart = dataConvert$gCDSstart2
		dataConvert$gCDSend = dataConvert$gCDSend2
	}

	Gene_pool = c(Gene_pool,as.character(dataConvert$Gene))
	Strand_pool = c(Strand_pool,as.character(dataConvert$Strand))
	gCDSstart_pool = c(gCDSstart_pool,dataConvert$gCDSstart)
	gCDSend_pool = c(gCDSend_pool,dataConvert$gCDSend)
	transcrit_pool = c(transcrit_pool,as.character(dataConvert$transcrit))
	Chr_pool = c(Chr_pool,as.character(dataConvert$Chr))
	idEx_pool = c(idEx_pool,dataConvert$idEx)
	lenEx_pool = c(lenEx_pool,dataConvert$lenEx)
	gStart_pool = c(gStart_pool,dataConvert$gStart)
	gEnd_pool = c(gEnd_pool,dataConvert$gEnd)
	cStart_pool = c(cStart_pool,dataConvert$cStart)
	cEnd_pool = c(cEnd_pool,dataConvert$cEnd)

	result <<-list(Gene_pool,Strand_pool,gCDSstart_pool,gCDSend_pool,transcrit_pool,Chr_pool,idEx_pool,lenEx_pool,gStart_pool,gEnd_pool,cStart_pool,cEnd_pool)
}

message("Convert to SpliceLauncher file...")
tmp <- mapply(convertcNomenIngNomen, transcrit = as.character(dataRefSeq[,"name"]), sens = as.character(dataRefSeq[,"strand"]), gene = as.character(dataRefSeq[,"name2"]), chr = as.character(dataRefSeq[,"chrom"]),
								exonStarts = dataRefSeq[,"exonStarts"], exonEnds =dataRefSeq[,"exonEnds"], gCDSstart = dataRefSeq[,"cdsStart"], gCDSend = dataRefSeq[,"cdsEnd"])

dataConvertPool = data.frame(Gene = unlist(tmp[1,]),
							Strand = unlist(tmp[2,]),
							gCDSstart = unlist(tmp[3,]),
							gCDSend = unlist(tmp[4,]),
							transcrit = unlist(tmp[5,]),
							Chr = unlist(tmp[6,]),
							idEx = unlist(tmp[7,]),
							lenEx = unlist(tmp[8,]),
							gStart = unlist(tmp[9,]),
							gEnd = unlist(tmp[10,]),
							cStart = unlist(tmp[11,]),
							cEnd = unlist(tmp[12,]))

message("File is saving...")
write.table(dataConvertPool,outputFile,sep="\t",quote=F,row.names=F)
message("Finish !!!")
