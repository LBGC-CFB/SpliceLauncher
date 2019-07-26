#!/usr/bin/Rscript
T1<-as.numeric(format(Sys.time(), "%s"))
options(scipen=50,stringsAsFactors=FALSE)

#Get SpliceLauncher databases

helpMessage="Usage: generateSpliceLauncherDB.r\n
    [Mandatory] \n
        -i, --input /path/to/inputFile
            RefSeq GFF annotation file, downloadable at UCSC: https://genome.ucsc.edu/cgi-bin/hgTables\n
        -o, --output /path/to/directory/
            Directory to output databases\n
    -h, --help
        print this help message and exit\n
   You could : Rscript generateRefSeqsjdb.r -i ./RefSeqAnnot.gff -o ./RefSeqAnnot.sjdb"

args <- commandArgs(trailingOnly = TRUE)

  if (length(args)<2){message(helpMessage);stop()}

i=1
while (i <= length(args)){
    if(args[i]=="-i"|args[i]=="--input"){
           inputFile=normalizePath(path=args[i+1]);i = i+2
       }else if(args[i]=="-o"|args[i]=="--output"){
           outputPath=args[i+1];i = i+2
           if (!file.exists(outputPath)){dir.create(outputPath, recursive = TRUE)}
           outputPath = normalizePath(outputPath)
       }else if(args[i]=="-h"|args[i]=="--help"){
           message(helpMessage);stop()
       }else{
           message(paste("********Unknown option:",args[i],"\n"));message(helpMessage);stop()
       }
}

############################
#Function used
############################

splitRawToTable <- function(raw, sep = "\t", columNames){
    nCol = length(columNames)
    splitRaw = unlist(strsplit(raw,sep,fixed=TRUE))
    data=as.data.frame(matrix(splitRaw, ncol = nCol,byrow = TRUE))
    colnames(data) <- columNames
    return(data)
}

extractHeader <- function(info,sep=";",sep2="="){
    tmp = unlist(strsplit(info,sep,fixed=TRUE))
    tmp = unlist(strsplit(tmp,sep2,fixed=TRUE))
    header = tmp[seq(from=1,to=length(tmp),by=2)]
    return(unique(header))
}

extractInfo <- function(info,sep="=",header){
    tmp = unlist(strsplit(info,sep,fixed=TRUE))
    idH = tmp[seq(from=1,to=length(tmp),by=2)]
    splitInfo = tmp[seq(from=2,to=length(tmp),by=2)]
    raw=rep(".",length(header))
    names(raw)=header
    raw[idH]=splitInfo
    return(raw)
}

convertNCtoChr <- function (NCname){
    tmp = unlist(strsplit(NCname,".",fixed=TRUE))
    tmp = tmp[seq(from=1,to=length(tmp),by=2)]
    Chrname = paste('chr',as.numeric(substr(tmp,8,9)),sep="")
    Chrname = sub(pattern = "chr23", replacement="chrX", Chrname,fixed = TRUE)
    Chrname = sub(pattern = "chr24", replacement="chrY", Chrname,fixed = TRUE)
    return(Chrname)
}

getBEDfile <- function(exonCoord){
    tmp = data.frame(V1 = exonCoord$chr,
        V2 = as.numeric(exonCoord$start)-1,
        V3 = exonCoord$end,
        V4 = paste(
                substr(exonCoord$transcript_id,rep(1,length(exonCoord$transcript_id)),as.numeric(regexpr(".",exonCoord$transcript_id,fixed =TRUE))-1),
                "exon",0,0,exonCoord$chr,exonCoord$start,exonCoord$strand,sep="_"),
        V5 = rep(0,nrow(exonCoord)),
        V6 = exonCoord$strand)
    tmp[,4] = sub(pattern = "_-", replacement="_r", tmp[,4],fixed = TRUE)
    tmp[,4] = sub(pattern = "_+", replacement="_f", tmp[,4],fixed = TRUE)
    return(tmp)
}

getSJDBfile <- function(exonCoord){
    tmp = data.frame(V1 = exonCoord$chr[-nrow(exonCoord)],
        V2 = rep(NA,nrow(exonCoord)-1),
        V3 = rep(NA,nrow(exonCoord)-1),
        V4 = exonCoord$strand[-nrow(exonCoord)])
    # print(summary(tmp))
    tmp$V2[tmp$V4=="+"] = as.numeric(exonCoord$end[which(tmp$V4=="+")])+1
    tmp$V3[tmp$V4=="+"] = as.numeric(exonCoord$start[which(tmp$V4=="+")+1])-1
    tmp$V2[tmp$V4=="-"] = as.numeric(exonCoord$end[which(tmp$V4=="-")+1])+1
    tmp$V3[tmp$V4=="-"] = as.numeric(exonCoord$start[which(tmp$V4=="-")])-1
    t1 = exonCoord$transcript_id[-nrow(exonCoord)]
    t2 = exonCoord$transcript_id[-1]
    tmp = tmp[-which(t1!=t2),]
    # print(summary(tmp))
    return(tmp)
}

extractCDS<-function(proteinInfo){
    sp1 = proteinInfo$Parent
    sp2 = c("toto",proteinInfo$Parent[-nrow(proteinInfo)])
    ep1 = proteinInfo$Parent
    ep2 = c(proteinInfo$Parent[-1],"toto")
    idstart <- sp1!=sp2
    idend <- ep1!=ep2
    tmp = proteinInfo[idstart,]
    tmp$start = as.numeric(tmp$start)-1
    tmp$end = proteinInfo$end[idend]
    tmp$start[tmp$strand=="-"] = as.numeric(proteinInfo$start[proteinInfo$strand=="-" & idend])-1
    tmp$end[tmp$strand=="-"] = proteinInfo$end[proteinInfo$strand=="-" & idstart]
    return(tmp)
}

getSpliceLauncherDB <- function(transcrit){
    tmp = tmpData[tmpData$transcript_id==transcrit,]
    tmp = tmp[order(tmp$start),]
    printInfo=TRUE
    if(nrow(tmp)<2){
        printInfo=FALSE
    }

    if(printInfo){
        Strand = tmp$strand[1]
        tailleExon = abs(as.numeric(tmp$start)- as.numeric(tmp$end))

        if(Strand =="+"){
            gCDSstart = tmp$start_CDS[1]
            gCDSend = tmp$end_CDS[1]
            idEx = c(1:length(tailleExon ))
            lenEx = tailleExon
            gStart = tmp$start
            gEnd = tmp$end
            cStart = rep(0,length(tailleExon ))
            cEnd = rep(0,length(tailleExon ))

            ExCDSstart=idEx[gStart<=gCDSstart & gEnd>=gCDSstart]
            ExCDSend=idEx[gStart<=gCDSend & gEnd>=gCDSend ]
            if(length(ExCDSstart)>1 | length(ExCDSend)>1){
                for(i in 2:length(gStart)){
                    if(gStart[i]==gEnd[i-1]){
                        gStart[i] = gStart[i]+1
                        ExCDSstart=idEx[gStart<=gCDSstart & gEnd>=gCDSstart]
                        ExCDSend=idEx[gStart<=gCDSend & gEnd>=gCDSend ]
                    }
                }
                if(length(ExCDSstart)>1 | length(ExCDSend)>1){
                    printInfo=FALSE
                }
            }
            if(printInfo){
                cStart[idEx==ExCDSstart]=gStart[idEx==ExCDSstart]-gCDSstart
                cEnd[idEx==ExCDSstart]=gEnd[idEx==ExCDSstart]-gCDSstart

                if(ExCDSstart<ExCDSend){
                    if(ExCDSstart < max(idEx)){
                        for (i in seq(from=ExCDSstart+1,to=ExCDSend,by=1)){
                            cStart[idEx==i]=cEnd[idEx==(i-1)]+1
                            cEnd[idEx==i]=cStart[idEx==i]+(lenEx[idEx==i]-1)
                        }
                    }
                }
                cEnd[idEx==ExCDSend]=gEnd[idEx==ExCDSend]-gCDSend +1

                if(ExCDSstart > 1){
                    for (i in seq(from=ExCDSstart-1,to=1,by=-1)){
                        cEnd[idEx==i]=cStart[idEx==(i+1)]-1
                        cStart[idEx==i]=cEnd[idEx==i]-(lenEx[idEx==i]-1)
                    }
                }
                if(ExCDSend < max(idEx)){
                    for (i in seq(from=ExCDSend+1,to=max(idEx),by=1)){
                        cStart[idEx==i]=cEnd[idEx==i-1]+1
                        cEnd[idEx==i]=cStart[idEx==i]+(lenEx[idEx==i]-1)
                    }
                }
            }
        }else if(Strand=="-"){
            gCDSstart2 = tmp$end_CDS[1]
            gCDSend2 = tmp$start_CDS[1]
            idEx = c(length(tailleExon):1)
            lenEx = tailleExon
            gStart = tmp$end
            gEnd = tmp$start
            cStart = rep(0,length(tailleExon))
            cEnd = rep(0,length(tailleExon))

            ExCDSstart=idEx[gStart>=gCDSstart2 & gEnd<=gCDSstart2]
            ExCDSend=idEx[gStart>=gCDSend2 & gEnd<=gCDSend2 ]
            if(length(ExCDSstart)>1 | length(ExCDSend)>1){
                for(i in 2:length(gStart)){
                    if(gEnd[i]==gStart[i-1]){
                        gEnd[i] = gEnd[i]-1
                        ExCDSstart=idEx[gStart>=gCDSstart2 & gEnd<=gCDSstart2]
                        ExCDSend=idEx[gStart>=gCDSend2 & gEnd<=gCDSend2 ]
                    }
                }
                if(length(ExCDSstart)>1 | length(ExCDSend)>1){
                    printInfo=FALSE
                }
            }
            if(printInfo){
                cStart[idEx==ExCDSstart]=gCDSstart2-gStart[idEx==ExCDSstart]
                cEnd[idEx==ExCDSstart]=gCDSstart2-gEnd[idEx==ExCDSstart]

                if(ExCDSstart<ExCDSend){
                    if(ExCDSstart < max(idEx)){
                        for (i in seq(from=ExCDSstart+1,to=ExCDSend,by=1)){
                            cStart[idEx==i]=cEnd[idEx==(i-1)]+1
                            cEnd[idEx==i]=cStart[idEx==i]+(lenEx[idEx==i]-1)
                        }
                    }
                }
                cEnd[idEx==ExCDSend]=gCDSend2-gEnd[idEx==ExCDSend]+1

                if(ExCDSstart > 1){
                    for (i in seq(from=ExCDSstart-1,to=1,by=-1)){
                        cEnd[idEx==i]=cStart[idEx==(i+1)]-1
                        cStart[idEx==i]=cEnd[idEx==i]-(lenEx[idEx==i]-1)
                    }
                }
                if(ExCDSend < max(idEx)){
                    for (i in seq(from=ExCDSend+1,to=max(idEx),by=1)){
                        cStart[idEx==i]=cEnd[idEx==i-1]+1
                        cEnd[idEx==i]=cStart[idEx==i]+(lenEx[idEx==i]-1)
                    }
                }
                gCDSstart = gCDSstart2
                gCDSend = gCDSend2
            }
        }
        if(printInfo){
            result <<-data.frame(Gene=tmp$gene,Strand=tmp$strand,gCDSstart=rep(gCDSstart,length(tailleExon)),
                gCDSend=rep(gCDSend,length(tailleExon)),transcrit=tmp$transcript_id,Chr = tmp$chr,
                idEx,lenEx,gStart,gEnd,cStart,cEnd)
        }
    }
}

###################
#RefSeq
###################
message("Read GFF file ...")
gffData =readLines(inputFile)

#Remove bad chr
message("Remove meta-header ...")
fileType = gffData[1]
if(substr(fileType,1,5)!="##gff"){
    stop(paste("Please select a GFF file because your file is:",fileType))
}
if(as.numeric(substr(fileType,nchar(fileType)-1,nchar(fileType)))<3){
    stop(paste("Import GFF file in version 3 or later, your version is:",as.numeric(substr(fileType,nchar(fileType)-1,nchar(fileType)))))
}
mHead = gffData[grep("#",substr(gffData,1,2),fixed=TRUE)]

# message(paste("Your database:",paste(mHead,collapse="\n"),sep="\n"))


gffData = gffData[-grep("#",substr(gffData,1,2),fixed=TRUE)]

#Change in table
message("Formate GFF file ...")
columNames = c("chr","status","seqType","start","end","X1","strand","X2","info")
gffData = splitRawToTable(gffData, sep = "\t", columNames)

#Convert chromosome annotation
message("Change chr annotation ...")
if(substr(gffData$chr[1],1,3)=="chr"){
    gffData = gffData[grep("chr",gffData$chr,fixed=TRUE),]
}else if(substr(gffData$chr[1],1,3)=="NC_"){
    gffData = gffData[grep("NC_000",gffData$chr,fixed=TRUE),]
}
if(length(grep('gene_name',gffData$info[gffData$seqType=='exon'][1]))>0){
    message("Change gene tag ...")
    gffData$info = sub('gene_name','gene',gffData$info)
}

# Extract the exon and CDS DB
message("Extract exon coordinates ...")
exonCoord = gffData[gffData$seqType=='exon',]
if(nrow(exonCoord)==0){
    stop("Your file has not exon information")
}
message("Extract exon annotation ...")
colNames = extractHeader(exonCoord$info)
print(colNames)
exonCoordtmp = strsplit(exonCoord$info, ";", fixed=TRUE)
tmp = unlist(lapply(exonCoordtmp,extractInfo,header = colNames))
exonCoordtmp = as.data.frame(matrix(tmp, ncol = length(colNames),byrow = TRUE))
names(exonCoordtmp) = colNames
message("Merge exon files ...")
exonCoord = cbind(exonCoord,exonCoordtmp)
exonCoord = exonCoord[exonCoord$transcript_id!=".",]

message("Extract CDS coordinates ...")
proteinInfo = gffData[gffData$seqType=='CDS',]
if(nrow(proteinInfo)==0){
    stop("Your file has not CDS information")
}
message("Extract CDS annotation ...")
colNames = extractHeader(proteinInfo$info)
print(colNames)
proteinInfotmp = strsplit(proteinInfo$info, ";", fixed=TRUE)
tmp = unlist(lapply(proteinInfotmp,extractInfo,header = colNames))
proteinInfotmp = as.data.frame(matrix(tmp, ncol = length(colNames),byrow = TRUE))
names(proteinInfotmp) = colNames
message("Merge CDS files ...")
proteinInfo = cbind(proteinInfo,proteinInfotmp)
proteinInfo = extractCDS(proteinInfo)

# get BED file
message("Generate BED annotation ...")
bedFile = getBEDfile(exonCoord)
write.table(bedFile,paste(outputPath,"BEDannotation.bed",sep="/"),sep="\t",quote=F,row.names=F,col.names=F)

# Get sjdb file
message("Generate SJDB file ...")
sjdbFile = getSJDBfile(exonCoord)
write.table(sjdbFile,paste(outputPath,"SJDBannotation.sjdb",sep="/"),sep="\t",quote=F,row.names=F,col.names=F)

# Get SpliceLauncher DB
message("Generate SpliceLauncher Reference file ...")
message("   Merge exon and CDS annotation ...")
mergeData = merge(exonCoord,proteinInfo[,c("start","end","Parent")],by="Parent",all.x=TRUE, suffixes = c("","_CDS"),sort=FALSE)
mergeData = mergeData[,c("Parent","chr","start","end","strand","ID","gene","transcript_id","exon_number","start_CDS","end_CDS")]
tmpTranscript = unlist(strsplit(mergeData$transcript_id,".",fixed=TRUE))
mergeData$transcript_id = tmpTranscript[nchar(tmpTranscript)>5]
mergeData$start = as.numeric(mergeData$start)-1
mergeData$end = as.numeric(mergeData$end)
mergeData$start_CDS = as.numeric(mergeData$start_CDS)
mergeData$end_CDS = as.numeric(mergeData$end_CDS)
mergeData$start_CDS[is.na(mergeData$start_CDS)] = mergeData$start[is.na(mergeData$start_CDS)]
mergeData$end_CDS[is.na(mergeData$end_CDS)] = mergeData$end[is.na(mergeData$end_CDS)]

message("   Get SpliceLauncher Reference file ...")
trans = unique(mergeData$transcript_id)
nIter = length(trans)%/%5000
print(length(trans))

j=1
for(i in 1:nIter){
    indStart = j
    indEnd = j + 4999
    message(paste("   Get transcript info for:",j,"to",j+4999))
    j = j+5000
    tmpTrans = trans[indStart:indEnd]
    tmpData = mergeData[which(mergeData$transcript_id%in%tmpTrans),]
    tmp <- unlist(mapply(getSpliceLauncherDB, transcrit = tmpTrans))
    dataConvertPool = data.frame(Gene = tmp[grep("Gene",names(tmp))],
                                Strand = tmp[grep("Strand",names(tmp))],
                                gCDSstart = tmp[grep("gCDSstart",names(tmp))],
                                gCDSend = tmp[grep("gCDSend",names(tmp))],
                                transcrit = tmp[grep("transcrit",names(tmp))],
                                Chr = tmp[grep("Chr",names(tmp))],
                                idEx = tmp[grep("idEx",names(tmp))],
                                lenEx = tmp[grep("lenEx",names(tmp))],
                                gStart = tmp[grep("gStart",names(tmp))],
                                gEnd = tmp[grep("gEnd",names(tmp))],
                                cStart = tmp[grep("cStart",names(tmp))],
                                cEnd = tmp[grep("cEnd",names(tmp))])
    if(i==1){
        write.table(dataConvertPool,paste(outputPath,"SpliceLauncherAnnot.txt",sep="/"),sep="\t",quote=F,row.names=F,col.names=TRUE,append=FALSE)
    }else{
        write.table(dataConvertPool,paste(outputPath,"SpliceLauncherAnnot.txt",sep="/"),sep="\t",quote=F,row.names=F,col.names=F,append=TRUE)
    }
}

rest = length(trans)%%5000

tmpTrans = trans[(nIter*5000+1):(nIter*5000+rest)]
tmpData = mergeData[which(mergeData$transcript_id%in%tmpTrans),]
tmp <- unlist(mapply(getSpliceLauncherDB, transcrit = unique(tmpData$transcript_id)))
message("   Formate SpliceLauncher file...")
dataConvertPool = data.frame(Gene = tmp[grep("Gene",names(tmp))],
                            Strand = tmp[grep("Strand",names(tmp))],
                            gCDSstart = tmp[grep("gCDSstart",names(tmp))],
                            gCDSend = tmp[grep("gCDSend",names(tmp))],
                            transcrit = tmp[grep("transcrit",names(tmp))],
                            Chr = tmp[grep("Chr",names(tmp))],
                            idEx = tmp[grep("idEx",names(tmp))],
                            lenEx = tmp[grep("lenEx",names(tmp))],
                            gStart = tmp[grep("gStart",names(tmp))],
                            gEnd = tmp[grep("gEnd",names(tmp))],
                            cStart = tmp[grep("cStart",names(tmp))],
                            cEnd = tmp[grep("cEnd",names(tmp))])

write.table(dataConvertPool,paste(outputPath,"SpliceLauncherAnnot.txt",sep="/"),sep="\t",quote=F,row.names=F,col.names=F,append=TRUE)

T2<-as.numeric(format(Sys.time(), "%s"))
print(T2-T1)
