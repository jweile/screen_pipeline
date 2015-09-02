library("hash")

read.stream <- function(result.dir) {
	data <- list()
	con <- pipe(paste("zcat ",result.dir,"all_counts.txt.gz",sep=""),open="rb")
	while(length(line <- readLines(con,1)) > 0) {
		if (substr(line,1,1)=="#") {
			sample.id <- substr(line,2,nchar(line))
			if (!(sample.id %in% names(data))) {
				data[[sample.id]] <- hash()
			}
		} else if (nchar(line) > 0) {
			key.val <- strsplit(line,":")[[1]]
			key <- key.val[[1]]
			val <- key.val[[2]]
			if (!has.key(key,data[[sample.id]])) {
				data[[sample.id]][[key]] <- as.integer(val)
			} else {
				data[[sample.id]][[key]] <- data[[sample.id]][[key]] + as.integer(val)
			}
		}
	}
	close(con)
	data
}

# result.dir <- "ube2i_compl_2015-07-21_09-44-33/UBE2I-Complementation_R1/"
# result.dir <- "ube2i_compl_sec_2015-07-30_08-41-21/UBE2I-Complementation_R1/"
# result.dir <- "ube2i_compl_sec_2015-08-19_10-46-33/UBE2I-Complementation_R1/"
result.dir <- "ube2i_y2h_2015-08-25_08-40-30/edgotyping_R1/"

# muxtag.db <- "res/muxtags"
muxtag.db <- "res/muxtags_y2h"

sample.table <- read.delim(paste(muxtag.db,"_samples.tsv",sep=""),stringsAsFactors=FALSE)

data <- read.stream(result.dir)
# data <- data[c(as.character(1:18),"invalid","undetermined")]
data <- data[c(as.character(1:72),"invalid","undetermined")]

demux.stats <- do.call(rbind,lapply(data,function(counts) {
	if (is.null(counts) || is.na(counts)) {
		return(c(unknown.BC.counts=0,known.BC.counts=0))
	}
	na <- counts[["NA"]]
	c(unknown.BC.counts=na,known.BC.counts=sum(values(counts))-na)
}))

clones <- Reduce(union,lapply(data,keys))
data.mat <- do.call(cbind,lapply(data,function(counts) {
	sapply(clones,function(clone) {
		if (has.key(clone,counts)) counts[[clone]] else 0
	})
}))
colnames(data.mat) <- names(data)
rownames(data.mat) <- clones

write.table(data.mat,"raw_counts.csv",sep=",")

data.mat2 <- data.mat[-1,1:(ncol(data.mat)-2)]

data.rel <- apply(data.mat2,2,function(x)x/sum(x))


quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
   
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
   
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}
#apply quantile normalization to replicate groups
invisible(lapply(with(sample.table,tapply(sample,paste(interactor,condition,sep="_"),c)), function(is) {
	data.rel[,is] <<- quantile_normalisation(data.rel[,is])
}))


#BOXPLOT
ord <- with(sample.table,order(interactor,condition,replicate))
plotcols <- do.call(c,(lapply(do.call(c,lapply(c("firebrick","darkolivegreen","steelblue","peachpuff","mediumpurple","aquamarine","darkgoldenrod","lightsteelblue"),function(x) paste(x,c(1,3,4),sep=""))),rep,3)))
boxplot(data.rel[,ord],ylim=c(0,.0003),col=plotcols,ylab="rel.freq.")
mtext(unique(sample.table$interactor[ord]),side=1,line=3,at=0:7*9+4.5)
mtext(c("-HIS","+HIS","+3AT"),side=1,line=2,at=seq(2,71,3))


#Load clones
clone.table <- read.csv("res/clones_y2h.csv",stringsAsFactors=FALSE)
clone.idx <- hash(clone.table$id,1:nrow(clone.table))
muts <- strsplit(clone.table$aa.calls,",")


delpos <- do.call(rbind,lapply(strsplit(clone.table$deletions,"-"),function(x) if (length(x)==1&&is.na(x)) c(NA,NA) else as.numeric(x)))
foo <- delpos[,1] < 82 & delpos[,2] > 545
foo[is.na(foo)] <- FALSE
null.clones <- clone.table$id[foo]
foo <- delpos[,1] < delpos[,2]
foo[is.na(foo)] <- FALSE
longdel.clones <- clone.table$id[foo]

permissive.cols <- with(sample.table,sample[interactor=="ZBED1" & condition=="+HIS"])
selective.cols <- with(sample.table,sample[interactor=="ZBED1" & condition=="-HIS"])
well.measured <- apply(data.mat2[,permissive.cols],1,median) > 10

logfc <- log(apply(data.rel[,selective.cols],1,mean) / apply(data.rel[,permissive.cols],1,mean))
# hist(logfc,breaks=50,col="darkolivegreen3",border="gray40")

all.data <- data.frame(
	row.names=rownames(data.rel),
	mut=clone.table[values(clone.idx,rownames(data.rel)),"aa.calls"],
	lfc=logfc,
	well.measured=well.measured,
	stringsAsFactors=FALSE
)
all.data[longdel.clones,"mut"] <- "longdel"
all.data[null.clones,"mut"] <- "null"
all.data$single <- regexpr("null|longdel|,",all.data$mut) < 0


sac.data <- all.data[all.data$mut != "longdel" & all.data$well.measured,]
lfcs <- do.call(rbind,tapply(sac.data$lfc,sac.data$mut,function(x) {
	if (length(x) > 1) t(combn(x,2)) else NULL
},simplify=FALSE))

plot(lfcs,xlim=c(-5,5),ylim=c(-5,5),pch=".",main="all well measured clones",col="steelblue3")
text(0,4,paste("R =",signif(cor(lfcs)[1,2],3)))
