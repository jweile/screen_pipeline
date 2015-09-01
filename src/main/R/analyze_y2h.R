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

# data.mat <- read.csv("raw_counts.csv")
# colnames(data.mat) <- c(as.character(1:72),"invalid","undetermined")

# primary <- which(substr(clones,1,1)=="U")
# secondary <- which(substr(clones,1,1)=="s")
# data.mat2 <- data.mat[primary,-c(19,20)]
data.mat2 <- data.mat[-1,1:(ncol(data.mat)-2)]

# data.mat2 <- data.mat[-1,-c(19,20)]
data.rel <- apply(data.mat2,2,function(x)x/sum(x))

# data.rel <- read.csv("data_rel")

#BOXPLOT
ord <- with(sample.table,order(interactor,condition,replicate))
plotcols <- do.call(c,(lapply(do.call(c,lapply(c("firebrick","darkolivegreen","steelblue","peachpuff","mediumpurple","aquamarine","darkgoldenrod","lightsteelblue"),function(x) paste(x,c(1,3,4),sep=""))),rep,3)))
boxplot(data.rel[,ord],ylim=c(0,.0003),col=plotcols,ylab="rel.freq.")
mtext(unique(sample.table$interactor[ord]),side=1,line=3,at=0:7*9+4.5)
mtext(c("-HIS","+HIS","+3AT"),side=1,line=2,at=seq(2,71,3))
