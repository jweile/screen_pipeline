library("hash")
options(stringsAsFactors=FALSE)

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

write.table(data.mat,"raw_counts_y2h.csv",sep=",")
# data.mat <- read.csv("raw_counts.csv")

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

interactors <- unique(sample.table$interactor)

pdf("techrep_y2h.pdf",12,4)
all.data <- do.call(cbind,lapply(interactors,function(ia) {

	permissive.cols <- with(sample.table,sample[interactor==ia & condition=="+HIS"])
	selective.cols <- with(sample.table,sample[interactor==ia & condition=="-HIS"])
	well.measured <- apply(data.mat2[,permissive.cols],1,median) > 10

	logfc <- log(apply(data.rel[,selective.cols],1,mean) / apply(data.rel[,permissive.cols],1,mean))

	logfc.rep <- log(data.rel[,selective.cols] / data.rel[,permissive.cols])
	logfc.rep <- do.call(rbind,lapply(1:nrow(logfc.rep),function(i) {
		if (any(is.na(logfc.rep[i,])) || any(is.infinite(logfc.rep[i,]))) NULL else logfc.rep[i,]
	}))

	op <- par(mfrow=c(1,3))
	invisible(apply(combn(3,2),2,function(is) {
		plot(logfc.rep[,is],
			xlim=c(-4,4),ylim=c(-4,4),pch=".",col="steelblue3",
			xlab=paste("lfc rep",is[[1]]),ylab=paste("lfc rep",is[[2]]),
			main=paste(ia,is[[1]],"vs",is[[2]])
		)
		text(0,4,paste("R =",signif(cor(na.omit(logfc.rep[,is]))[1,2],3)))
	}))
	par(op)
	# hist(logfc,breaks=50,col="darkolivegreen3",border="gray40")

	out <- data.frame(lfc=logfc,wm=well.measured)
	colnames(out) <- paste(ia,c("lfc","wm"),sep=".")
	out

}))
dev.off()

# all.data <- data.frame(
# 	row.names=rownames(data.rel),
# 	mut=clone.table[values(clone.idx,rownames(data.rel)),"aa.calls"],
# 	lfc=logfc,
# 	well.measured=well.measured,
# 	stringsAsFactors=FALSE
# )

all.data$mut <- clone.table[values(clone.idx,rownames(data.rel)),"aa.calls"]
all.data[longdel.clones,"mut"] <- "longdel"
all.data[null.clones,"mut"] <- "null"
all.data$single <- regexpr("null|longdel|,",all.data$mut) < 0

pdf("biorep_y2h.pdf",16,8)
op <- par(mfrow=c(2,4))
invisible(lapply(interactors,function(ia) {
	
	sac.data <- all.data[all.data$mut != "longdel" & all.data[,paste(ia,"wm",sep=".")],]

	#lfc pairs for biological replicates
	br.pairs <- do.call(rbind,tapply(sac.data[,paste(ia,"lfc",sep=".")],sac.data$mut,function(x) {
		if (length(x) > 1) t(combn(x,2)) else NULL
	},simplify=FALSE))
	plot(br.pairs,xlim=c(-5,5),ylim=c(-5,5),pch=".",main=ia,col="steelblue3")
	text(0,4,paste("R =",signif(cor(br.pairs)[1,2],3)))
}))
par(op)
dev.off()



muts <- strsplit(all.data$mut,",")
single.idx <- hash()
for (i in 1:nrow(all.data)) {
	if (all.data$single[[i]]) {
		single.idx[[muts[[i]]]] <- c(single.idx[[muts[[i]]]],i)
	}
}



ube2i.prot <- scan("res/ube2i_aa.fa",what="character")[[2]]
wt.aa <- c(sapply(1:nchar(ube2i.prot),function(i)substr(ube2i.prot,i,i)),"*")
# insig.singles <- setdiff(single.muts,keys(hit.single.idx))


aas <- c("A","V","L","I","M","F","Y","W","R","H","K","D","E","S","T","N","Q","G","C","P")
hmap <- matrix(NA,nrow=length(aas),ncol=160,dimnames=list(aas,1:160))

#fill in single mutants and wt
for (m in keys(single.idx)) {
	# from.aa <- substr(m,1,1)
	to.aa <- substr(m,nchar(m),nchar(m))
	pos <- as.integer(substr(m,2,nchar(m)-1))
	lfc <- mean(all.data[single.idx[[m]],"ZBED1.lfc"])
	hmap[to.aa,pos] <- lfc
}


pdf("comp_map_zbed1.pdf",width=16,height=5)
layout(cbind(c(3,1),c(4,2)),widths=c(9.5,.5),heights=c(2,9))
op <- par(cex=.6,las=1,mar=c(5,4,0,0)+.1)
plot(
	0,type="n",
	xlim=c(-0.5,ncol(hmap)),ylim=c(0,nrow(hmap)),
	axes=FALSE,xlab="AA position",ylab="Amino acid"
)
axis(1,at=c(1,seq(10,160,10))-.5,labels=c(1,seq(10,160,10)))
axis(2,at=20:1-.5,labels=aas)
text(-1,c(16,7.5),c("hydrophobic","polar"),srt=90)
# text(-3,c(8,10.5),c("\U2296","\U2295"),srt=90)
text(-3,c(8,10.5),c("-","+"),srt=90)
arrows(-2,c(3.1,12.1),-2,c(11.9,19.9),length=.02,angle=90,code=3)
arrows(-4,c(7.1,9.1),-4,c(8.9,11.9),length=.02,angle=90,code=3)
xy <- expand.grid(x=1:ncol(hmap),y=1:nrow(hmap))
cp <- colorRampPalette(c("royalblue3","white","firebrick3"))(9)
col.fill <- apply(xy,1,function(.xy) {
	v <- hmap[.xy[["y"]],.xy[["x"]]]
	wt <- wt.aa[[.xy[["x"]]]] == aas[[.xy[["y"]]]]
	if (wt) "lightgoldenrod1" 
	else if (is.nan(v)) "white"
	else if (is.na(v)) "gray80" 
	else {
		i <- round(v+5)
		if (i < 1) i <- 1
		if (i > 9) i <- 9
		cp[[i]]
	}
})
rect(xy$x-1,20-xy$y,xy$x,21-xy$y,col=col.fill,border=NA)
# cross out insignificant values
# for (m in insig.singles) {
# 	to.aa <- substr(m,nchar(m),nchar(m))
# 	x <- as.integer(substr(m,2,nchar(m)-1))
# 	y <- 21-which(aas==to.aa)
# 	segments(c(x-1,x-1),c(y-1,y),c(x,x),c(y,y-1),col="gray")
# }
segments(0,c(3,7,9,12),160,lty="dotted")
#legend
par(mar=c(5,0,0,4)+.1)
plot(0,type="n",xlim=c(0,1),ylim=c(0,11),xlab="",ylab="",axes=FALSE)
rect(0,0:8,1,1:9,col=cp,border=NA)
rect(0,9,1,10,col="lightgoldenrod1",border=NA)
rect(0,10,1,11,col="gray80",border=NA)
# rect(0,11,1,12,col="white",border=NA)
# segments(c(0,0),c(11,12),c(1,1),c(12,11),col="gray")
axis(4,at=0:10+.5,labels=c(seq(-4,4),"wt","n/d"),tick=FALSE)
mtext("log(fold-change)",side=4,line=2,las=3,cex=.6)

par(mar=c(0,4,1,0)+.1)
# barvals <- apply(hmap,2,mean,na.rm=TRUE)
barvals <- apply(apply(hmap,2,function(x)c(
	sum(x >= 0.5,na.rm=TRUE),
	sum(x < 0.5 & x > -0.5,na.rm=TRUE)+sum(is.nan(x)),
	sum(x <= -0.5,na.rm=TRUE)
)),2,function(x)x/sum(x))
plot(
	0,type="n",
	xlim=c(-0.5,ncol(hmap)),
	# ylim=c(-2.2,2.2),
	ylim=c(0,1),
	axes=FALSE,xlab="",
	# ylab="mean log(fc)"
	ylab="pos/neutral/neg"
)
# rect(1:ncol(hmap)-1,0,1:ncol(hmap),barvals,col=cp[round(barvals+5)],border="gray")
rect(1:ncol(hmap)-1,0,1:ncol(hmap),barvals[3,],col=cp[[2]],border="gray")
rect(1:ncol(hmap)-1,barvals[3,],1:ncol(hmap),barvals[3,]+barvals[2,],col="white",border="gray")
rect(1:ncol(hmap)-1,barvals[3,]+barvals[2,],1:ncol(hmap),1,col=cp[[8]],border="gray")
grid(NA,NULL)
axis(2)
par(op)
dev.off()

