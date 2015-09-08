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
result.dir <- "ube2i_compl_sec_2015-08-19_10-46-33/UBE2I-Complementation_R1/"
# result.dir <- "ube2i_y2h_2015-08-25_08-40-30/edgotyping_R1/"

muxtag.db <- "res/muxtags"

sample.table <- read.delim(paste(muxtag.db,"_samples.tsv",sep=""),stringsAsFactors=FALSE)

data <- read.stream(result.dir)
data <- data[c(as.character(1:18),"invalid","undetermined")]
# data <- data[c(as.character(1:72),"invalid","undetermined")]

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

write.table(data.mat,"raw_counts_comp.csv",sep=",")

# data.mat <- read.csv("raw_counts.csv")
# colnames(data.mat) <- c(as.character(1:18),"invalid","undetermined")

primary <- which(substr(clones,1,1)=="U")
secondary <- which(substr(clones,1,1)=="s")
data.mat2 <- data.mat[primary,-c(19,20)]
# data.mat2 <- data.mat[-1,1:(ncol(data.mat)-2)]

# data.mat2 <- data.mat[-1,-c(19,20)]
data.rel <- apply(data.mat2,2,function(x)x/sum(x))

# data.rel <- read.csv("data_rel")

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
# for (i in seq(1,16,3)) {
# 	data.rel[,i:(i+2)] <- quantile_normalisation(data.rel[,i:(i+2)])
# }


#BOXPLOT

plotcols <- do.call(c,(lapply(do.call(c,lapply(c("firebrick","darkolivegreen","steelblue"),function(x) paste(x,c(1,4),sep=""))),rep,3)))
boxplot(data.rel[-1,-c(19,20)],ylim=c(0,.0003),col=plotcols,ylab="rel.freq.")
mtext(c("GFP","wt UBE2I","no plasmid"),side=1,line=3,at=c(3.5,9.5,15.5))
mtext(c("25C","37C","25C","37C","25C","37C"),side=1,line=2,at=c(2,5,8,11,14,17))

#CORRELATION SCATTERPLOTS

op <- par(mfrow=c(3,6))
apply(unique(sample.table[,3:4]),1,function(set) {
	idxs <- with(sample.table, which(plasmid==set[[1]] & temp==set[[2]]))
	apply(combn(idxs,2),2,function(is) {
		plot(
			data.rel[,is],
			log="xy",
			pch=".",
			col="steelblue3",
			main=paste(set,collapse=" "),
			xlab=paste("rep.",is[[1]]),
			ylab=paste("rep.",is[[2]])
		)
		text(500,500,paste("R =",signif(cor(data.mat[,is])[1,2],digits=4)))
	})
})
par(op)




# posCtrl25 <- apply(data.mat[,with(sample.table,which(plasmid=="wt UBE2I" & temp==25))],1,mean)
# posCtrl37 <- apply(data.mat[,with(sample.table,which(plasmid=="wt UBE2I" & temp==37))],1,mean)
# norm25 <- apply(data.mat[,with(sample.table,which(plasmid=="none" & temp==25))],2,`/`,posCtrl25)
# norm37 <- apply(data.mat[,with(sample.table,which(plasmid=="none" & temp==37))],2,`/`,posCtrl37)
# logfc <- log(apply(norm37,1,mean)) - log(apply(norm25,1,mean))
# p.values <- sapply(1:nrow(norm37),function(i) {
# 	tryCatch (
# 		t.test(norm37[i,],norm25[i,],alternative="two.sided")$p.value
# 		# wilcox.test(norm37[i,],norm25[i,],alternative="two.sided")$p.value
# 		,error=function(x)NA

# 	)
# })
# p.adj <- p.adjust(p.values,"fdr")


#Remove clones that were not well-measured
missing <- which(apply(data.mat2[,with(sample.table,plasmid=="none" & temp==25)],1,median) <= 10)

data.mat2 <- data.mat2[-missing,]
data.rel <- data.rel[-missing,]

c25 <- data.rel[,with(sample.table,which(plasmid=="none" & temp==25))]
c37 <- data.rel[,with(sample.table,which(plasmid=="none" & temp==37))]

#Standard-deviation of technical replicates determines well-measuredness
well.measured <- apply(c25,1,sd) < 0.000025

lfc.rep <- log(c37) - log(c25)
lfc.rep.fin <- do.call(rbind,lapply(1:nrow(lfc.rep), function(i) {
	if (!any(is.na(lfc.rep[i,])) && all(is.finite(lfc.rep[i,]))) lfc.rep[i,] else NULL
}))
op <- par(mfrow=c(1,3))
apply(combn(1:3,2),2,function(rep.is) {
	plot(
		lfc.rep[,rep.is],
		# log="xy",
		pch=".",
		col="steelblue3",
		xlab=paste("log(fc) rep.",rep.is[[1]]),
		ylab=paste("log(fc) rep.",rep.is[[2]]),
		xlim=c(-6,6),
		ylim=c(-6,6)
	)
	text(0,2,paste("R =",signif(cor(lfc.rep.fin[,rep.is])[1,2],digits=2)))
})
par(op)


logfc <- log(apply(c37,1,mean)) - log(apply(c25,1,mean))
p.values <- sapply(1:nrow(c37),function(i) {
	tryCatch (
		t.test(c37[i,],c25[i,],alternative="two.sided")$p.value, error=function(x)NA
	)
})
p.adj <- p.adjust(p.values,"fdr")

hits <- which(p.adj < .05)
hitnames <- rownames(data.rel)[hits]

pdf("volcano.pdf",10,8)
layout(rbind(1,2,3),heights=c(.2,.13,1))
op <- par(mar=c(0,4,4,1)+.1)
hist(logfc,
	main="No plasmid",axes=FALSE,col="orange",
	breaks=seq(-5,5,length.out=50),
	ylab="Freq.(all)",ylim=c(0,1100)
)
grid(NA,NULL)
axis(2)
par(op)
op <- par(mar=c(0,4,0,1)+.1)
hist(logfc[hits],main="",axes=FALSE,col="firebrick3",
	breaks=seq(-5,5,length.out=50),
	ylab="Freq.(signif.)",ylim=c(0,1100)
)
grid(NA,NULL)
axis(2)
par(op)
op <- par(mar=c(5,4,0,1)+.1)
plot(
	logfc,-log(p.adj),
	xlab="log(fold change)",
	ylab=expression(-log(p[BH])),
	pch=20,
	col="steelblue3",
	xlim=c(-5,5)
)
grid(NULL,NA)
abline(h=-log(0.05),col=2)
text(-4,-log(0.05),"p=0.05",pos=3,cex=.7,col=2)
# text(logfc[hits],-log(p.adj[hits]),rownames(data.rel)[hits],pos=1,cex=.7)
par(op)
dev.off()



# c25 <- data.rel[,with(sample.table,which(plasmid=="wt UBE2I" & temp==25))]
# c37 <- data.rel[,with(sample.table,which(plasmid=="wt UBE2I" & temp==37))]
# logfc <- log(apply(c37,1,mean)) - log(apply(c25,1,mean))
# p.values <- sapply(1:nrow(c37),function(i) {
# 	tryCatch (
# 		t.test(c37[i,],c25[i,],alternative="two.sided")$p.value, error=function(x)NA
# 	)
# })
# p.adj <- p.adjust(p.values,"fdr")

# hits <- which(p.adj < .05)

# layout(rbind(1,2),heights=c(.3,1))
# op <- par(mar=c(0,4,4,1)+.1)
# hist(logfc,main="wt UBE2I",axes=FALSE,col="orange",breaks=50)
# grid(NA,NULL)
# axis(2)
# par(op)
# op <- par(mar=c(5,4,0,1)+.1)
# plot(
# 	logfc,-log(p.adj),
# 	xlab="log(fold change)",
# 	ylab=expression(-log(p[adj])),
# 	pch=20,
# 	col="steelblue3"
# )
# grid(NULL,NA)
# abline(h=-log(0.05),col=2)
# text(-4,-log(0.05),"p=0.05",pos=3,cex=.7,col=2)
# text(logfc[hits],-log(p.adj[hits]),rownames(data.rel)[hits],pos=1,cex=.7)
# par(op)



as.df <- function(x) {
	df <- as.data.frame(lapply(1:length(x[[1]]),function(i) sapply(x,`[[`,i) ), stringsAsFactors=FALSE)
	colnames(df) <- names(x[[1]])
	df
}


clone.table <- read.csv("res/clones.csv",stringsAsFactors=FALSE)
clone.idx <- hash(clone.table$id,1:nrow(clone.table))
muts <- strsplit(clone.table$aa.calls,",")

delpos <- do.call(rbind,lapply(strsplit(clone.table$deletions,"-"),function(x) if (length(x)==1&&is.na(x)) c(NA,NA) else as.numeric(x)))
foo <- delpos[,1] < 82 & delpos[,2] > 545
foo[is.na(foo)] <- FALSE
null.clones <- clone.table$id[foo]
foo <- delpos[,1] < delpos[,2]
foo[is.na(foo)] <- FALSE
longdel.clones <- clone.table$id[foo]



all.data <- as.df(lapply(1:nrow(data.rel), function(i) {
	id <- rownames(data.rel)[[i]]
	cid <- clone.idx[[id]]
	null <- id %in% null.clones
	longdel <- id %in% longdel.clones
	list(
		id = id,
		mut = if (null) "null" else if (longdel) "longdel" else clone.table$aa.calls[[cid]],
		lfc = logfc[[i]],
		single = !null & !longdel & regexpr(",",clone.table$aa.calls[[cid]]) < 0,
		signif = id %in% hitnames,
		well.measured = well.measured[[i]]
	)
}))
rownames(all.data) <- all.data$id

big.insig <- rownames(all.data)[abs(all.data$lfc) > .5 & !all.data$signif]
sd25 <- apply(c25,1,sd)
sd37 <- apply(c37,1,sd)
op <- par(mfrow=c(4,1))
hist(log(sd25[big.insig]),breaks=seq(-20,0,.5))
hist(log(sd25[hitnames]),breaks=seq(-20,0,.5))
hist(log(sd37[big.insig]),breaks=seq(-20,0,.5))
hist(log(sd37[hitnames]),breaks=seq(-20,0,.5))
par(op)

sp <- sqrt(apply(c25,1,var)/3 + apply(c37,1,var)/3)
op <- par(mfrow=c(2,1))
hist(log(sp[big.insig]),breaks=seq(-20,0,.5))
hist(log(sp[hitnames]),breaks=seq(-20,0,.5))
par(op)

# all.data[null.clones,"mut"] <- "null"
# all.data[longdel.clones,"mut"] <- "longdel"

# hit.data <- as.df(lapply(hits, function(i) {
# 	id <- rownames(data.rel)[[i]]
# 	cid <- clone.idx[[id]]
# 	list(
# 		id = id,
# 		mut = clone.table$aa.calls[[cid]],
# 		lfc = logfc[[i]]
# 	)
# }))
# hit.data <- hit.data[order(hit.data$lfc),]

library("Biostrings")
data(BLOSUM62)

# min.blosum <- sapply(strsplit(hit.data$mut,","),function(ms) {
# 	min(sapply(ms, function(m) {
# 		from <- substr(m,1,1)
# 		to <- substr(m,nchar(m),nchar(m))
# 		BLOSUM62[from,to]
# 	}))
# })
all.data$min.blosum <- sapply(1:nrow(all.data), function(i) {
	# if (all.data$single[[i]]) {
		ms <- all.data$mut[[i]]
		if (!(ms %in% c("null","longdel"))) {
			min(sapply(ms, function(m) {
				from <- substr(m,1,1)
				to <- substr(m,nchar(m),nchar(m))
				BLOSUM62[from,to]
			}))
		} else NA
	# } else NA
})

# plot(all.data$lfc,all.data$min.blosum)

# smooth.blosum <- sapply(1:length(min.blosum),function(i) {
# 	w <- 20
# 	n <- length(min.blosum)
# 	left <- if (i < 1+w) 1 else i - w
# 	right <- if (i > n-w) n else i + w
# 	mean(min.blosum[left:right])
# })

# plot(smooth.blosum,type="l",ylab="min(BLOSUM62)",xlab="LFC Rank")

#Do positive LFC clones have hither BLOSUM scores?
wilcox.test(
	all.data$min.blosum[which(all.data$lfc >= 0)],
	all.data$min.blosum[which(all.data$lfc < 0)],
	alternative="greater"
)

single.data <- all.data[all.data$single,]
with(single.data,boxplot(tapply(lfc,min.blosum,c),xlab="BLOSUM",ylab="log(fc)"))
# single.data <- single.data[order(single.data$single.blosum,single.data$lfc),]



pp.in <- read.delim("res/polyphen_ube2i.txt",comment.char="#",stringsAsFactors=FALSE,header=FALSE)
polyphen <- with(pp.in,data.frame(row.names=paste(V8,V7,V9,sep=""),pp=V12))
single.data$pp <- polyphen[single.data$mut,"pp"]

single.data <- single.data[order(single.data$pp,single.data$lfc),]
# plot(single.data$pp,single.data$lfc,xlab="Polyphen",ylab="log(fc)")
pp.classes <- factor(sapply(single.data$pp,
	function(p) if (p < .1) "pp < 0.1" else if (p > .9) "pp > 0.9" else "0.1 < pp < 0.9"
),levels=c("pp < 0.1","0.1 < pp < 0.9","pp > 0.9"))
boxplot(tapply(single.data$lfc,pp.classes,c),ylab="log(fold-change)")



# #Are non-significant hits enriched for low barcode counts in permissive condition?
# perm.barcount <- data.mat2[,with(sample.table,plasmid=="none",temp=25)]
# sig.perm.bc <- apply(perm.barcount[hits,],1,median)
# nonsig.perm.bc <- apply(perm.barcount[-hits,],1,median)
# op <- par(mfrow=c(2,1))
# hist(sig.perm.bc,breaks=0:40000,xlim=c(0,200),ylim=c(0,300),)
# hist(nonsig.perm.bc,breaks=0:40000,xlim=c(0,200),ylim=c(0,300))
# par(op)


#Correlation between clones with same mutations
sac.data <- all.data[all.data$mut != "longdel" & all.data$well.measured,]
lfcs <- do.call(rbind,tapply(sac.data$lfc,sac.data$mut,function(x) {
	if (length(x) > 1) t(combn(x,2)) else NULL
},simplify=FALSE))
# uv.is <- lapply(unique(clone.table$aa.calls),function(uv)which(clone.table$aa.calls==uv))
# repl.is <- uv.is[which(sapply(uv.is,length) > 1)]

# lfcs <- do.call(rbind,lapply(repl.is, function(is) {
# 	cids <- clone.table$id[is]

# 	cids <- cids[cids %in% hitnames]
# 	if (length(cids) < 2) return(NULL)

# 	lfc <- logfc[cids]
# 	t(combn(lfc,2))
# }))
# colnames(lfcs) <- c("log(fc) A","log(fc) B")
plot(lfcs,xlim=c(-5,5),ylim=c(-5,5),pch=20,
	main="all signif. clones",col="steelblue3",
	xlab="log(fc) rep 1",ylab="log(fc) rep 2"
)
# .lfcs <- do.call(rbind,lapply(1:nrow(lfcs), function(i) {
# 	if (any(is.na(lfcs[i,])) || !all(is.finite(lfcs[i,]))) NULL else lfcs[i,]
# }))
text(0,4,paste("R =",signif(cor(lfcs)[1,2],3)))

#standard deviations of biological replicate sets
biorep.sd <- do.call(c,tapply(sac.data$lfc,sac.data$mut,function(x) { 
	if (length(x) > 1) sd(x) else NULL
},simplify=FALSE))


wt.ish.clones <- with(single.data,id[pp < 0.1 & min.blosum > 2 & lfc > 0])

null.lfc <- all.data[null.clones,"lfc"]
wt.lfc <- all.data[wt.ish.clones,"lfc"]

pdf("volcano.pdf",10,8)
layout(rbind(1,2,3),heights=c(.2,.13,1))
op <- par(mar=c(0,4,4,1)+.1)
hist(logfc,
	main="No plasmid",axes=FALSE,col="orange",
	breaks=seq(-5,5,length.out=50),
	ylab="Freq.(all)",ylim=c(0,1100)
)
grid(NA,NULL)
axis(2)
par(op)
op <- par(mar=c(0,4,0,1)+.1)
hist(logfc[hits],main="",axes=FALSE,col="firebrick3",
	breaks=seq(-5,5,length.out=50),
	ylab="Freq.(signif.)",ylim=c(0,1100)
)
grid(NA,NULL)
axis(2)
par(op)
op <- par(mar=c(5,4,0,1)+.1)
plot(
	logfc,-log(p.adj),
	xlab="log(fold change)",
	ylab=expression(-log(p[BH])),
	pch=20,
	col="steelblue3",
	xlim=c(-5,5)
)
grid(NULL,NA)
abline(h=-log(0.05),col=2)
text(-4,-log(0.05),"p=0.05",pos=3,cex=.7,col=2)
abline(v=wt.lfc,col="darkolivegreen3",lty="dashed")
abline(v=null.lfc,col="firebrick3",lty="dashed")
# text(logfc[hits],-log(p.adj[hits]),rownames(data.rel)[hits],pos=1,cex=.7)
par(op)
dev.off()

.lfc <- all.data$lfc
.lfc[!all.data$signif] <- 0

#Correlations between single-mut cones with same mutations
# single.muts <- unique(clone.table$aa.calls[sapply(muts,length)==1])
# uv.is <- lapply(single.muts,function(uv)which(clone.table$aa.calls==uv))
# repl.is <- uv.is[which(sapply(uv.is,length) > 1)]

# repl.cids <- lapply(repl.is,function(is) clone.table$id[is])

# lfcs <- do.call(rbind,lapply(repl.is, function(is) {
# 	# is <- is[clone.table$freq[is] > .6]
# 	# if (length(is) < 2) return(NULL)
# 	cids <- clone.table$id[is]

# 	cids <- cids[cids %in% hitnames]
# 	if (length(cids) < 2) return(NULL)

# 	lfc <- logfc[cids]
# 	t(combn(lfc,2))
# }))
# colnames(lfcs) <- c("log(fc) A","log(fc) B")
# plot(lfcs,xlim=c(-5,5),ylim=c(-5,5),pch=20,main="signif. single mutants",col="steelblue3")

# .lfcs <- do.call(rbind,lapply(1:nrow(lfcs), function(i) {
# 	if (any(is.na(lfcs[i,])) || !all(is.finite(lfcs[i,]))) NULL else lfcs[i,]
# }))
# text(0,4,paste("R =",signif(cor(.lfcs)[1,2],3)))

#Generate an index of single mutant rows
muts <- strsplit(all.data$mut,",")
single.idx <- hash()
for (i in 1:nrow(all.data)) {
	if (all.data$single[[i]] && all.data$well.measured[[i]]) {
		single.idx[[muts[[i]]]] <- c(single.idx[[muts[[i]]]],i)
	}
}

#sum of single LFC (product of FC)
doubleVsingles <- do.call(rbind,lapply(1:nrow(all.data), function(i) {
	ms <- muts[[i]]
	if (length(ms)==2 && all.data$well.measured[[i]]) {
		lfc <- all.data[i,"lfc"]
		pred <- sum(sapply(ms, function(m) {
			idxs <- single.idx[[m]]
			if (is.null(idxs)) return(NA)
			mean(all.data[idxs,"lfc"])
		}))
		if (!is.na(pred)) {
			c(real=lfc,predicted=pred)
		} else NULL
	} else NULL
}))
plot(doubleVsingles,xlim=c(-4,4),ylim=c(-4,4),col="steelblue3",pch=20)
abline(h=0,v=0,lty="dashed",col="gray")
text(0,4,paste("R =",signif(cor(doubleVsingles)[1,2],3)))


#Predict singles from doubles
to.df <- function(x) {
	out <- do.call(data.frame,lapply(1:ncol(x),function(i)unlist(x[,i])))
	colnames(out) <- colnames(x)
	out
}
plfcs <- to.df(do.call(rbind,lapply(1:nrow(all.data), function(i) {
	ms <- muts[[i]]
	if (length(ms)==2 && all.data$well.measured[[i]]) {
		dlfc <- all.data[i,"lfc"]
		do.call(rbind,lapply(ms, function(m) {
			other.m <- setdiff(ms,m)
			idxs <- single.idx[[m]]
			if (is.null(idxs)) return(NULL)
			slfc <- mean(all.data[idxs,"lfc"])
			list(mut=other.m,plfc=dlfc-slfc)
		}))
	} else NULL
})))
plfcs <- tapply(plfcs$plfc,plfcs$mut,mean)

singleVsPred <- na.omit(do.call(rbind,lapply(keys(single.idx),function(sm) {
	c(
		real=mean(all.data[single.idx[[sm]],"lfc"]),
		pred=if(sm %in% names(plfcs)) plfcs[[sm]] else NA
	)
})))
plot(singleVsPred,xlim=c(-4,4),ylim=c(-4,4),col="steelblue3",pch=20)
abline(h=0,v=0,lty="dashed",col="gray")
text(0,4,paste("R =",signif(cor(singleVsPred)[1,2],3)))


#make row index over mutations occurring anywhere in a clone
anymut.idx <- hash()
for (i in 1:length(muts)) {
	ms <- muts[[i]]
	for (m in ms) {
		anymut.idx[[m]] <- c(anymut.idx[[m]],i)
	}
}
#predict singles from averages
plfcs2 <- sapply(keys(anymut.idx), function(m) {
	is <- anymut.idx[[m]]
	if (length(is) >=3 ) mean(all.data[is,"lfc"]) else NA
})
singleVsPred2 <- do.call(rbind,lapply(keys(single.idx), function(sm) {
	# idxs <- which(sapply(muts,function(ms)sm%in%ms) & !all.data$single)
	idxs <- anymut.idx[[sm]]
	idxs <- idxs[!all.data$single[idxs]]
	if (length(idxs)==0) return(NULL)
	c(
		real=mean(all.data[single.idx[[sm]],"lfc"]),
		pred=mean(all.data[idxs,"lfc"])
	)
}))
plot(singleVsPred2,xlim=c(-4,4),ylim=c(-4,4),col="steelblue3",pch=20)
abline(h=0,v=0,lty="dashed",col="gray")
text(0,4,paste("R =",signif(cor(singleVsPred2)[1,2],3)))

table(sapply(keys(anymut.idx),function(key)length(anymut.idx[[key]])))

all.combos <- function(x) {
	if (length(x) == 0) {
		return(list(NULL))
	}
	if (length(x) == 1) {
		return(list(x))
	}
	do.call(c,lapply(1:length(x), function(n) {
		do.call(c,apply(combn(length(x),n),2,function(is) {
			a <- x[is]
			b <- x[-is]
			lapply(all.combos(b),function(cs)c(list(a),cs))
		}))
	}))
}



#Plot all doubles versus singles
dvs <- na.omit(do.call(rbind,lapply(1:nrow(hit.data),function(i) {
	ms <- hit.muts[[i]]
	if (length(ms) == 2) {
		lfc <- hit.data[i,"lfc"]
		pred <- sapply(ms, function(m) {
			idxs <- hit.single.idx[[m]]
			if (is.null(idxs)) return(NA)
			mean(hit.data[idxs,"lfc"])
		})
		if (any(is.na(pred))) return(NULL) else c(lfc,sort(pred))
	} else return(NULL)
})))
z <- lm(dvs[,1]~dvs[,-1])
plot(dvs[,1],predict(z),xlim=c(-3,5),ylim=c(-3,5),pch=20,col="firebrick3")
cor(dvs[,1],predict(z))

plot(0,type="n",xlim=c(0,nrow(dvs)),ylim=c(-3,5),xlab="double mutant rank",ylab="log(fold-change)")
rect(1:nrow(dvs)-1,dvs[,1]-.05,1:nrow(dvs),dvs[,1]+.05,col="firebrick3",border=NA)
rect(1:nrow(dvs)-1,dvs[,2]-.05,1:nrow(dvs),dvs[,2]+.05,col="royalblue3",border=NA)
rect(1:nrow(dvs)-1,dvs[,3]-.05,1:nrow(dvs),dvs[,3]+.05,col="darkolivegreen3",border=NA)
abline(h=0,col="gray")
legend("topleft",c("double mutant","lower single","higher single"),fill=c("firebrick3","royalblue3","darkolivegreen3"))



# mut.idx <- hash()
# for (i in 1:length(hit.muts)) {
# 	if (length(hit.muts[[i]]) > 1) {
# 		for (m in hit.muts[[i]]) {
# 			mut.idx[[m]] <- c(mut.idx[[m]],i)
# 		}
# 	}
# }
# sva <- do.call(rbind,lapply(keys(hit.single.idx),function(sm) {
# 	sng <- mean(hit.data[hit.single.idx[[sm]],"lfc"])
# 	avg <- mean(hit.data[mut.idx[[sm]],"lfc"])
# 	c(single=sng,average=avg)
# }))
# plot(sva,xlim=c(-3,5),ylim=c(-3,5),pch=20,col="firebrick3")


ube2i.prot <- scan("res/ube2i_aa.fa",what="character")[[2]]
wt.aa <- c(sapply(1:nchar(ube2i.prot),function(i)substr(ube2i.prot,i,i)),"*")
# insig.singles <- setdiff(single.muts,keys(hit.single.idx))



aas <- c("A","V","L","I","M","F","Y","W","R","H","K","D","E","S","T","N","Q","G","C","P")
hmap <- matrix(NA,nrow=length(aas),ncol=160,dimnames=list(aas,1:160))
#Fill in predictions
for (i in 1:length(plfcs2)) {
	m <- names(plfcs2)[[i]]
	if (m %in% c("longdel","null")) next
	lfc <- plfcs2[[i]]
	to.aa <- substr(m,nchar(m),nchar(m))
	pos <- as.integer(substr(m,2,nchar(m)-1))
	hmap[to.aa,pos] <- lfc
}
#Override product-rule predictions
for (i in 1:length(plfcs)) {
	m <- names(plfcs)[[i]]
	lfc <- plfcs[[i]]
	to.aa <- substr(m,nchar(m),nchar(m))
	pos <- as.integer(substr(m,2,nchar(m)-1))
	hmap[to.aa,pos] <- lfc
}
#Override in single mutants and wt
for (m in keys(single.idx)) {
	# from.aa <- substr(m,1,1)
	to.aa <- substr(m,nchar(m),nchar(m))
	pos <- as.integer(substr(m,2,nchar(m)-1))
	lfc <- mean(all.data[single.idx[[m]],"lfc"])
	hmap[to.aa,pos] <- lfc
}
# for (m in insig.singles) {
# 	to.aa <- substr(m,nchar(m),nchar(m))
# 	pos <- as.integer(substr(m,2,nchar(m)-1))
# 	hmap[to.aa,pos] <- NaN
# }

pdf("comp_map_inferred.pdf",width=16,height=5)
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
	else  cp[[round(v+5)]]
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

