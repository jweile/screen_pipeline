#!$Rbin

####
# demuxer.R is a worker script that runs on the cluster nodes
# It uses k-mer search to identify multiplexing tags (muxtags)
# which are then used in pairs to identify samples. It further
# uses k-mer search to identify clones according to barcodes 
# which are split into groups according to their samples and then 
# counted
# 
# Written by Jochen Weile <jochenweile@gmail.com>
# library("Biostrings")

# source("../src/main/R/liblogging.R")   #Logger
# source("../src/main/R/cliargs.R")      #Command-line argument processing
# source("../src/main/R/libyogitools.R") #Helper functions
# source("../src/main/R/libyogiseq.R")   #FASTQ and bowtie

source("lib/liblogging.R")   #Logger
source("lib/cliargs.R")      #Command-line argument processing
source("lib/libyogitools.R") #Helper functions
source("lib/libyogiseq.R")   #FASTQ and bowtie

# R1 reads chunk file.
r1.file <- getArg("r1",required=TRUE)
# R2 reads chunk file.
r2.file <- getArg("r2",required=TRUE)
# Working directory
dir.name <- getArg("dir",required=TRUE)
# Job ID
job.id <- getArg("id",required=TRUE)
# Location of well tag DB
muxtag.db <- getArg("muxtags",default="res/muxtags")
# A csv file containing information about clones and their barcodes
clone.db <- getArg("cloneDB",default="res/clones")

# Turns on debug mode
debug.mode <- as.logical(getArg("debug",default=FALSE))

#Create logger
log.file <- paste(dir.name,"demuxer_",job.id,".log",sep="")
logger <- new.logger(log.file)

#Function for loading sequences
read.fastq <- function(f) {
	tryCatch({
		con <- gzfile(f, open="r")
		p <- new.fastq.parser(con)
		out <- list()
		while (length(s <- p$parse.next(1)) > 0) {
			out[[length(out)+1]] <- s[[1]]
		}
		out
	},
	error = function(ex) {
		logger$fatal(paste("Error reading file",f," :\n",ex))
		stop(ex)
	},
	finally = {
		if (exists("con") && isOpen(con)) {
			close(con)
		}
	})
}

#Function that returns the names of the sequences
seqnames <- function(seqs) sapply(seqs,function(s)s$getID())

#Load reads
logger$info("Loading sequence data...")
r1.seq <- read.fastq(r1.file)
r2.seq <- read.fastq(r2.file)

if (!all(seqnames(r2.seq)==seqnames(r1.seq))) {
	logger$fatal("R1 and R2 reads do not correspond to each other!")
	stop()
}

#####
# STEP 1: identify muxtags
#####
logger$info("Identifying multiplexing tags...")

#Get muxtag combos from sample table
sample.table <- read.delim(paste(muxtag.db,"_samples.tsv",sep=""),stringsAsFactors=FALSE)
mux.combos <- with(sample.table,paste(tag1,tag2,sep="-"))

#K-mer search implemented in libyogiseq.R
ks <- new.kmer.search()
kmer.file <- paste(muxtag.db,"_index.rdata",sep="")
ks$load.index(kmer.file)
#Trunkate sequence information around the muxtag position plusminus 2bp before k-mer search
r1.mux.ids <- ks$search(substr(sapply(r1.seq,function(x)x$toString()),4,16),useAlignment=FALSE)
r2.mux.ids <- ks$search(substr(sapply(r2.seq,function(x)x$toString()),4,16),useAlignment=FALSE)
r12.mux.combos <- paste(r1.mux.ids,r2.mux.ids,sep="-")

logger$info("Calling samples...")
#Call samples according to muxtag combinations
sample.ids <- sapply(r12.mux.combos,function(x) {
	if (regexpr("NA",x) > 0) {
		"undetermined"
	} else {
		hits <- mux.combos == x
		if (any(hits)) {
			which(hits)
		} else {
			"invalid"
		}
	}
})

#####
# STEP 2: Extract Barcodes
#####
logger$info("Identifying barcodes and clones...")

ks <- new.kmer.search()
ks$load.index(paste(clone.db,"_index.rdata",sep=""))
#Trunkate sequence information around the bc position plusminus 2bp before k-mer search
# clone.ids <- ks$search(substr(sapply(r2.seq,function(x)x$toString()),32,60),max.d=7)
clone.ids <- ks$search(substr(sapply(r2.seq,function(x)x$toString()),32,60),useAlignment=FALSE,min.hits=8)

#####
# STEP 3: GROUP ACCORDING TO SAMPLES AND COUNT
#####
logger$info("Grouping by sample and counting...")

out <- tapply(clone.ids, sample.ids, function(cids) {
	table(cids,useNA="ifany")
},simplify=FALSE)

#Function for safely writing sequences
write.fastq <- function(f,seqs,gz=TRUE) {
	tryCatch({
		con <- if (gz) gzfile(f,open="w") else file(f, open="w")
		writeFASTQ(con,seqs)
	},
	error = function(ex) {
		logger$fatal(paste("Error while writing file",f," :\n",ex))
		stop(ex)
	},
	finally = {
		if (exists("con") && isOpen(con)) {
			close(con)
		}
	})
}

#Sort reads in to subdirectories according to samples
invisible(tapply(1:length(r1.seq), sample.ids, function(idx) {
	sample.id <- sample.ids[idx[[1]]]
	sub.dir <- paste(dir.name,sample.id,"/",sep="")
	if (!file.exists(sub.dir)) dir.create(sub.dir,showWarnings=FALSE)

	r1.out.file <- paste(sub.dir,"R1_",job.id,".fastq.gz",sep="")
	r2.out.file <- paste(sub.dir,"R2_",job.id,".fastq.gz",sep="")
	write.fastq(r1.out.file,r1.seq[idx])
	write.fastq(r2.out.file,r2.seq[idx])
}))

logger$info("Writing output...")

out.file <- paste(dir.name,"counts_",job.id,".txt.gz",sep="")
con <- gzfile(out.file,open="w")
invisible(lapply(1:length(out),function(i) {
	sample.id <- names(out)[[i]]
	writeLines(paste("#",sample.id,sep=""),con)
	writeLines(paste(names(out[[i]]),out[[i]],sep=":"),con)
}))
close(con)


#####
# CLEANUP: Delete input files when we're done, to save HD space.
####
if (!debug.mode) {
	logger$info("Cleaning up temporary files...")
	file.remove(r1.file)
	file.remove(r2.file)
}

logger$info("Done!")
