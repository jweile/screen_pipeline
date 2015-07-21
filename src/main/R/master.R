#!$Rbin

####
# master.R is the top-level controller script that chops large fastq
# files into smaller chunks on-the-fly and spawns worker jobs on the 
# SGE cluster to process them. Afterwards it collates the results for
# further processing.
# 
# Written by Jochen Weile <jochenweile@gmail.com>

library("hash",warn.conflicts=FALSE,quietly=TRUE,verbose=FALSE)

source("lib/sge2.R")         #SUN Grid Engine
source("lib/cliargs.R")      #Command-line argument processing
source("lib/liblogging.R")   #Writing log files
source("lib/libyogitools.R") #Some handy tools
source("lib/libyogiseq.R")   #k-mer search

###
# Get command line arguments
#

# Files containing the R1 and R2 reads. Example file:
r1.files <- getArg("r1",required=TRUE)
r2.files <- getArg("r2",required=TRUE)
if (length(strsplit(r1.files,",")[[1]]) != length(strsplit(r2.files,",")[[1]])) {
	stop("There needs to be an equal number of r1 and r2 files!")
}
rfile.table <- cbind(r1=strsplit(r1.files,",")[[1]],r2=strsplit(r2.files,",")[[1]])


# The session tag is used to name the ouptut directory together with a timestamp.
session.tag <- getArg("session",default="screen")

# The chunk size determines how many reads are processed by each slave script.
chunk.size <- as.integer(getArg("chunksize",default=20000))
if (is.na(chunk.size)) stop("chunksize must be integer number!")

# Turns on debug mode
debug.mode <- as.logical(getArg("debug",default=FALSE))
if (is.na(debug.mode)) stop("debug must be TRUE or FALSE !")

# This option determines how many jobs will be submitted to the cluster at maximum.
# When more jobs are available, they get buffered internally until enough room exists in the queue.
max.queue <- as.integer(getArg("maxQueue",default=30))
if (is.na(max.queue)) stop("maxQueue must be integer number!")

# Location of well tag DB
muxtag.db <- getArg("muxtags",default="res/muxtags")
if (!file.exists(paste(muxtag.db,".fa",sep=""))) stop("muxtag DB not found!")

clone.db <- getArg("cloneDB",default="res/clones")
if (!file.exists(paste(clone.db,".fa",sep=""))) stop ("clone DB not found!")


###
# Create output directory
#
timestamp <- format(Sys.time(),format='%Y-%m-%d_%H-%M-%S')
out.dir <- paste(session.tag,"_",timestamp,"/", sep="")
dir.create(out.dir, mode="0755")

# Set up log file
logger <- new.logger(paste(out.dir,"master.log",sep=""))

# Build k-mer search index
logger$info("Building indices...")
ks <- new.kmer.search()
ks$build.index(paste(muxtag.db,".fa",sep=""))
ks <- new.kmer.search()
ks$build.index(paste(clone.db,".fa",sep=""))


###
# PHASE 1: CREATE CHUNKS OF READS ON-THE-FLY AND DEMUX THEM
#

###
# This function uses an already open file connection
# to create a new FASTQ chunk file. R
#
# incon = incoming file connection
# out.dir = output directory
# i = chunk number
# direction = R1 or R2
#
make.chunk <- function(incon, out.dir, i, direction) {

	outfile <- paste(out.dir,direction,"-",i,".fastq",sep="")
	outcon <- file(outfile,open="w")

	lines <- readLines(incon, chunk.size*4)
	writeLines(lines,outcon)

	close(outcon)

	list(file=outfile, last=(length(lines) < chunk.size*4))
}

###
# Function to find out whether file is GZip file
#
is.gz <- function(f) {
	substr(f,nchar(f)-2,nchar(f))==".gz" || 
		regexpr("gzip compressed data", system(paste("file",f),intern=TRUE) ) > 0
}

###
# work through read file pairs
#
# Construct SunGridEngine object with designated maximum queue size.
sge <- new.sge(max.queue.length=max.queue, logger=logger, debug=debug.mode)

logger$info("Demultiplexing...")

#Iterate over file pairs, create chunks and submit jobs on the go.
#The return result directories.
result.dirs <- apply(rfile.table, 1, function(rfiles) {

	r1.file <- rfiles[["r1"]]
	r2.file <- rfiles[["r2"]]

	sub.dir <- gsub(".+/|\\.fastq(\\.gz)?$","",r1.file)

	#TODO: name directory according to SWIM well
	dir.name <- paste(out.dir,sub.dir,"/",sep="")
	dir.create(dir.name)

	###
	# Open connections to sequencing results files
	#
	con.r1 <- file(r1.file,open="r")
	if (is.gz(r1.file)) {
		con.r1 <- gzcon(con.r1)
	}
	con.r2 <- file(r2.file,open="r")
	if (is.gz(r2.file)) {
		con.r2 <- gzcon(con.r2)
	}


	done <- FALSE
	i <- 0
	while (!done) {
		i <- i+1
		# Make chunks for R1 and R2
		r1.chunk <- make.chunk(con.r1,dir.name,i,"R1")
		r2.chunk <- make.chunk(con.r2,dir.name,i,"R2")

		#Create a job id
		job.id <- paste(session.tag,sub.dir,timestamp,i,sep="_")


		#Submit Slave job to SunGridEngine
		sge$enqueue(
			id=job.id,
			command="$Rbin",
			arguments=list(
				"lib/demuxer.R",
				paste("r1=",r1.chunk$file,sep=""),
				paste("r2=",r2.chunk$file,sep=""),
				paste("dir=",dir.name,sep=""),
				paste("id=",job.id,sep=""),
				paste("muxtags=",muxtag.db,sep=""),
				paste("cloneDB=",clone.db,sep=""),
				paste("debug=",debug.mode,sep="")
			)
		)

		#We're done if we run out of reads to process
		done <- r1.chunk$last || r2.chunk$last
	}
	#Close file connections
	close(con.r1)
	close(con.r2)

	dir.name
})
#Wait for the remaining jobs to finish
sge$wait(verbose=TRUE)


# ####
# # PHASE 2: CONSOLIDATE JOB RESULTS
# #
logger$info("Consolidating results...")

read.stream <- function(result.dir) {
	data <- list()
	con <- gzcon(pipe(paste("cat ",result.dir,"*.txt.gz",sep=""),open="rb"))
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

logger$info("Writing results to file...")

invisible(lapply(result.dirs, function(result.dir) {

	out <- read.stream(result.dir)

	out.file <- paste(result.dir,"all_counts.txt.gz",sep="")
	con <- gzfile(out.file,open="w")
	invisible(lapply(1:length(out),function(i) {
		sample.id <- names(out)[[i]]
		writeLines(paste("#",sample.id,sep=""),con)
		keys <- keys(out[[i]])
		writeLines(paste(keys,values(out[[i]],keys),sep=":"),con)
	}))
	close(con)

	#Clean up slave output
	if (!debug.mode) {
		file.remove(list.files(dir.name,pattern="counts_.+\\.txt\\.gz"))
	}

}))


logger$info("Done!")
