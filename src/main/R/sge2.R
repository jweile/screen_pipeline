
library("XML")

###
# This class manages the submission of jobs to the 
# SunGridEngine cluster
#
new.sge <- function(max.queue.length=20, logger=NULL, debug=FALSE) {

	.debug <- debug
	#queue of jobs yet to submit
	.queue <- list()
	#jobs already submitted
	.submitted <- list()
	#the maximal permitted queue length
	.ql <- max.queue.length

	#returns a list of the IDs of all previously submitted jobs
	job.ids <- function() union(
		sapply(.submitted,"[[",1),
		sapply(.queue,"[[",1)
	)

	running <- function() {
		#query the queue state in XML format and parse
		qstat <- pipe("qstat -xml",open="r")
		xml <- xmlParse(paste(readLines(qstat),collapse="\n"))
		close(qstat)
		#extract job ids of queued jobs
		all.running <- do.call(c, xpathApply(xml, "//JB_name", xmlValue))
		#intersection of queued jobs with jobs submitted via this interface
		intersect(job.ids(), all.running)
	}

	submit <- function() {
		#while queue is shorter than maximum, submit more of them
		
		while ((nrunning <- length(running())) < .ql && length(.queue) > 0) {
			# cat(nrunning,"jobs running.\n")
			job <- .queue[[1]]
			#add job to job list
			.submitted[[length(.submitted)+1]] <<- job
			#submit job to sge
			submission.cmd <- paste(
				"qsub -V -N",job$id,
				ifelse(!.debug,"-e /dev/null",""),
				"-o /dev/null",
				"-wd `pwd` -b y",job$command,
				paste(job$arguments,collapse=" ")
			)
			# cat(submission.cmd,"\n")
			system(submission.cmd)
			#remove job from queue
			.queue[[1]] <<- NULL
		}
		
	}

	enqueue <- function(id, command, arguments) {
		#check if job exists already
		if (id %in% job.ids()) {
			stop("Job ID already exists")
		}

		#enqueue
		.queue[[length(.queue)+1]] <<- list(id=id,command=command,arguments=arguments)

		submit()
		report()
	}

	.lastSL <- -1
	.lastQL <- -1
	report <- function(verbose=FALSE) {
		currSL <- length(running())
		currQL <- length(.queue)
		if (currSL != .lastSL || currQL != .lastQL) {
			if (!is.null(logger)) {
				logger$info(paste("Submitted:",currSL,"\tQueued:",currQL))	
			} else if (verbose) {
				cat("\rSubmitted:",currSL,"\tQueued:",currQL,"        ")
			}
			.lastSL <<- currSL
			.lastQL <<- currQL
		}
	}

	wait <- function(verbose=FALSE) {
		
		if (verbose) cat("\n")

		while (length(.queue) > 0 || length(running()) > 0) {
			submit()
			# if (verbose) cat("\rSubmitted:",length(running()),"\tQueued:",length(.queue),"        ")	
			# if (!is.null(logger)) logger$info(paste("\rSubmitted:",length(running()),"\tQueued:",length(.queue)))	
			report(verbose)
			Sys.sleep(1)
		}

		if (verbose) cat("\n")
	}

	structure(
		list(
			enqueue = enqueue,
			wait = wait,
			running = running
		),class="sge"
	)
}