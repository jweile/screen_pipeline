###
# Create a new logger instance
#
new.logger <- function(logfile="R.log") {

	.logfile <- logfile

	if (!file.exists(.logfile)) {
		file.create(.logfile)
	}

	if (file.access(.logfile, mode=2) < 0) {
		stop("Logfile cannot be written!")
	}

	.log <- function(type,msg) {
		tryCatch({
			timestamp <- format(Sys.time(),format='%Y-%m-%d_%H:%M:%S')
			.con <- file(.logfile,open="a")
			out <- paste(timestamp,type,msg)
			cat(out,"\n")
			writeLines(out, .con)
		},
		error = function(ex) {
			print(ex)
		},
		finally = {
			if (exists(".con") && isOpen(.con)) {
				close(.con)
			}
		})
		

	}

	info <- function(msg) .log("INFO:",msg)
	warn <- function(msg) .log("WARNING:",msg)
	fatal <- function(msg) .log("FATAL:",msg)

	structure(list(
		info=info,
		warn=warn,
		fatal=fatal
	),class="logger")
}
