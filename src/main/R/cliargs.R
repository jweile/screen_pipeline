
###
# get the user-supplied argument with the given name
# if not defined return default value
getArg <- function(name, default=NULL, required=FALSE) {

	if (length(commandArgs(TRUE)) == 0) {
		if (required) {
			stop("Required argument:",name)
		} else {
			return(default)
		}
	}

	#tabulate arguments by name
	argTable <- do.call(rbind,strsplit(commandArgs(TRUE),"="))
	#get index of argument with given name
	i <- which(argTable[,1] == name)


	if (length(i) == 0) {
		#return default value if no matching arguments are found.
		if (required) {
			stop("Required argument:",name)
		} else {
			return(default)
		}
	} else if (length(i) > 1) {
		#if multiple matches are found, throw error message.
		stop("Multiple values for", name, "found!")
	} else {
		#if everything checks out, return the argument value
		return(argTable[i,2])
	}
}
