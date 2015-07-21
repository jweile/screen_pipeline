####
## THIS LIBRARY IS A COLLECTION OF VARIOUS USEFUL TOOLS I WROTE
#


###
# This object can count occurrences of different items
#
new.counter <- function() {

	a <- list()

	###
	# Add x occurrences to item id
	#
	add <- function(id,x) {
		if (!is.character(id)) {
			stop("Illegal argument:",id)
		}
		if (is.null(a[[id]])) {
			a[[id]] <<- x
		} else {
			a[[id]] <<- a[[id]] + x
		}
	}

	# increase counter for item id by 1
	inc <- function(id) add(id,1)

	# get the counter state for id
	get <- function(id) a[[id]]

	# list counts for all ids
	ls <- function() a

	# export counter state as a string
	export <- function() {
		paste(lapply(names(a), function(id) paste(id,"=",a[[id]],sep="") ), collapse=",")
	}

	# import counter state from string
	import.add <- function(strs) { 
		lapply(strsplit(strs,","), function(eqs) {
			lapply(strsplit(eqs,"="), function(vals) {
				add(vals[[1]],as.numeric(vals[[2]]))
			})
		})
		invisible()
	}

	structure(
		list(
			inc = inc,
			add = add,
			get = get,
			ls = ls,
			export = export,
			import.add = import.add
		),
		class="yogicounter"
	)
}

###
# Function to *locally* excise regex groups from string vectors.
# I.e. only extract the first occurrence of each group within each string.
# x = string vector
# re = regular expression with groups
#
extract.groups <- function(x, re) {
	matches <- regexpr(re,x,perl=TRUE)
	start <- attr(matches,"capture.start")
	end <- start + attr(matches,"capture.length") - 1
	do.call(cbind,lapply(1:ncol(start), function(i) {
		sapply(1:nrow(start),function(j){
			if (start[j,i] > -1) substr(x[[j]],start[j,i],end[j,i]) else NA
		})
	}))
}
###
# Function to *globally* excise regex groups from string vectors.
# x = string vector
# re = regular expression with groups
#
global.extract.groups <- function(x,re) {
    all.matches <- gregexpr(re,x,perl=TRUE)
    mapply(function(matches,x) {
        start <- attr(matches,"capture.start")
        end <- start + attr(matches,"capture.length") - 1
        apply(zbind(start,end),c(1,2),function(pos) substr(x,pos[[1]],pos[[2]]) )
    },matches=all.matches,x=x,SIMPLIFY=FALSE)
}




#Function for returning the i'th ranked item in a list
ith.rank <- function(values, i) sort(values,decreasing=TRUE)[i]

###
# Matthew's correlation coefficient (MCC)
# 
mcc <- function(t, scores, truth) {

	# exclude <- is.na(scores) | is.na(truth)
	# scores <- scores[-exclude]
	# truth <- truth[-exclude]

	.truth <- truth == 1
	.calls <- scores >= t

	tp <- sum(.truth & .calls)
	tn <- sum(!.truth & !.calls)
	fp <- sum(.calls & !.truth)
	fn <- sum(.truth & !.calls)

	# mcc <- (tp * tn - fp * fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
	#this formula prevents integer overflow errors
	mcc <- exp( log(tp*tn - fp*fn) - ( log(tp+fp) + log(tp+fn) + log(tn+fp) + log(tn+fn) )/2 )
	prec <- tp/(tp+fp)
	recall <- tp/(tp+fn)
	c(mcc=mcc,prec=prec,recall=recall)
}



# ###
# # assemble list of matrices into 3d matrix
# #
# matrix.3d <- function(matrices) {
# 	nr <- nrow(matrices[[1]])
# 	nc <- ncol(matrices[[1]])
# 	nm <- length(matrices)
# 	mat <- rep(0,nr*nc*nm)
# 	dim(mat) <- c(nr,nc,nm)
# 	rownames(mat) <- rownames(matrices[[1]])
# 	colnames(mat) <- colnames(matrices[[1]])
# 	for(i in 1:nm) {
# 		mat[,,i] <- matrices[[i]]
# 	}
# 	mat
# }


###
# turn rbound lists in to a dataframe
#
to.df <- function(x) {
	if (is.null(x)) return(NULL)
	y <- lapply(1:ncol(x), function(col) {
		unlist(x[,col])
	})
	names(y) <- colnames(x)
	as.data.frame(y,stringsAsFactors=FALSE)
}

###
# Binds matrices of same size together to a 3D matrix, analogously
# to cbind and rbind.
#
zbind <- function(...) {
	x <- list(...)
	y <- array(0,dim=c(nrow(x[[1]]),ncol(x[[1]]),length(x)),dimnames=dimnames(x[[1]]))
	for (i in 1:length(x)) y[,,i] <- x[[i]]
	y
}

####
# converting between plate coordinate systems
#
q2c <- function(x) {
	q <- which(LETTERS==substr(x,1,1))
	r <- which(LETTERS==substr(x,3,3))
	c <- as.numeric(substr(x,4,nchar(x)))
	.r <- r*2 - (q<3)
	.c <- c*2 - q%%2
	paste(LETTERS[[.r]],sprintf("%02d",.c),sep="")
}
c2q <- function(x) {
	r <- which(LETTERS==substr(x,1,1))
	c <- as.numeric(substr(x,2,nchar(x)))
	.r <- ceiling(r/2)
	.c <- ceiling(c/2)
	.q <- ((r-1)%%2)*2 + (c-1)%%2+1
	paste(LETTERS[[.q]],"_",LETTERS[[.r]],sprintf("%02d",.c),sep="")
}


###
# This function forms all possible subsets of a given set
#
combo <- function(l) {
	do.call(c,lapply(1:length(l),function(n){
		tab <- combn(l,n)
		lapply(1:ncol(tab),function(i)tab[,i])
	}))
}


###
# This object can be used to create cluster maps
#
new.cluster.map <- function(n) {
	
	.clusters <- as.list(1:n)

	.getIdx <- function(i) which(sapply(.clusters,function(x) i %in% x))

	addLink <- function(i,j) {
		i.idx <- .getIdx(i)
		j.idx <- .getIdx(j)
		joint <- union(.clusters[[i.idx]],.clusters[[j.idx]])
		.clusters[c(i.idx,j.idx)] <<- NULL
		.clusters[[length(.clusters)+1]] <<- joint
	}

	getClusters <- function() .clusters

	list(addLink=addLink, getClusters=getClusters, getIdxOf=.getIdx)
}
