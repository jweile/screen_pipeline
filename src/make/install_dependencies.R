#!/home/jweile/bin/Rscript


list.of.packages <- c("parallel","hash","XML","bitops","RJSONIO","RCurl")
missing <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if (length(missing) > 0) {
	install.packages(missing,repos="http://cran.rstudio.com/")
}

list.of.packages <- c("Biostrings")
missing <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if (length(missing) > 0) {
	source("http://bioconductor.org/biocLite.R")
	biocLite("Biostrings")
}
