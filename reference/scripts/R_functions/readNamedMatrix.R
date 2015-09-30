readNamedMatrix <- function(matFile, what=double()) {
	if (!file.exists(matFile)) {
		stop(paste("Cannot find file '", matFile, "'", sep=""))
	}

	#
	#m = read.table(matFile, header=TRUE, check.names=FALSE, row.names=1)
	#return(m)
	#

        use.matFile = matFile
        # Escape file name just in case:
        if (system("which printf", ignore.stdout=TRUE, ignore.stderr=TRUE) == 0 && system("which sh", ignore.stdout=TRUE, ignore.stderr=TRUE) == 0 && system("echo 'printf %q \"a\"' | sh", ignore.stdout=TRUE, ignore.stderr=TRUE) == 0) {
            use.matFile = system(paste("echo 'printf %q \"", matFile, "\"' | sh", sep=""), intern=TRUE)
        }

        READ_MAT_FILE_PREFIX = ""
        READ_MAT_FILE_SUFFIX = paste(" ", use.matFile, sep="")
        
        if (regexpr("\\.gz$", matFile) > 0) {
            READ_MAT_FILE_PREFIX = paste("gzip -cd ", use.matFile, " | ", sep="")
            READ_MAT_FILE_SUFFIX = ""
        }

	############################
	# AS IN XHMM C++ CODE:
	# Instead of splitting by whitespace, use ONLY tab as delimiter (to allow for sample names with space in them):
	############################


	getDims = paste(READ_MAT_FILE_PREFIX, "awk -F'\t' '{print NF}'", READ_MAT_FILE_SUFFIX, " | sort | uniq -c | awk '{print $1,$2}'", sep="")
	dimsVec = system(getDims, intern=TRUE)
	if (length(dimsVec) != 1) {
		stop(paste("Cannot read jagged matrix: ", matFile, sep=""))
	}
	# Subtract 1 for row and column names:
	rows_cols = as.numeric(strsplit(dimsVec, "\\s+")[[1]]) - 1
	rows = rows_cols[1]
	cols = rows_cols[2]

	writeLines(paste("Reading ", rows, " x ", cols, " named matrix", sep=""))

	n = rows * cols
	readMat = paste(READ_MAT_FILE_PREFIX, "awk -F'\t' 'BEGIN{OFS=\"\t\"} {$1=\"\"; gsub(\"^\"OFS\"+\", \"\", $0); print $0;}'", READ_MAT_FILE_SUFFIX, sep="")

	if (log2(n) < 31) {
		con <- pipe(readMat)
		m = matrix(scan(con, what=what, n=n, skip=1, quiet=TRUE, sep="\t"), rows, cols, byrow = TRUE)
		close(con)
	}
	else {
		m = data.frame()

		for (r in 1:rows) {
			con <- pipe(readMat)
			m = rbind(m, scan(con, what=what, n=cols, skip=(1+(r-1)), quiet=TRUE, sep="\t"))
			close(con)
		}
	}

	readColNames = paste(READ_MAT_FILE_PREFIX, "awk -F'\t' 'BEGIN{OFS=\"\t\"} {if (NR == 1) {$1=\"\"; gsub(\"^\"OFS\"+\", \"\", $0); print $0;", ifelse(READ_MAT_FILE_PREFIX == "", "exit", ""), "}}'", READ_MAT_FILE_SUFFIX, sep="")
	con <- pipe(readColNames)
	colnames(m) = scan(con, what="", quiet=TRUE, sep="\t")
	close(con)

	readRowNames = paste(READ_MAT_FILE_PREFIX, "awk -F'\t' 'BEGIN{OFS=\"\t\"} {if (NR > 1) print $1}'", READ_MAT_FILE_SUFFIX, sep="")
	con <- pipe(readRowNames)
	rownames(m) = scan(con, what="", quiet=TRUE, sep="\t")
	close(con)

	return(m)
}
