readTableNoRowNames <- function(tableFile, what=character(), setNumeric=NULL) {
	if (!file.exists(tableFile)) {
		stop(paste("Cannot find file '", tableFile, "'", sep=""))
	}

	#
	#m = read.table(tableFile, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
	#return(m)
	#

        READ_TABLE_FILE_PREFIX = ""
        READ_TABLE_FILE_SUFFIX = paste(" ", tableFile, sep="")
        
        if (regexpr("\\.gz$", tableFile) > 0) {
            READ_TABLE_FILE_PREFIX = paste("gzip -cd ", tableFile, " | ", sep="")
            READ_TABLE_FILE_SUFFIX = ""
        }

	############################
	# AS IN XHMM C++ CODE:
	# Instead of splitting by whitespace, use ONLY tab as delimiter (to allow for sample names with space in them):
	############################


	getDims = paste(READ_TABLE_FILE_PREFIX, "awk -F'\t' '{print NF}'", READ_TABLE_FILE_SUFFIX, " | sort | uniq -c | awk '{print $1,$2}'", sep="")
	dimsVec = system(getDims, intern=TRUE)
	if (length(dimsVec) != 1) {
		stop(paste("Cannot read jagged table: ", tableFile, sep=""))
	}
	# Subtract 1 for column names:
	rows_cols = as.numeric(strsplit(dimsVec, "\\s+")[[1]])
	rows = rows_cols[1] - 1
	cols = rows_cols[2]

	writeLines(paste("Reading ", rows, " x ", cols, " table", sep=""))

	n = rows * cols
	readMat = paste(READ_TABLE_FILE_PREFIX, "awk -F'\t' 'BEGIN{OFS=\"\t\"} {print $_}'", READ_TABLE_FILE_SUFFIX, sep="")

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

	readColNames = paste(READ_TABLE_FILE_PREFIX, "awk -F'\t' 'BEGIN{OFS=\"\t\"} {if (NR == 1) {print $_; ", ifelse(READ_TABLE_FILE_PREFIX == "", "exit", ""), "}}'", READ_TABLE_FILE_SUFFIX, sep="")
	con <- pipe(readColNames)
	colnames(m) = scan(con, what="", quiet=TRUE, sep="\t")
	close(con)

        m = as.data.frame(m, check.names=FALSE, stringsAsFactors=FALSE)
        if (!is.null(setNumeric)) {
            for (col in setNumeric) {
                if (col %in% colnames(m)) {
                    m[, col] = as.numeric(m[, col])
                }
            }
        }

	return(m)
}
