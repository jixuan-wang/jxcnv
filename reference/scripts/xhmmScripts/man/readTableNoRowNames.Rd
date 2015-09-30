\name{readTableNoRowNames}
\alias{readTableNoRowNames}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Quickly read a rectangular table that has column names.
}
\description{
Quickly read a rectangular table that has column names.
}
\usage{
readTableNoRowNames(tableFile, what = character(), setNumeric = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tableFile}{
If \code{tableFile} ends in '.gz', then this file will be read as a gzip-ed text file.
}
  \item{what}{
%%     ~~Describe \code{what} here~~
}
  \item{setNumeric}{
%%     ~~Describe \code{what} here~~
}
}
\examples{
\dontrun{readTableNoRowNames("tableWithHeader.txt")}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
