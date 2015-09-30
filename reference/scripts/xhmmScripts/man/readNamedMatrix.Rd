\name{readNamedMatrix}
\alias{readNamedMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Quickly read a rectangular matrix that has row and column names.
}
\description{
Quickly read a rectangular matrix that has row and column names.
}
\usage{
readNamedMatrix(matFile, what = double())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matFile}{
If \code{matFile} ends in '.gz', then this file will be read as a gzip-ed text file.
}
  \item{what}{
%%     ~~Describe \code{what} here~~
}
}
\examples{
\dontrun{readNamedMatrix("table.txt")}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
