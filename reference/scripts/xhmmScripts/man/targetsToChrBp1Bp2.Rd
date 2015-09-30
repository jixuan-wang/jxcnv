\name{targetsToChrBp1Bp2}
\alias{targetsToChrBp1Bp2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parse a list of genomic intervals into their component chromosome, start, and stop.
}
\description{
Parse a list of genomic intervals into their component chromosome, start, and stop.
}
\usage{
targetsToChrBp1Bp2(targets, startStopDelim = "-")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{targets}{
%%     ~~Describe \code{targets} here~~
}
  \item{startStopDelim}{
%%     ~~Describe \code{startStopDelim} here~~
}
}
\examples{
targetsToChrBp1Bp2(c("1:23-29", "X:100-200", "Y:10-50"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
