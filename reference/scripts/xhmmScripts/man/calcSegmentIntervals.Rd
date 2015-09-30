\name{calcSegmentIntervals}
\alias{calcSegmentIntervals}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Determine all of the indices where a transition in value occurs.
}
\description{
Determine all of the vector indices where a transition in value occurs
along the input vector.
}
\usage{
calcSegmentIntervals(valuesVectorList)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{valuesVectorList}{
%%     ~~Describe \code{valuesVectorList} here~~
}
}
\examples{
ints = c(rep(1, 100), rep(2, 200), rep(3, 300))
breaks = calcSegmentIntervals(list(ints))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
