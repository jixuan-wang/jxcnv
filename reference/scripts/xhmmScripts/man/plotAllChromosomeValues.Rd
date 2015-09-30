\name{plotAllChromosomeValues}
\alias{plotAllChromosomeValues}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot to a png file the values for the given intervals.
}
\description{
Plot to a png file the values for the given intervals, alternating color
by chromosome.
}
\usage{
plotAllChromosomeValues(chr, values, ylab, title,
outputName, specialColors = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chr}{
%%     ~~Describe \code{chr} here~~
}
  \item{values}{
%%     ~~Describe \code{values} here~~
}
  \item{ylab}{
%%     ~~Describe \code{ylab} here~~
}
  \item{title}{
%%     ~~Describe \code{title} here~~
}
  \item{outputName}{
%%     ~~Describe \code{outputName} here~~
}
  \item{specialColors}{
%%     ~~Describe \code{specialColors} here~~
}
}
\examples{
\dontrun{plotAllChromosomeValues(c(rep(1, 200), rep(2, 100), rep(3, 50)),
c(sample(1:200), sample(1:100), sample(1:50)), "Values", "Genomic values", "./out")}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
