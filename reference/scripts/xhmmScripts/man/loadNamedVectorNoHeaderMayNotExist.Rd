\name{loadNamedVectorNoHeaderMayNotExist}
\alias{loadNamedVectorNoHeaderMayNotExist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Load a vector with named rows.
}
\description{
  Load a vector with named rows.
  If the file does not exist, warn the user and return an empty vector.
}
\usage{
loadNamedVectorNoHeaderMayNotExist(fileName)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fileName}{
%%     ~~Describe \code{fileName} here~~
}
}
\examples{
\dontrun{vec = loadNamedVectorNoHeaderMayNotExist("named_vector.txt")}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
