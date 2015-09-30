\name{loadTargetsToGenes}
\alias{loadTargetsToGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Load the exome-target to gene mappings.
}
\description{
Load the exome-target to gene mappings.
}
\usage{
loadTargetsToGenes(targToGenesFile, targRegexpReplace =
list(c("\\\\.\\\\.", "-")), geneDelim = "\\\\|",
targColumn = 1, geneColumn = 3, CHROMOSOMES_START_WITH_CHR = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{targToGenesFile}{
%%     ~~Describe \code{targToGenesFile} here~~
}
  \item{targRegexpReplace}{
%%     ~~Describe \code{targRegexpReplace} here~~
}
  \item{geneDelim}{
%%     ~~Describe \code{geneDelim} here~~
}
  \item{targColumn}{
%%     ~~Describe \code{targColumn} here~~
}
  \item{geneColumn}{
%%     ~~Describe \code{geneColumn} here~~
}
  \item{CHROMOSOMES_START_WITH_CHR}{
%%     ~~Describe \code{CHROMOSOMES_START_WITH_CHR} here~~
}
}
\examples{
\dontrun{targsToGenes = loadTargetsToGenes("targs_to_genes.txt")}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
