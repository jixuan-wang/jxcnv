\name{XHMM_plots}
\alias{XHMM_plots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make a series of plots documenting the output of XHMM.
}
\description{
Make a series of plots documenting the output of XHMM.
}
\usage{
XHMM_plots(PLOT_PATH, JOB_PREFICES, JOB_TARGETS_TO_GENES = NULL,
SAMPLE_FEATURES = NULL, SQ_THRESH = 60, NUM_ADD_TARGS = 2,
PLOT_READ_DEPTHS = TRUE, PLOT_PC_CORRS = TRUE,
PLOT_ALL_CNVS = TRUE, USE_XCNV_TO_PLOT = NULL,
INCLUDE_PEDIGREE_SAMPLES = NULL, PLOT_ONLY_PNG = TRUE, LIMIT_MEMORY = FALSE,
MIN_PC_SD_TO_PLOT = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PLOT_PATH}{
%%     ~~Describe \code{PLOT_PATH} here~~
}
  \item{JOB_PREFICES}{
%%     ~~Describe \code{JOB_PREFICES} here~~
}
  \item{JOB_TARGETS_TO_GENES}{
%%     ~~Describe \code{JOB_TARGETS_TO_GENES} here~~
}
  \item{SAMPLE_FEATURES}{
%%     ~~Describe \code{SAMPLE_FEATURES} here~~
}
  \item{SQ_THRESH}{
%%     ~~Describe \code{SQ_THRESH} here~~
}
  \item{NUM_ADD_TARGS}{
%%     ~~Describe \code{NUM_ADD_TARGS} here~~
}
  \item{PLOT_READ_DEPTHS}{
%%     ~~Describe \code{PLOT_READ_DEPTHS} here~~
}
  \item{PLOT_PC_CORRS}{
%%     ~~Describe \code{PLOT_PC_CORRS} here~~
}
  \item{PLOT_ALL_CNVS}{
%%     ~~Describe \code{PLOT_ALL_CNVS} here~~
}
  \item{USE_XCNV_TO_PLOT}{
%%     ~~Describe \code{USE_XCNV_TO_PLOT} here~~
}
  \item{INCLUDE_PEDIGREE_SAMPLES}{
%%     ~~Describe \code{INCLUDE_PEDIGREE_SAMPLES} here~~
}
  \item{PLOT_ONLY_PNG}{
%%     ~~Describe \code{PLOT_ONLY_PNG} here~~
}
  \item{LIMIT_MEMORY}{
%%     ~~Describe \code{LIMIT_MEMORY} here~~
}
  \item{MIN_PC_SD_TO_PLOT}{
%%     ~~Describe \code{MIN_PC_SD_TO_PLOT} here~~
}
}
\examples{
\dontrun{XHMM_plots(PLOT_PATH, JOB_PREFICES, JOB_TARGETS_TO_GENES, SAMPLE_FEATURES)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
