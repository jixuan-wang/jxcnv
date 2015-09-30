\name{plot_XHMM_targets}
\alias{plot_XHMM_targets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the XHMM CNV spanning the input target indices.
}
\description{
Plot the XHMM CNV spanning the input target indices.
}
\usage{
plot_XHMM_targets(BASE_OUTPUT_NAME, xhmm_data, allTargsToGenes,
binarySampleFeatures, SQ_THRESH, startTarg, stopTarg,
PLOT_ONLY_PNG = TRUE, PLOT_LINE_SEGMENTS = TRUE, DEFAULT_WIDTH = 14,
DEFAULT_HEIGHT = 7, UPDATE_FIG_HEIGHT_FOR_MARGINS = TRUE, GENE_NAME_cex
= 0.8,
PLOT_PIPELINE_TRANSITIONS = TRUE, DIFFERENTIATE_SAMPLE_GROUPS = NULL,
COLOR_DIFFERENTIAL_TARGETS = TRUE, MARK_SAMPLES = NULL,
MARK_SAMPLES_COLORS = DEFAULT_MARK_SAMPLES_COLORS,
MARK_SAMPLES_ONLY_IN_CALLS = FALSE,
EXCLUDE_SAMPLES = NULL, MARK_INTERVALS = NULL, APPEND_REGION_NAME = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BASE_OUTPUT_NAME}{
%%     ~~Describe \code{BASE_OUTPUT_NAME} here~~
}
  \item{xhmm_data}{
%%     ~~Describe \code{xhmm_data} here~~
}
  \item{allTargsToGenes}{
%%     ~~Describe \code{allTargsToGenes} here~~
}
  \item{binarySampleFeatures}{
%%     ~~Describe \code{binarySampleFeatures} here~~
}
  \item{SQ_THRESH}{
%%     ~~Describe \code{SQ_THRESH} here~~
}
  \item{startTarg}{
%%     ~~Describe \code{startTarg} here~~
}
  \item{stopTarg}{
%%     ~~Describe \code{stopTarg} here~~
}
  \item{PLOT_ONLY_PNG}{
%%     ~~Describe \code{PLOT_ONLY_PNG} here~~
}
  \item{PLOT_LINE_SEGMENTS}{
%%     ~~Describe \code{PLOT_LINE_SEGMENTS} here~~
}
  \item{DEFAULT_WIDTH}{
%%     ~~Describe \code{DEFAULT_WIDTH} here~~
}
  \item{DEFAULT_HEIGHT}{
%%     ~~Describe \code{DEFAULT_HEIGHT} here~~
}
  \item{UPDATE_FIG_HEIGHT_FOR_MARGINS}{
%%     ~~Describe \code{UPDATE_FIG_HEIGHT_FOR_MARGINS} here~~
}
  \item{GENE_NAME_cex}{
%%     ~~Describe \code{GENE_NAME_cex} here~~
}
  \item{PLOT_PIPELINE_TRANSITIONS}{
%%     ~~Describe \code{PLOT_PIPELINE_TRANSITIONS} here~~
}
  \item{DIFFERENTIATE_SAMPLE_GROUPS}{
%%     ~~Describe \code{DIFFERENTIATE_SAMPLE_GROUPS} here~~
}
  \item{COLOR_DIFFERENTIAL_TARGETS}{
%%     ~~Describe \code{COLOR_DIFFERENTIAL_TARGETS} here~~
}
  \item{MARK_SAMPLES}{
%%     ~~Describe \code{MARK_SAMPLES} here~~
}
  \item{MARK_SAMPLES_COLORS}{
%%     ~~Describe \code{MARK_SAMPLES_COLORS} here~~
}
  \item{MARK_SAMPLES_ONLY_IN_CALLS}{
%%     ~~Describe \code{MARK_SAMPLES_ONLY_IN_CALLS} here~~
}
  \item{EXCLUDE_SAMPLES}{
%%     ~~Describe \code{EXCLUDE_SAMPLES} here~~
}
  \item{MARK_INTERVALS}{
%%     ~~Describe \code{MARK_INTERVALS} here~~
}
  \item{APPEND_REGION_NAME}{
%%     ~~Describe \code{APPEND_REGION_NAME} here~~
}
}
\examples{
\dontrun{plot_XHMM_targets("./out", xhmm_data, allTargsToGenes,
binarySampleFeatures, SQ_THRESH, 150, 180)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
