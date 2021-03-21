\name{BEAST.null.simu}
\alias{BEAST.null.simu}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
BEAST Null Distribution
}
\description{
\code{BEAST.null.simu} gives a vector of the null distribution of the BEAST statistic.
}
\usage{
BEAST.null.simu(
  n,
  p,
  d,
  subsample.percent = 1/2,
  B = 100,
  lambda = NULL,
  test.independence = FALSE,
  index = NULL,
  test.uniformity = TRUE,
  p.value.method = "p",
  num.permutations = NULL
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{n}{sample size.}

\item{p}{dimension.}

\item{d}{depth of the binary expansion.}

\item{subsample.percent}{sample size for subsampling.}

\item{B}{times of subsampling.}

\item{lambda}{soft thresholding.}

\item{test.independence}{
  logicals. If \code{TRUE}, test independence among the variables listed in \code{index}. Default to be \code{FALSE}
}
\item{index}{
  a list of indices. If provided, test the independence among two or more groups of variables, for example, \eqn{(X_1, X_2)} and \eqn{X_3}.
}
\item{test.uniformity}{
  logicals. If \code{TRUE}, test mutual independence among all variables.
}
\item{p.value.method}{
  If \code{"p"}, then compute null distribution with permutations. If \code{"s"}, then compute null distribution with simulations. Otherwise, do not provide null distribution or p-value.
}
\item{num.permutations}{
  number of permutations if p.value.method == "p" (default to be 100), or simulations if p.value.method == "s" (default to be 1000).
}
}

\value{
\code{BEAST.null.simu} returns a vector of length \code{num.permutations} that simulates the null distribution of the BEAST for given sample size n, dimension p, and depth D
}

\examples{
## Elapsed times 2.42499 secs
## Measured in R 4.0.2, 32 bit, on a processor 3.3 GHz 6-Core Intel Core i5 under MacOS.
\dontrun{BEAST.null.simu(128, 2, 3)}
}

