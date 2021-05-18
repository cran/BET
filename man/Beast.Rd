\name{BEAST}
\alias{BEAST}
\title{Binary Expansion Approximate Symmetry Test}
\usage{
BEAST(
  X,
  d,
  subsample.percent = 1/2,
  B = 100,
  unif.margin = FALSE,
  lambda = NULL,
  test.independence = FALSE,
  index = NULL,
  test.uniformity = TRUE,
  null.simu = NULL,
  p.value.method = "p",
  num.permutations = NULL
)
}
\arguments{
\item{X}{a matrix to be tested.}

\item{d}{depth of the binary expansion for the \code{BEAST}.}

\item{subsample.percent}{sample size for subsampling.}

\item{B}{times of subsampling.}

\item{unif.margin}{logicals. If \code{TRUE} the marginal distribution of data is known. Default to be \code{FALSE} and the function uniforms the data.}

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
\item{null.simu}{
  a vector of the null distribution of the BEAST statistic.
}
\item{p.value.method}{
  If \code{"p"}, then compute null distribution with permutations. If \code{"s"}, then compute null distribution with simulations. Otherwise, do not provide null distribution or p-value.
}
\item{num.permutations}{
  number of permutations if p.value.method == "p" (default to be 100), or simulations if p.value.method == "s" (default to be 1000).
}
}

\value{
\item{Interaction}{the most frequent interaction among all subsampling.}
\item{BEAST.Statistic}{BEAST statistic.}
\item{p.value}{simulated p-value.}
}

\description{
\code{BEAST} (Binary Expansion Approximate Symmetry Test) is used for nonparametric detection of nonuniformity or dependence.
}
\examples{
## Elapsed times 2.975584 secs
## Measured in R 4.0.2, 32 bit, on a processor 3.3 GHz 6-Core Intel Core i5 under MacOS.
\dontrun{
x1 = runif(128)
x2 = runif(128)
y = sin(4*pi*(x1 + x2)) + 0.8*rnorm(128)
BEAST(cbind(x1, x2, y), 3, test.independence = TRUE,
index = list(c(1,2), c(3)), test.uniformity = FALSE, p.value.method = "p")
}
}
