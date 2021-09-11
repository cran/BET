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
  dep,
  subsample.percent = 1/2,
  B = 100,
  lambda = NULL,
  index = NULL,
  method = "p",
  num = NULL
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{n}{sample size.}

\item{p}{dimension.}

\item{dep}{depth of the binary expansion.}

\item{subsample.percent}{sample size for subsampling.}

\item{B}{times of subsampling.}

\item{lambda}{tuning parameter for soft-thresholding, default to be \eqn{sqrt(log(2^(p * dep)) / (8*n))}.}

\item{index}{
  a list of indices. If provided, test the independence among two or more groups of variables. For example, \code{index = list(c(1,2), c(3)))} refers to test the independence between \eqn{(X_1, X_2)} and \eqn{X_3}. Default to be \code{list(c(1:p))} to test if the data follow the uniform distribution over [0,1]^p, where \code{p = ncol(X)}.
}
\item{method}{
  If \code{"p"}, then compute null distribution with permutations. If \code{"s"}, then compute null distribution with simulations.
}
\item{num}{
  number of permutations if method == "p" (default to be 100), or simulations if method == "s" (default to be 1000).
}
}

\value{
\code{BEAST.null.simu} returns a vector of length \code{num.permutations} that simulates the null distribution of the BEAST for given sample size n, dimension p, and depth D
}

\examples{
## Elapsed times 2.4 secs
## Measured in R 4.0.2, 32 bit, on a processor 3.3 GHz 6-Core Intel Core i5 under MacOS, 2021/5/30.
\dontrun{BEAST.null.simu(128, 2, 3)}

## power study example
## Elapsed times 36.8 secs
## Measured in R 4.0.2, 32 bit, on a processor 3.3 GHz 6-Core Intel Core i5 under MacOS, 2021/5/30.
\dontrun{
nperm = 100
nsim = 1000
BEAST.res = rep(0, nsim)
## simulate null distribution: only need once
BEAST.null.dist = BEAST.null.simu(128, 3, 3, index = list(c(1,2), c(3)), num = nperm)
for(i in 1:nsim){
  x1 = runif(128)
  x2 = runif(128)
  y = sin(4*pi*(x1 + x2)) + 0.8*rnorm(128)
  BEAST.stat = BEAST(cbind(x1, x2, y), 3, index = list(c(1,2), c(3)),
    method = "stat")$BEAST.Statistic
  BEAST.pvalue = sum(BEAST.null.dist >= BEAST.stat) / nperm
  BEAST.res[i] = BEAST.pvalue
}
## compute power
level = 0.1
power = mean(BEAST.res < level)

}


}

