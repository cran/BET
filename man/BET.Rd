\name{BET}
\alias{BET}
\title{Binary Expansion Testing at a Certain Depth}
\usage{
BET(X, d)
}
\arguments{
\item{X}{a matrix to be tested. When \code{X} has only one column, \code{BET} will test whether \code{X} is [0, 1]-uniformly distributed (an error will be given if data is out of range [0, 1]). When \code{X} has two or more columns, \code{BET} tests the independence among those column vectors.}

\item{d}{depth of the binary expansion for the \code{BET}. }
}
\details{
\code{BET} tests the independence or uniformity by considering the maximal magnitude of the symmetry statistics in the \eqn{sigma}-field generated from marginal binary expansions at the depth \code{d}.
}
\value{
\item{Interaction}{a dataframe with \eqn{p} columns, where \eqn{p} is the number of columns of \code{X}. It displays the interactions where the extreme symmetry statistics happens. For each column in \code{X}, we use a binary index to indicate binary variables involved in the extreme symmetry statistic.}
\item{Extreme.Asymmetry}{the extreme asymmetry statistics.}
\item{p.value}{p-value of the test with Bonferroni adjustment.}
\item{z.statistic}{normal approximation of the test statistic.}
}
\description{
\code{BET} stands for Binary Expansion Testing. It is used for nonparametric detection of nonuniformity or dependence. It can be used to test whether a column vector is [0, 1]-uniformly distributed. It can also be used to detect dependence between columns of a matrix \code{X}, if \code{X} has more than one column.
}
\examples{
v <- runif(128, -pi, pi)
X1 <- cos(v) + 2.5 * rnorm(128, 0, 1/20)
X2 <- sin(v) + 2.5 * rnorm(128, 0, 1/20)
BET(cbind(X1, X2), 3)
}
