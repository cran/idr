 \name{idr-package}
\alias{idr-package}
\alias{idr}
\docType{package}
\title{
Irreproducible discovery rate
}
\description{
This package estimates the reproducibility of observations on
a pair of replicate rank lists. It consists of three components:
(1) plotting the correspondence curve to visualize reproducibility, (2)
quantifying reproducibility using a copula mixture model and estimating
the posterior probability for each obsrvation to be irreproducible
(local irreproducible discovery rate), and (3) ranking and selecting
observations by  their irreproducibility.
}
\details{
\tabular{ll}{
Package: \tab idr\cr
Type: \tab Package\cr
Version: \tab 1.2\cr
Date: \tab 2012-10-26\cr
Updates: \tab Improve the convergence of est.IDR (2014-08-15) \cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
The main functions are est.IDR(), get.correspondence() and
select.IDR(). est.IDR estimates the copula mixture model and the
posterior probability for each observation to be irreproducible. get.correspondence generates the values for plotting the correspondence curve. select.IDR ranks obervations by their reproducibility and reports
the number of observations passing the specified IDR thresholds.
%~~ An overview of how to use the package, including the most important ~~
%~~ functions ~~
}
\author{
Qunhua Li

Maintainer: Qunhua Li <qunhua.li@gmail.com>
%~~ The author and/or maintainer of the package ~~
}
\references{
   Q. Li, J. B. Brown, H. Huang and P. J. Bickel. (2011)
 Measuring reproducibility of high-throughput experiments. Annals of Applied Statistics, Vol. 5, No. 3, 1752-1779.


%~~ Literature or other references for background information ~~
}
\keyword{ package }
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
\examples{

data("simu.idr")
x <- cbind(-simu.idr$x, -simu.idr$y)

mu <- 2.6
sigma <- 1.3
rho <- 0.8
p <- 0.7

idr.out <- est.IDR(x, mu, sigma, rho, p, eps=0.001, max.ite=20)

names(idr.out)
}
