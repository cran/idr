\name{m.step.2normal}
\alias{m.step.2normal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
M-step for parameterized bivariate 2-component Gaussian mixture models
}
\description{
Maximization step in the EM algorithm for parameterized bivariate 2-component
Gaussian mixture models with $(1-p)N(0, 0, 1, 1, 0) + pN(mu, mu, sigma^2,
sigma^2, rho)$. 
}
\usage{
m.step.2normal(z.1, z.2, e.z)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z.1}{ a numerical data vector on coordinate 1.}
  \item{z.2}{ a numerical data vector on coordinate 2.}
  \item{e.z}{ a vector of expected conditional probability that the $i$th
    observation is reproducible. }
}
\details{
  This function is used in the EM algorithm for estimating the parameters
of the Gaussian mixture model at the latent copula space. 
}
\value{
  Estimated parameters, basically a list including elements 
  \item{p}{ the estimated mixing proportion of the reproducible component.}
  \item{mu}{ the estimated mean for the reproducible component.}
  \item{sigma}{ the estimated standard deviation of the reproducible component.}
  \item{rho}{ the estimated correlation coefficient of the reproducible component.}
}
\references{
  Q. Li, J. B. Brown, H. Huang and P. J. Bickel. (2011)
   Measuring reproducibility of high-throughput experiments. Annals of Applied Statistics, Vol. 5, No. 3, 1752-1779.
}
\author{
Qunhua Li
}
%\note{
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{e.step.2normal}, \code{loglik.2binormal}, \code{est.IDR}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

 
z.1 <- c(rnorm(500, 0, 1), rnorm(500, 3, 1))
rho <- 0.8

##The component with higher values is correlated with correlation coefficient=0.8 
z.2 <- c(rnorm(500, 0, 1), rnorm(500, 3 + 0.8*(z.1[501:1000]-3), (1-rho^2)))
e.z <- c(rep(0, 500) + abs(rnorm(0, 0.05)), rep(1, 500)-abs(rnorm(0, 0.05)))

para <- m.step.2normal(z.1, z.2, e.z) 

para
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
