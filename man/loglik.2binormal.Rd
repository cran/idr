\name{loglik.2binormal}
\alias{loglik.2binormal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Compute log-likelihood of parameterized bivariate 2-component Gaussian mixture models
}
\description{
  Compute the log-likelihood for parameterized bivariate 2-component
Gaussian mixture models with (1-p)N(0, 0, 1, 1, 0) + pN(mu, mu, sigma,
sigma, rho).
}
\usage{
loglik.2binormal(z.1, z.2, mu, sigma, rho, p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z.1}{ a numerical data vector on coordinate 1. }
  \item{z.2}{ a numerical data vector on coordinate 1. }
  \item{mu}{ mean for the reproducible component. }
  \item{sigma}{ standard deviation of the reproducible component. }
  \item{rho}{ correlation coefficient of the reproducible component. }
  \item{p}{ mixing proportion of the reproducible component. }  
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Log-likelihood of the bivariate 2-component
Gaussian mixture models (1-p)N(0, 0, 1, 1, 0) + N(mu, mu, sigma,
sigma, rho)$.
}
\references{
  Q. Li, J. B. Brown, H. Huang and P. J. Bickel. (2011)
   Measuring reproducibility of high-throughput experiments. Annals of Applied Statistics, Vol. 5, No. 3, 1752-1779.
}
\author{
Qunhua Li
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
    \code{m.step.2normal}, \code{e.step.normal}, \code{est.IDR}
}
\examples{
z.1 <- c(rnorm(500, 0, 1), rnorm(500, 3, 1))
rho <- 0.8

## The component with higher values is correlated with correlation coefficient=0.8 
z.2 <- c(rnorm(500, 0, 1), rnorm(500, 3 + 0.8*(z.1[501:1000]-3), (1-rho^2)))

## Starting values
mu <- 3
sigma <- 1
rho <- 0.85
p <- 0.55

## The function is currently defined as
loglik <- loglik.2binormal(z.1, z.2, mu, sigma, rho, p) 

loglik
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
