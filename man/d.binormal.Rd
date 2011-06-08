\name{d.binormal}
\alias{d.binormal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Log density of bivariate Gaussian distribution with symmetric marginals
}
\description{
Compute the log-density for parameterized bivariate Gaussian
distribution N(mu, mu, sigma, sigma, rho).
}
\usage{
d.binormal(z.1, z.2, mu, sigma, rho)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z.1}{ a numerical data vector on coordinate 1. }
  \item{z.2}{ a numerical data vector on coordinate 1. }
  \item{mu}{ mean }
  \item{sigma}{ standard deviation }
  \item{rho}{ correlation coefficient }

}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  Log density of bivariate Gaussian distribution N(mu, mu, sigma, sigma, rho).
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

%\seealso{
  
%}
\examples{

z.1 <- rnorm(500, 3, 1)
rho <- 0.8

## The component with higher values is correlated with correlation coefficient=0.8 
z.2 <- rnorm(500, 3 + 0.8*(z.1-3), (1-rho^2))
mu <- 3
sigma <- 1
den.z <- d.binormal(z.1, z.2, mu, sigma, rho)

den.z
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
