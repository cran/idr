\name{e.step.2normal}
\alias{e.step.2normal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  E-step for parameterized bivariate 2-component Gaussian mixture models
}
\description{
Expectation step in the EM algorithm for parameterized bivariate 2-component
Gaussian mixture models with (1-p)N(0, 0, 1, 1, 0) + pN(mu, mu, sigma,
sigma, rho).
}
\usage{
e.step.2normal(z.1, z.2, mu, sigma, rho, p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{z.1}{ a numerical data vector on coordinate 1.}
  \item{z.2}{ a numerical data vector on coordinate 2.}
  \item{mu}{ mean for the reproducible component.}
  \item{sigma}{ standard deviation of the reproducible component.}
  \item{rho}{ correlation coefficient of the reproducible component.}
  \item{p}{ mixing proportion of the reproducible component.}

}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  \item{e.z}{ a numeric vector, where each entry represents the estimated expected conditional probability
    that an observation is in the reproducible component. }
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
  \code{m.step.2normal}, \code{loglik.2binormal}, \code{est.IDR}
}
\examples{

z.1 <- c(rnorm(500, 0, 1), rnorm(500, 3, 1))
rho <- 0.8

## The component with higher values is correlated with correlation coefficient=0.8 
z.2 <- c(rnorm(500, 0, 1), rnorm(500, 3 + 0.8*(z.1[501:1000]-3), (1-rho^2)))

## Starting values
mu0 <- 3
sigma0 <- 1
rho0 <- 0.85
p0 <- 0.55

e.z <- e.step.2normal(z.1, z.2, mu0, sigma0, rho0, p0)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
