\name{get.pseudo.mix}
\alias{get.pseudo.mix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute the pseudo values of a mixture model from the empirical CDF
}
\description{
Compute the pseudo values of a mixture model from the empirical CDF  
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get.pseudo.mix(x, mu, sigma, rho, p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A vector of values of empirical CDF
}
  \item{mu}{
Mean of the reproducible component in the mixture model on the latent space
}
  \item{sigma}{
Standard deviation of the reproducible component in the mixture model on the latent space
}
  \item{rho}{
Correlation coefficient of the reproducible component in the mixture model on the latent space
}
  \item{p}{
Mixing proportion of the reproducible component in the mixture model on the latent space
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  The values of a mixture model corresponding to the empirical CDF
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
 Q. Li, J. B. Brown, H. Huang and P. J. Bickel. (2011)\\
 Measuring reproducibility of high-throughput experiments. Annals of Applied Statistics. In press.  
%% ~put references to the literature/web site here ~
}
\author{
  Qunhua Li
%%  ~~who you are~~
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{

x <- seq(0.1, 0.9, by=0.1)

mu <- 2.6
sigma <- 1.3
rho <- 0.8
p <- 0.7

pseudo.x <- get.pseudo.mix(x, mu, sigma, rho, p)

pseudo.x
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
