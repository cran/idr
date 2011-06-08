\name{get.uri.2d}
\alias{get.uri.2d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Compute the correspondence profile (Psi(t, v)) and the derivative of
  the correspondence profile (Psi'(t, v)) 
}
\description{
An internal function to compute the correspondence profile (Psi(t, v)) and the derivative of the correspondence profile (Psi'(t, v)) for t and v between 0 and 1 
}
\usage{
get.uri.2d(x1, x2, tt, vv, spline.df = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x1}{
%%     ~~Describe \code{x1} here~~
}
  \item{x2}{
%%     ~~Describe \code{x2} here~~
}
  \item{tt}{
A vector of right-tail percentages for x1. It is a numeric vector between 0 and 1, e.g. seq(0.01, 0.99, by=0.01) 
}
  \item{vv}{
A vector of right-tail percentages for x2. It is a numeric vector between 0 and 1, e.g. seq(0.01, 0.99, by=0.01)  
}
 \item{spline.df}{
The degree of freedom to control the smoothness of the spline to be
superimposed on the points  
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
\references{
  %% ~put references to the literature/web site here ~
  Q. Li, J. B. Brown, H. Huang and P. J. Bickel. (2011)
   Measuring reproducibility of high-throughput experiments. Annals of Applied Statistics, Vol. 5, No. 3, 1752-1779.
}
\author{
  %%  ~~who you are~~
  Qunhua Li
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  comp.uri, get.correspondence
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

# salmon data from Kallenberg and Ledwina, 1999
data(salmon)

# get.correspondence() needs large ranks have high correlation and
# small ranks have low correlation. In this dataset, small values
# have high correlation and large values have low correlation.
# Ranking negative values makes the data follow the structure required
# by get.correspondence().
# There are 28 observations in this data set.

rank.x <- rank(-salmon$spawners)
rank.y <- rank(-salmon$recruits)
t <- seq(0.01, 0.99, by=1/28)
psi.dpsi <- get.uri.2d(rank.x, rank.y, t, t, spline.df=6.4)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
