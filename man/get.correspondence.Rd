\name{get.correspondence}
\alias{get.correspondence}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Compute correspondence profiles
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Compute the correspondence profiles (Psi and Psi') and the
  corresponding smoothed curve using spline
}
\usage{
get.correspondence(x1, x2, t, spline.df = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x1}{
    %%     ~~Describe \code{x1} here~~
    Data values or ranks of the data values on list 1, a vector of
    numeric values. Large values need to be significant signals. If small
    values represent significant signals, rank the signals reversely
    (e.g. by ranking negative values) and use the rank as x1.
}
  \item{x2}{
    %%     ~~Describe \code{x2} here~~
    Data values or ranks of the data values on list 2, a vector of
    numeric values. Large values need to be significant signals. If small
    values represent significant signals, rank the signals reversely
    (e.g. by ranking negative values) and use the rank as x1.
}
  \item{t}{
    %%     ~~Describe \code{t} here~~
    A numeric vector between 0 and 1 in ascending order. t is the
    right-tail percentage.
}
  \item{spline.df}{
    %%     ~~Describe \code{spline.df} here~~
    Degree of freedom for spline, to control the smoothness of the
    smoothed curve. 
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
  \item{psi}{the correspondence profile in terms of the scale of percentage, i.e. between (0, 1)}
  \item{dpsi}{the derivative of the correspondence profile in terms of the scale of percentage, i.e. between (0, 1)}
  \item{psi.n}{the correspondence profile in terms of the scale of the number of observations}
  \item{dpsi.n}{the derivative of the correspondence profile in terms of
    the scale of the number of observations}

  Each object above is a list consisting of the following items:
  t: upper percentage (for psi and dpsi) or number of top ranked
  observations (for psi.n and dpsi.n) 
  value: psi or dpsi
  smoothed.line: smoothing spline
  ntotal: the number of observations
  jump.point: the index of the vector of t such that psi(t[jump.point])
              jumps up due to ties at the low values. This only
	      happends when data consists of a large number of discrete
	      values, e.g. values imputed for observations
	      appearing on only one replicate. 
}
\references{
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

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{

# salmon data 
data(salmon)

# get.correspondence() needs the observations with high ranks have 
# high correlation and the observations with low ranks have low correlation. 
# In this dataset, small values have high correlation and large values 
# have low correlation.
# Ranking negative values makes the data follow the structure required
# by get.correspondence().
# There are 28 observations in this data set.

rank.x <- rank(-salmon$spawners)
rank.y <- rank(-salmon$recruits)
uv <- get.correspondence(rank.x, rank.y, seq(0.01, 0.99, by=1/28)) 


# plot correspondence curve on the scale of percentage
plot(uv$psi$t, uv$psi$value, xlab="t", ylab="psi", xlim=c(0, max(uv$psi$t)), 
ylim=c(0, max(uv$psi$value)), cex.lab=2)
lines(uv$psi$smoothed.line, lwd=4)
abline(coef=c(0,1), lty=3)

# plot the derivative of correspondence curve on the scale of percentage
plot(uv$dpsi$t, uv$dpsi$value, xlab="t", ylab="psi'", xlim=c(0, max(uv$dpsi$t)),
ylim=c(0, max(uv$dpsi$value)), cex.lab=2)
lines(uv$dpsi$smoothed.line, lwd=4)
abline(h=1, lty=3)

# plot correspondence curve on the scale of the number of observations
plot(uv$psi.n$t, uv$psi.n$value, xlab="t", ylab="psi", xlim=c(0, max(uv$psi.n$t)), 
ylim=c(0, max(uv$psi.n$value)), cex.lab=2)
lines(uv$psi.n$smoothed.line, lwd=4)
abline(coef=c(0,1), lty=3)

# plot the derivative of correspondence curve on the scale of the number
# of observations
plot(uv$dpsi.n$t, uv$dpsi.n$value, xlab="t", ylab="psi'", xlim=c(0, max(uv$dpsi.n$t)), 
ylim=c(0, max(uv$dpsi.n$value)), cex.lab=2)
lines(uv$dpsi.n$smoothed.line, lwd=4)
abline(h=1, lty=3)

# If the rank lists consist of a large number of ties at the bottom
# (e.g. the same low value is imputed to the list for the observations
# that appear on only one list), it may be desirable to plot only
# observations before hitting the ties. Then it can be plotted using the
# following
plot(uv$psi$t[1:uv$psi$jump.point], uv$psi$value[1:uv$psi$jump.point], xlab="t",
ylab="psi", xlim=c(0, max(uv$psi$t[1:uv$psi$jump.point])), 
ylim=c(0, max(uv$psi$value[1:uv$psi$jump.point])), cex.lab=2)
lines(uv$psi$smoothed.line, lwd=4)
abline(coef=c(0,1), lty=3)
 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Statistical models }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
