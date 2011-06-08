\name{comp.uri}
\alias{comp.uri}
\title{
  Compute the value of Psi and Psi' for a given t
}
\description{
An internal function to compute Psi for a given (t, v)
}
\usage{
comp.uri(tv, x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tv}{A vector of two numeric values, t and v. Both t and
    v are in [0,1]
}
  \item{x}{A numeric vector of x, sorted by the order of y
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A numeric value of Psi(t, v)
}
\references{
    Q. Li, J. B. Brown, H. Huang and P. J. Bickel. (2011)
Measuring reproducibility of high-throughput experiments. Annals of Applied Statistics. In press.
}
\author{
  Qunhua Li
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
get.uri.2d, get.correspondence
}
\examples{
x <- seq(1, 10, by=1)
tv <- c(0.1, 0.5)

# opposite order
y1 <- seq(10, 1, by=-1)
o1 <- order(y1)
x.o1 <- x[o1]

comp.uri(tv, x.o1)

# same order
y2 <- seq(1, 10, by=1)
o2 <- order(y2)
x.o2 <- x[o2]

comp.uri(tv, x.o2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
