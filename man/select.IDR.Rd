\name{select.IDR}
\alias{select.IDR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Select observations according to IDR
}
\description{
  Select observations that exceeding a given IDR level
}
\usage{
select.IDR(x, IDR.x, IDR.level)
}
\arguments{
  \item{x}{a n by m numeric matrix, where m= num of replicates,
    n=num of observations. Numerical values representing the
    significance of the observations, where larger values represent
    higher significance, for example, -log(p-value).   Currently, m=2. }
  \item{IDR.x}{Irreproducibile discovery rate for each entry of x. It is
  computed from est.IDR().}
  \item{IDR.level}{IDR cutoff, a numerical value between [0,1].}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  \item{x}{Observations that are selected.}
  \item{n}{Number of observations that are selected.}
  \item{IDR.level}{IDR cutoff, a numerical value between [0,1].} 
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
  \code{\link{est.IDR}}
}
\examples{

data("simu.idr")
x <- cbind(-simu.idr$x, -simu.idr$y)

mu <- 2.6
sigma <- 1.3
rho <- 0.8
p <- 0.7

idr.out <- est.IDR(x, mu, sigma, rho, p, eps=0.001, max.ite=20)
# select observations exceeding IDR threshold=0.01 
IDR.level <- 0.01
x.selected <- select.IDR(x, idr.out$IDR, IDR.level)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Statistical Models}
% __ONLY ONE__ keyword per line
