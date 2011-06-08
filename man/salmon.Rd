\name{salmon}
\alias{salmon}
\docType{data}
\title{
  Salmon data
}
\description{
  This data is from Simonoff (1990, p 161). It concerns the
  size of the annual spawning stock and its production of new
  catchable-sized fish for 1940 through 1967 for the Skeena river
  sockeye salmon stock (in thousands of fish). It has three columns,
  year, spawners and recruits. It was speculated to consist of two
  different populations.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(salmon)}
\format{
  A data frame with 28 observations on the following 3 variables.
  \describe{
    \item{\code{year}}{a numeric vector of the year}
    \item{\code{spawners}}{a numeric vector of the annual spawning stock}
    \item{\code{recruits}}{a numeric vector of the production of new
      catchable-sized fish}
  }
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
  Data is salmon.dat in Simonoff (1990). It can be downloaded from the
  book's website. 
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
  Simonoff, J. S. (1990), Smoothing Methods in Statistics, New York: Spinger-Verlag. 
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(salmon)
plot(rank(salmon$spawners), rank(salmon$recruits)) 
}
\keyword{datasets}
