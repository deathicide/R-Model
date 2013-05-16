\name{plotTMD}
\alias{plotTMD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotTMD(TMD, endpt, totalTMD, adder, datacol, fitcol)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{TMD}{
%%     ~~Describe \code{TMD} here~~
}
  \item{endpt}{
%%     ~~Describe \code{endpt} here~~
}
  \item{totalTMD}{
%%     ~~Describe \code{totalTMD} here~~
}
  \item{adder}{
%%     ~~Describe \code{adder} here~~
}
  \item{datacol}{
%%     ~~Describe \code{datacol} here~~
}
  \item{fitcol}{
%%     ~~Describe \code{fitcol} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (TMD, endpt, totalTMD, adder, datacol, fitcol) 
{
    if (missing(adder)) {
        adder = T
    }
    if (missing(datacol)) {
        datacol = "blue"
    }
    if (missing(fitcol)) {
        fitcol = "red"
    }
    if (as.logical(adder)) {
        plot3d(TMD[, 1], TMD[, 2], TMD[, 3], col = datacol, xlim = range(totalTMD[, 
            1]), ylim = range(totalTMD[, 2]), zlim = range(totalTMD[, 
            3]), xlab = "\305", ylab = "\305", zlab = "\305", 
            add = F, type = "s", size = 1)
    }
    else {
        plot3d(TMD[, 1], TMD[, 2], TMD[, 3], col = datacol, xlim = range(totalTMD[, 
            1]), ylim = range(totalTMD[, 2]), zlim = range(totalTMD[, 
            3]), xlab = "\305", ylab = "\305", zlab = "\305", 
            add = T, type = "s", size = 1)
    }
    plot3d(endpt[, 1], endpt[, 2], endpt[, 3], type = "l", col = fitcol, 
        lwd = 2, add = T)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
