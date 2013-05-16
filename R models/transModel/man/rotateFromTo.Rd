\name{rotateFromTo}
\alias{rotateFromTo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rotateFromTo(from, to, tolerance)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{from}{
%%     ~~Describe \code{from} here~~
}
  \item{to}{
%%     ~~Describe \code{to} here~~
}
  \item{tolerance}{
%%     ~~Describe \code{tolerance} here~~
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
function (from, to, tolerance) 
{
    if (missing(tolerance)) {
        tolerance <- 0.001
    }
    fromendpts <- fitLine(from)
    toendpts <- fitLine(to)
    fromCenter <- centerTMD(from, fromendpts)
    toCenter <- centerTMD(to, toendpts)
    fromCenterendpts <- fitLine(fromCenter)
    toCenterendpts <- fitLine(toCenter)
    while (angleBetweenVec(fromCenterendpts, toCenterendpts) > 
        tolerance) {
        Q <- 0 + (toCenterendpts[2, 1]/2 + fromCenterendpts[2, 
            1]/2) * Hi + (toCenterendpts[2, 2]/2 + fromCenterendpts[2, 
            2]/2) * Hj + (toCenterendpts[2, 3]/2 + fromCenterendpts[2, 
            3]/2) * Hk
        Q <- Q/sqrt(Norm(Q))
        fromCenter <- rotate(fromCenter, Q)
        fromCenterendpts <- fitLine(fromCenter)
        if (abs(fromCenterendpts[2, 1]) < abs(fromCenterendpts[1, 
            1]) && abs(fromCenterendpts[2, 2]) < abs(fromCenterendpts[1, 
            2]) && abs(fromCenterendpts[2, 3]) < abs(fromCenterendpts[1, 
            3])) {
            for (i in 1:3) {
                temp <- fromCenterendpts[2, i]
                fromCenterendpts[2, i] <- fromCenterendpts[1, 
                  i]
                fromCenterendpts[1, i] <- temp
            }
        }
    }
    return(matrix(c(fromCenter[, 1] + fromendpts[1, 1], fromCenter[, 
        2] + fromendpts[1, 2], fromCenter[, 3] + fromendpts[1, 
        3]), ncol = 3))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
