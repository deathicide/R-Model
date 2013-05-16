\name{centerTMD}
\alias{centerTMD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
centerTMD(TMD, TMDendpts)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{TMD}{
%%     ~~Describe \code{TMD} here~~
}
  \item{TMDendpts}{
%%     ~~Describe \code{TMDendpts} here~~
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
function (TMD, TMDendpts) 
{
    if (missing(TMDendpts)) {
        TMDendpts <- fitLine(TMD)
    }
    return(matrix(c(TMD[, 1] - TMDendpts[1, 1], TMD[, 2] - TMDendpts[1, 
        2], TMD[, 3] - TMDendpts[1, 3]), ncol = 3))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
