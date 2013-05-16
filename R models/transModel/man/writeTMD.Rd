\name{writeTMD}
\alias{writeTMD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
writeTMD(tmd, filename, tmdNum, domains, totalAtoms)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tmd}{
%%     ~~Describe \code{tmd} here~~
}
  \item{filename}{
%%     ~~Describe \code{filename} here~~
}
  \item{tmdNum}{
%%     ~~Describe \code{tmdNum} here~~
}
  \item{domains}{
%%     ~~Describe \code{domains} here~~
}
  \item{totalAtoms}{
%%     ~~Describe \code{totalAtoms} here~~
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
function (tmd, filename, tmdNum, domains, totalAtoms) 
{
    write.table(matrix(c(replicate(domains$end[tmdNum] - domains$start[tmdNum] + 
        1, "ATOM"), (1:(domains$end[tmdNum] - domains$start[tmdNum] + 
        1)), as.character(totalAtoms$name[c(domains$start[tmdNum]:domains$end[tmdNum])]), 
        as.character(totalAtoms$AA[c(domains$start[tmdNum]:domains$end[tmdNum])]), 
        totalAtoms$AANum[c(domains$start[tmdNum]:domains$end[tmdNum])], 
        totalAtoms$x[tmd[, 1]], totalAtoms$y[tmd[, 2]], totalAtoms$z[tmd[, 
            3]], replicate(domains$end[tmdNum] - domains$start[tmdNum] + 
            1, "1.00"), replicate(domains$end[tmdNum] - domains$start[tmdNum] + 
            1, "0.00"), as.character(totalAtoms$symbol[c(domains$start[tmdNum]:domains$end[tmdNum])])), 
        ncol = 11, byrow = FALSE), file = filename, row.names = FALSE, 
        quote = FALSE, sep = ",", col.names = FALSE)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
