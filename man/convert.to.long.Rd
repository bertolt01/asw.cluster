\name{convert.to.long}
\alias{convert.to.long}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{convert.to.long
%%  ~~function to do ... ~~
}
\description{convert.to.long
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
convert.to.long(res,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res}{
%%     ~~Describe \code{res} here~~
}
  \item{...}{Other parameters.
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
function (res) 
{
    mts <- res$mbr_to_subgroups
    mbr <- res$fl.mbr
    res.out <- c()
    teams <- as.numeric(as.factor(unlist(res$team)))
    for (team in unique(teams)) {
        for (mb in 1:length(mbr[[team]])) {
            for (sg in 1:length(rbind(mts[[team]])[mb, ])) {
                res.out <- rbind(res.out, c(res$team[[team]], 
                  res$teamsize[[team]], res$fl.value[[team]][mb], 
                  res$fl.mbr[[team]][mb], rbind(mts[[team]])[mb, 
                    sg], res$number_of_subgroups[[team]][mb]))
            }
        }
    }
    res.out <- as.data.frame(res.out)
    names(res.out) <- names(res)[1:length(res.out)]
    return(res.out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
