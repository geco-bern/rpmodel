\name{calc_ftemp_inst_jmax}
\alias{calc_ftemp_inst_jmax}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{calc_ftemp_inst_jmax
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
calc_ftemp_inst_jmax(tcleaf, tcgrowth = tcleaf, tchome = NA, tcref = 25, kuma_par = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tcleaf}{
%%     ~~Describe \code{tcleaf} here~~
}
  \item{tcgrowth}{
%%     ~~Describe \code{tcgrowth} here~~
}
  \item{tchome}{
%%     ~~Describe \code{tchome} here~~
}
  \item{tcref}{
%%     ~~Describe \code{tcref} here~~
}
  \item{kuma_par}{
%%     ~~Describe \code{kuma_par} here~~
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
function (tcleaf, tcgrowth = tcleaf, tchome = NA, tcref = 25, 
    kuma_par = T) 
{
    Hd <- 2e+05
    tkref <- tcref + 273.15
    tkleaf <- tcleaf + 273.15
    Ha <- 49884
    Rgas <- 8.3145
    a_ent <- 659.7
    b_ent <- 0.75
    dent <- a_ent - b_ent * tcgrowth
    if (kuma_par) {
        Ha = 40710
        a_ent = 658.77
        b_ent = 0.84
        c_ent = 0.52
        dent = a_ent - (b_ent * tchome) - c_ent * (tcgrowth - 
            tchome)
    }
    fva <- calc_ftemp_arrh(tkleaf, Ha, tkref = tkref)
    fvb <- (1 + exp((tkref * dent - Hd)/(Rgas * tkref)))/(1 + 
        exp((tkleaf * dent - Hd)/(Rgas * tkleaf)))
    fv <- fva * fvb
    return(fv)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
