% Copyright (C) 1998
% Berwin A. Turlach <bturlach@stats.adelaide.edu.au>
% Bill Venables <wvenable@stats.adelaide.edu.au>
% --> ../COPYRIGHT for more details
\name{merge.formula}
\alias{merge.formula}
\title{Merge Formula With Right Hand Side of Second Formula}
\usage{
merge.formula(x, y, \dots)
}
\arguments{
  \item{x,y}{formulas.}
  \item{\dots}{potentially further arguments passed to methods.}
}
\description{
  This is method for formulas of the \code{\link[base]{merge}} generic
  function.  Here it is support for the function \code{\link{l1ce}} and
  not intended to be called directly by users.
}
\examples{
merge.formula(y ~ x1, ~ x2) ## ->  y ~ x1 + x2

f2 <- merge.formula(y ~ x1*x2, z ~ (x2+x4)^3)
f. <- merge.formula(y ~ x1*x2,   ~ (x2+x4)^3) # no LHS for 2nd term
f2
stopifnot(f2 == f.)
}
\keyword{utilities}
