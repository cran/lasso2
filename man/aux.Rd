% Copyright (C) 1998
% Berwin A. Turlach <bturlach@stats.adelaide.edu.au>
% Bill Venables <wvenable@stats.adelaide.edu.au>
% $Id: aux.d,v 1.1 1998/09/18 11:12:48 bturlach Exp $
% --> ../COPYRIGHT for more details
\name{aux}
\alias{aux}
\title{Extract Auxiliary Information From an Object}
\description{
  Generic function for extracting auxiliary information from fitted
  model objects.
}
\usage{
aux(object, \dots)
}
\arguments{
  \item{object}{fitted model object (here typically of class
    \code{l1celist}, see \code{\link{aux.l1celist}}).}
  \item{\dots}{potentilly further arguments passed to methods.}
}
\value{
  a matrix with the bound(s) (relative [if used] and absolute) and
  the Lagrangian(s) for the fitted model(s).
}
\details{
  See documentation (technical reports).
}
\keyword{utilities}
