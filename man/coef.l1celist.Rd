% Copyright (C) 1998
% Berwin A. Turlach <bturlach@stats.adelaide.edu.au>
% Bill Venables <wvenable@stats.adelaide.edu.au>
% $Id: coef.l1celist.d,v 1.4 1999/12/01 05:08:55 bturlach Exp $
% --> ../COPYRIGHT for more details
\name{coef.l1celist}
\alias{coef.l1celist}
\title{Coefficients of an `l1celist' Object}
\usage{
coef.l1celist(object, all=TRUE, constrained=FALSE, \dots)
}
\description{
  This is a method for \code{\link{coef}()} for objects inheriting from
  class \code{l1ce}.
  See \code{\link{coef}} or \code{\link{coef.default}} for the general
  behavior of this function and for the interpretation of \code{object}.
}
\arguments{
  \item{object}{an object of class \code{l1celist}, see help on
    \code{\link{l1celist.object}}.}
  \item{all}{logical; if false, then the coefficients that are zero in
    all fitted models of the list are not returned.}
  \item{constrained}{logical; if true, then only the coefficients that
    were constrained are returned.}
  \item{\dots}{possibly further arguments (none at the moment).}
}
\keyword{regression}
\keyword{models}
