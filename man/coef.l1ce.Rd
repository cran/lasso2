% Copyright (C) 1998
% Berwin A. Turlach <bturlach@stats.adelaide.edu.au>
% Bill Venables <wvenable@stats.adelaide.edu.au>
% $Id: coef.l1ce.d,v 1.5 1999/12/01 05:08:40 bturlach Exp $
% --> ../COPYRIGHT for more details
\name{coef.l1ce}
\alias{coef.l1ce}
\title{Coefficients of an `l1ce' Object}
\usage{
coef.l1ce(object, all=TRUE, constrained=FALSE, \dots)
}
\description{
  This is a method for the function \code{coef()} for objects inheriting from
  class \code{l1ce}.
  See \code{\link{coef}} or \code{\link{coef.default}} for the general
  behavior of this function and for the interpretation of \code{object}.
}
\arguments{
  \item{object}{an object of class \code{l1ce}, see help on
    \code{\link{l1ce.object}}.}
  \item{all}{logical; if false, then only the non-zero coefficients are
    returned.}
  \item{constrained}{logical; if true, then only the coefficients that
    were constrained are returned.}
  \item{\dots}{possibly further arguments (none at the moment).}
}
\keyword{regression}
\keyword{models}
