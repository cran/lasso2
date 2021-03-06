% Copyright (C) 1998
% Berwin A. Turlach <bturlach@stats.adelaide.edu.au>
% Bill Venables <wvenable@stats.adelaide.edu.au>
% $Id: coef.l1ce.d,v 1.5 1999/12/01 05:08:40 bturlach Exp $
% --> ../COPYRIGHT for more details
\name{coef.l1ce}
\alias{coef.l1ce}
\title{Coefficients of an `l1ce' Object}
\usage{
\method{coef}{l1ce}(object, all=TRUE, constrained=FALSE, \dots)
}
\description{
%% For R >= 1.9.0, replace  \link{   by   \link[stats]{:
  This is a method for \code{\link{coef}()} for objects inheriting from
  class \code{l1ce}.  See \code{\link{coef}} for the general behavior of
  this function and for the interpretation of \code{object}.
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
