% Copyright (C) 1998
% Berwin A. Turlach <bturlach@stats.adelaide.edu.au>
% Bill Venables <wvenable@stats.adelaide.edu.au>
% $Id: print.l1ce.d,v 1.2 1999/12/01 05:11:22 bturlach Exp $
\name{print.l1ce}
\alias{print.l1ce}
\alias{print.l1celist}
\alias{print.gl1ce}
\title{Print Methods for `l1ce', `l1celist' and `gl1ce' Objects}
\usage{
\method{print}{l1ce}(x, \dots)
\method{print}{l1celist}(x, \dots)
\method{print}{gl1ce}(x, \dots)
}
\arguments{
  \item{x}{fitted model of class \code{"l1ce"} or
    \code{"l1celist"}, respectively.}
  \item{\dots}{potentially further arguments passed to method.}
}
\description{
  These are methods of the generic \code{\link{print}()} function for objects
  inheriting from class \code{l1ce}, \code{l1celist} or \code{gl1ce} (see
  \code{help(\link{l1ce.object})}, \code{help(\link{l1celist.object})}
  or \code{help(\link{gl1ce.object})}).  See \code{\link{print}} or
  \code{\link{print.default}} for the general behavior of this
  function.
}
\keyword{print}
