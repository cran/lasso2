% Copyright (C) 1998
% Berwin A. Turlach <bturlach@stats.adelaide.edu.au>
% Bill Venables <wvenable@stats.adelaide.edu.au>
% $Id: labels.l1ce.d,v 1.2 1999/12/01 05:11:03 bturlach Exp $
% --> ../COPYRIGHT for more details
\name{labels.l1ce}
\alias{labels.l1ce}
\alias{labels.l1celist}
\title{`Labels' Method for `l1ce' and `l1celist' Objects}
\usage{
labels.l1ce(object, \dots)
labels.l1celist(object, \dots)
}
\arguments{
  \item{object}{fitted model of class \code{"l1ce"} or
    \code{"l1celist"}, respectively.}
  \item{\dots}{potentially further arguments passed to method.}
}
\description{
  This is a method for the \code{\link{labels}()} function for objects
  inheriting from class \code{l1ce} or \code{l1celist} (see
  \code{help(\link{l1ce.object})} and \code{help(\link{l1celist.object})}).
  See \code{\link{labels}} or \code{\link{labels.default}} for the
  general behavior of this function and for the interpretation of
  \code{object}.
}
\keyword{utilities}
