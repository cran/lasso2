% Copyright (C) 1998
% Berwin A. Turlach <bturlach@stats.adelaide.edu.au>
% Bill Venables <wvenable@stats.adelaide.edu.au>
% $Id: residuals.l1ce.d,v 1.2 1999/12/01 05:11:45 bturlach Exp $
% --> ../COPYRIGHT for more details
\name{residuals.l1ce}
\alias{residuals.l1ce}
\alias{residuals.l1celist}
\title{Residuals of `l1ce' or `l1celist' Objects}
\usage{
residuals.l1ce(object, type, \dots)
residuals.l1celist(object, type, \dots)
}
\arguments{
  \item{object}{
    an object inheriting from class \code{l1ce} representing a fitted model.
  }
  \item{type}{
    type of residuals, with choices "deviance",  "pearson", "working" or
    "response"; the first is the default.
  }
 \item{\dots}{possibly further arguments (none at the moment).}
}
\description{
  This is a method for the function \code{\link{residuals}} for objects
  inheriting from class \code{l1ce} or \code{l1celist} (see
  \code{help(\link{l1ce.object})} and \code{help(\link{l1celist.object})}).
  See \code{\link{residuals}} or \code{\link{residuals.default}} for the
  general behavior of this function and for the interpretation of
  \code{object} and \code{type}.
}
\keyword{models}
% Converted by Sd2Rd version 1.21.
