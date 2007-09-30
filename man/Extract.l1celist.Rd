% Copyright (C) 1998
% Berwin A. Turlach <bturlach@stats.adelaide.edu.au>
% Bill Venables <wvenable@stats.adelaide.edu.au>
% $Id: [[.l1celist.d,v 1.1 1998/09/18 11:23:09 bturlach Exp $
% --> ../COPYRIGHT for more details
\name{Extract.l1celist}
\alias{[.l1celist}
\alias{[[.l1celist}
\title{Extract Parts of a `l1celist' Object}
\description{
  Allows the user to extract values from a \code{l1celist} object by using
  subscripts.
}
\usage{
\method{[}{l1celist}(x, \dots, drop = TRUE)
\method{[[}{l1celist}(x, \dots, drop = TRUE)
}
\arguments{
  \item{x}{an object inheriting from class \code{"l1celist"}.}
  \item{\dots}{a specification of indices -- see \code{\link{Extract}}.}
  \item{drop}{logical defaulting to \code{TRUE}.  If only one model is
    subscribed, then it is returned as an object of class \code{"l1ce"}.
    If \code{drop=F}, then an object of class \code{"l1celist"} is
    always returned.}
}
\value{
  an object of class \code{"l1celist"} or class \code{"l1ce"}
  extracted from the original list.
}
\keyword{manip}
% Converted by Sd2Rd version 1.21.
