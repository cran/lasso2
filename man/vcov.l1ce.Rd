% Copyright (C) 1998
% Berwin A. Turlach <bturlach@stats.adelaide.edu.au>
% Bill Venables <wvenable@stats.adelaide.edu.au>
% $Id: vcov.l1ce.d,v 1.2 1999/12/01 05:11:56 bturlach Exp $
% --> ../COPYRIGHT for more details
\name{vcov.l1ce}
\alias{vcov.l1ce}
\alias{vcov.l1celist}
\title{Variance-Covariance Matrix of `l1ce' or `l1celist' Objects}
\usage{
\method{vcov}{l1ce}(object, type = c("OPT", "Tibshirani"),
    gen.inverse.diag = 0, \dots)
\method{vcov}{l1celist}(object, type = c("OPT", "Tibshirani"),
    gen.inverse.diag = 0, \dots)
}
\arguments{
  \item{object}{an object of class \code{l1ce} or \code{l1celist}.}
  \item{type}{character indicating whether to use the covariance formula
    of Osborne, Presnell and Turlach or the formula of Tibshirani.}
  \item{gen.inverse.diag}{if Tibshirani's formula for the covariance
    matrix is used, this value is used for the diagonal elements of the
    generalised inverse that appears in the formula that corresponds to
    parameters estimated to be zero.  The default is 0, i.e. use the
    Moore-Penrose inverse.  Tibshirani's code uses
    \code{gen.inverse.diag = 1e11}.}
  \item{\dots}{further potential arguments passed to methods.}
}
\description{
  This is a method for the function \code{vcov()} for objects inheriting
  from class \code{l1ce} or \code{l1celist} (see
  \code{help(\link{l1ce.object})} and
  \code{help(\link{l1celist.object})}).  See \code{\link{vcov}} for the
  general behavior of this function.
}
\keyword{models}
% Converted by Sd2Rd version 1.21.
