% Copyright (C) 1998
% Berwin A. Turlach <bturlach@stats.adelaide.edu.au>
% Bill Venables <wvenable@stats.adelaide.edu.au>
% $Id: plot.l1celist.d,v 1.2 1999/11/30 06:14:16 bturlach Exp $
% --> ../COPYRIGHT for more details
\name{plot.l1celist}
\alias{plot.l1celist}
\title{Plot Method for `l1celist' Objects}
\description{
   Plots a \code{l1celist} object on the current graphics device.
}
\usage{
plot.l1celist(x, plot=TRUE, all=TRUE, constrained=FALSE, \dots)
}
\arguments{
  \item{x}{fitted model object of class \code{l1celist}.}
  \item{plot}{logical; if \code{TRUE} a \code{matplot()} of all the
    coefficients in the list against the relative bound (absolute bound if
    relative is not present) is plotted.  Otherwise no plot is done.
  }
  \item{all}{logical; if \code{FALSE}, then only the non-zero
    coefficients are returned.}
  \item{constrained}{
    if \code{TRUE} then only the coefficients that were constrained
    are returned.}
 \item{\dots}{possibly further arguments (none at the moment).}
}
\value{
  A matrix with the bound(s) (relative [if used] and absolute),
  the Lagrangian(s) and coefficints of the fitted model(s).
}
\examples{
data(Prostate)
l1c.P <- l1ce(lpsa ~ ., Prostate, bound=(1:20)/20)
length(l1c.P)# 20 l1ce models
plot(l1c.P)
}
\keyword{hplot}
% Converted by Sd2Rd version 1.21.
