% Copyright (C) 1998
% Berwin A. Turlach <bturlach@stats.adelaide.edu.au>
% Bill Venables <wvenable@stats.adelaide.edu.au>
% $Id: gcv.d,v 1.2 1998/09/18 07:03:35 bturlach Exp $
% --> ../COPYRIGHT for more details
\name{gcv}
\alias{gcv}
\title{Generalised Cross-Validation Score}
\description{
  Extracts the generalised cross-validation score(s) from fitted
  model objects.
}
\usage{
gcv(object, \dots)
}
\arguments{
  \item{object}{fitted model object; see gcv methods for details.}
  \item{\dots}{arguments passed to methods.}
}
\value{
  A vector (or matrix) with the bound(s) (relative [if used] and absolute),
  the Lagrangian(s) and the generalised cross-validation score(s) for
  the fitted model(s).
}
\details{
  See documentation.
}
\keyword{models}
% Converted by Sd2Rd version 1.21.
