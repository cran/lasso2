% Copyright (C) 1998
% Berwin A. Turlach <bturlach@stats.adelaide.edu.au>
% Bill Venables <wvenable@stats.adelaide.edu.au>
% $Id: tr.d,v 1.1 1998/12/18 06:41:47 bturlach Exp $
% --> ../COPYRIGHT for more details
\name{tr}
\alias{tr}
\title{Trace of a Matrix}
\description{Calculates the trace of a matrix}
\usage{
tr(mat)
}
\arguments{
  \item{mat}{a square matrix.}
}
\value{
The trace of the matrix, i.e. the sum of its diagonal elements, is returned.
}
\examples{
tr(cbind(1,1:3,4:2)) # 5
}
\keyword{math}
% Converted by Sd2Rd version 1.21.
