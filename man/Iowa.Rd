% Copyright (C) 1998
% Berwin A. Turlach <bturlach@stats.adelaide.edu.au>
% Bill Venables <wvenable@stats.adelaide.edu.au>
% $Id: Iowa.d,v 1.3 1999/12/01 05:07:43 bturlach Exp $
% --> ../COPYRIGHT for more details
\name{Iowa}
\alias{Iowa}
\title{The Iowa Wheat Yield Data}
\description{
  The data gives the pre-season and three growing months' precipitation,
  the mean temperatures for the three growing months and harvest month,
  the year, and the yield of wheat for the USA state of Iowa, for the
  years 1930--1962.
}
\usage{
data(Iowa)
}
\format{
  The data frame has the following components:
  \describe{
    \item{\code{Year}}{Year of measurement (surrogate for variety improvements)}
    \item{\code{Rain0}}{Pre-season rainfall (in.)}
    \item{\code{Temp1}}{Mean temperature for the first growing month (deg. F)}
    \item{\code{Rain1}}{Rainfall for the first growing month (in.)}
    \item{\code{Temp2}}{Mean temperature for the second growing month (deg. F)}
    \item{\code{Rain2}}{Rainfall for the second growing month (in.)}
    \item{\code{Temp3}}{Mean temperature for the third growing month (deg. F)}
    \item{\code{Rain3}}{Rainfall for the third growing month (in.)}
    \item{\code{Temp4}}{Mean temperature for the harvest month (deg. F)}
    \item{\code{Yield}}{Yield of wheat in Iowa for the given year (bush./acre)}
  }
}
\section{CATEGORY}{
Multiple regression; diagnostics.
}
\source{
  CAED Report, 1964.  Quoted in Draper and Smith, Applied Regression
  Analysis.
}
\examples{
data(Iowa)
pairs(Iowa)
}
\keyword{datasets}
% Converted by Sd2Rd version 1.21.
