% Copyright (C) 1998
% Berwin A. Turlach <bturlach@stats.adelaide.edu.au>
% Bill Venables <wvenable@stats.adelaide.edu.au>
% $Id: predict.l1ce.d,v 1.1 1998/09/18 07:27:59 bturlach Exp $
%
% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Library General Public
% License as published by the Free Software Foundation; either
% version 2 of the License, or (at your option) any later version.
%
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Library General Public License for more details.
%
% You should have received a copy of the GNU Library General Public
% License along with this library; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston,
% MA 02111-1307 USA
\name{predict.l1ce}
\alias{predict.l1ce}
\title{Predict Method for `l1ce' Objects}
\description{
  This is a method for the generic function \code{\link{predict}}
  for class \code{"l1ce"}, typically produced from
  \code{\link{l1ce}()}.  When \code{newdata} is missing, the
  fitted values are extracted, otherwise returns \emph{new} predictions.
}
\usage{
\method{predict}{l1ce}(object, newdata,
           type = c("response"), se.fit = FALSE, \dots)
}
\arguments{
  \item{object}{a fitted \code{l1ce} object.  }
  \item{newdata}{
    a data frame containing the values at which predictions are required.
    This argument can be missing, in which case predictions are made at
    the same values used to compute the object.
    Only those predictors referred to in the right side of the formula in
    object need be present by name in \code{newdata}.}
  \item{type}{currrently only \code{"response"}.}
  \item{se.fit}{logical indicating if standard errors should be returned
    as well. \bold{Not yet available}.}
  \item{\dots}{further potential arguments passed to methods.}
}
\value{
  a vector of predictions.
}
\section{Warning}{
  \code{predict} can produce incorrect predictions when the \code{newdata}
  argument is used if the formula in \code{object} involves
  transformations, such as \code{poly(Age, 3)} or \code{sqrt(Age - min(Age))}.

  MM: Not sure this is true for \R \ldots.
}
\examples{
data(Iowa)
l1c.I <- l1ce(Yield ~ ., Iowa, bound = 10, absolute.t=TRUE)
p10 <- predict(l1c.I, newdata =   Iowa[10:19,])
stopifnot(all.equal(p10, fitted(l1c.I)[10:19]))
}
\keyword{models}





