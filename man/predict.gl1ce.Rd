\name{predict.gl1ce}
\alias{predict.gl1ce}
\title{
Use `predict()' on a `gl1ce' object
}
\description{
Extracts the fitted values from a \code{gl1ce()} object and returns a matrix of
predictions.
}
\usage{
predict.gl1ce(object, newdata, type=c("link", "response"),
              se.fit = FALSE, \dots)
}
\arguments{
  \item{object}{
    a fitted \code{gl1ce} object.
  }
  \item{newdata}{
    a data frame containing the values at which predictions are required.
    This argument can be missing, in which case predictions are made at
    the same values used to compute the object.
    Only those predictors referred to in the right side of the formula in
    object need be present by name in \code{newdata}.
  }
  \item{type}{
    type of predictions, with choices "link" (the default), or
    "response".  The default produces predictions on the scale of the
    additive predictors, and with \code{newdata} missing,
    \code{predict()} is simply an extractor function for this component
    of a \code{gl1ce} object.  If "response" is selected, the
    predictions are on the scale of the response, and are monotone
    transformations of the additive predictors, using the inverse link
    function.}
  \item{se.fit}{logical indicating if standard errors should be returned
    as well. \bold{Not yet available}.}
  \item{\dots}{further potential arguments passed to methods.}
}
\value{
  a vector of predictions.
}
\details{
This function is a method for the generic function predict
for  class  \code{gl1ce}.   It can be invoked by calling predict
for an object of the appropriate  class,  or  directly  by
calling  predict.gl1ce  regardless  of the class of the object.
}
\section{Warning}{
\code{predict} can produce incorrect predictions when the \code{newdata}
argument is used if the formula in \code{object} involves
data-dependent transformations, such as \code{poly(Age, 3)}
or \code{sqrt(Age - min(Age))}.
}
\keyword{models}
% Converted by Sd2Rd version 1.21.
