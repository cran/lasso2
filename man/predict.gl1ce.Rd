\name{predict.gl1ce}
\alias{predict.gl1ce}
\title{Prediction Method for a `gl1ce' Object}
\description{
  This is a method for the generic function \code{\link{predict}}
  for class \code{"gl1ce"}, typically produced from
  \code{\link{gl1ce}()}.  When \code{newdata} is missing, the
  fitted values are extracted, otherwise returns \emph{new} predictions.
}
\usage{
\method{predict}{gl1ce}(object, newdata, type=c("link", "response"),
              se.fit = FALSE, \dots)
}
\arguments{
  \item{object}{a fitted \code{\link{gl1ce}} object.}
  \item{newdata}{
    a data frame containing the values at which predictions are required.
    This argument can be missing, in which case predictions are made at
    the same values used to compute the object.
    Only those predictors referred to in the right side of the formula in
    object need be present by name in \code{newdata}.}
  \item{type}{type of predictions, with choices "link" (the default), or
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
\section{Warning}{
  \code{predict} can produce incorrect predictions when the \code{newdata}
  argument is used if the formula in \code{object} involves
  data-dependent transformations, such as \code{poly(Age, 3)}
%% MM: this was for S-plus -- which has such problems with glm() -- maybe not R
  or \code{sqrt(Age - min(Age))}.
}
\examples{
  ## start with
  example(gl1ce)
  predict(gl1c.E, new = esoph[1:7,])# type 'link'
  predict(gl1c.E, new = esoph[1:7,], type = "response")

  ## identities / consistency checks :
  stopifnot(predict(gl1c.E, type = "response") == fitted(gl1c.E),
            all.equal(predict(gl1c.E)[1:7],
                      as.vector(predict(gl1c.E, new = esoph[1:7,]))),
            all.equal(fitted(gl1c.E)[1:7],
                      as.vector(predict(gl1c.E, new = esoph[1:7,], type = "response")))
           )
}
\keyword{models}

