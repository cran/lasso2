\name{residuals.gl1ce}
\alias{residuals.gl1ce}
\title{Compute Residuals for `gl1ce' Objects}
\description{
  Computes one of the four types of residuals available for \code{gl1ce}
  objects.  This is a method for the function residuals() for  objects
  inheriting from class \code{gl1ce}.  As several types of residuals are
  available for \code{gl1ce} objects, there is an additional optional
  argument type.
}
\usage{
\method{residuals}{gl1ce}(object,
      type=c("deviance", "pearson", "working", "response"), \dots)
}
\arguments{
  \item{object}{
    an object inheriting from class \code{gl1ce} representing a fitted model.
  }
  \item{type}{
    type of residuals, with choices "deviance",  "pearson", "working" or
    "response"; the first is the default.
  }
 \item{\dots}{possibly further arguments (none at the moment).}
}
\value{
  A numeric vector of residuals. See Statistical Models in S for detailed
  definitions of each type of residual. The sum of squared deviance residuals
  add up to the deviance.  The pearson residuals are standardized residuals on
  the scale of the response.  The  working residuals reside on the
  object, and are the residuals from the final fit. The response residuals are
  simply  y-fitted(object).  The summary() method for \code{gl1ce}
  objects produces deviance residuals.  The residuals component of a
  \code{gl1ce} object contains the working residuals.
}
\references{
  Chambers,  J.M.,  and  Hastie,  T.J.  (1991).  Statistical Models in S,
  pp. 204--206.
}
\keyword{models}
% Converted by Sd2Rd version 1.21.
