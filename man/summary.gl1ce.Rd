\name{summary.gl1ce}
\alias{summary.gl1ce}
\alias{print.summary.gl1ce}
\title{Summary Method for Generalized L1 Constrained Regression Models}
\description{
  Returns (and prints) a summary list for fitted generalized L1
  constrained regression models.
}
\usage{
\method{summary}{gl1ce}(object, dispersion = NULL, correlation = FALSE, \dots)
\method{print}{summary.gl1ce}(x, digits = max(3, getOption("digits") - 3), \dots)
}
\arguments{
  \item{object}{
    fitted \code{gl1ce} model object.  This is assumed to be the
    result of some fit that produces an object inheriting from
    the class \code{gl1ce}, in the sense that the components returned
    by the gl1ce() function will be available.
  }
  \item{dispersion}{prescribed dispersion, see \code{\link{summary.glm}}.}
  \item{correlation}{logical indicating if the correlation matrix should
    be returned as well.}
  \item{x}{object of class \code{summary.gl1ce}.}
  \item{digits, \dots}{further arguments for the \code{print} method.}
}
\value{
  an object of class \code{summary.gl1ce} (for which there's a
  \code{print} method).
  It is basically a list with the following components:

  \item{call}{as contained on object}
  \item{residuals}{the deviance residuals, as produced by
    \code{residuals(object, type = "deviance")}.}
  \item{coefficients}{the coefficients of the model.}
  \item{family}{
    the family of models to which object belongs, along with the
    variance and link functions for that model.}
  \item{bound}{the bound used in fitting this model}
  \item{Lagrangian}{the Lagrangian of the model}
}
\details{
  This function is a method for the generic function summary for  class
  \code{gl1ce}.
  It can be invoked by calling summary for an object of the appropriate class,
  or directly by calling  summary.gl1ce  regardless of the class of the object.
}
\keyword{models}
% Converted by Sd2Rd version 1.21.
