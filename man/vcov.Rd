\name{vcov}
\alias{vcov}
\title{Variance-Covariance Matrix}
\description{
  Extract the variance-covariance matrix for (regression coefficients
  from) fitted model objects.
}
\usage{
vcov(object, \dots)
}
\arguments{
  \item{object}{an \R\ object for which a \code{vcov()} method exists,
    e.g. of class \code{lm}, \code{mlm}, \code{aov}, \code{aovlist},
    \code{glm}, \code{nls}, \code{l1ce} or \code{l1celist}.}
  \item{\dots}{further potential arguments passed to methods.}
}
\value{
  variance-covariance matrix of the corresponding model parameter vector.
}
\keyword{models}
% Converted by Sd2Rd version 1.21.
