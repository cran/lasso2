\name{l1ce}
\alias{l1ce}
\title{Regression Fitting With L1-constraint on the Parameters}
\description{
  Returns an object of class \code{"l1ce"} or \code{"licelist"} that represents
  fit(s) of linear models while imposing L1 constraint(s) on the parameters.
}
\usage{
l1ce(formula, data = sys.parent(), weights, subset, na.action,
     sweep.out = ~ 1, x = FALSE, y = FALSE,
     contrasts = NULL, standardize = TRUE,
     trace = FALSE, guess.constrained.coefficients = double(p),
     bound = 0.5, absolute.t = FALSE)
}
\arguments{
  \item{formula}{a formula object, with the response on the left of a
    \code{ ~ } operator, and the terms, separated by \code{+} operators,
    on the right.}
  \item{data}{
    a \code{data.frame} in which to interpret the variables named in the
    formula, the \code{weights}, the \code{subset} and the
    \code{sweep.out} argument.  If this is missing, then the variables
    in the formula should be globally available.}
  \item{weights}{
    vector of observation weights.  The length of \code{weights} must be
    the same as the number of observations.  The weights must be
    nonnegative and it is strongly recommended that they be strictly
    positive, since zero weights are ambiguous, compared to use of the
    \code{subset} argument.
  }
  \item{subset}{
    expression saying which subset of the rows of the data should be
    used in the fit.  This can be a logical vector (which is replicated
    to have length equal to the number of observations), or a numeric
    vector indicating which observation numbers are to be included, or a
    character vector of the row names to be included.  All observations
    are included by default.
  }
  \item{na.action}{
    a function to filter missing data.  This is applied to the
    \code{\link{model.frame}} after any \code{subset} argument has been used.
    The default (with \code{na.fail}) is to create an error if any
    missing values are found.  A possible alternative is \code{na.omit},
    which deletes observations that contain one or more missing values.
  }
  \item{sweep.out}{
    a formula object, variables whose parameters are not put under the
    constraint are swept out first.  The variables should appear on the
    right of a \code{ ~ } operator and be separated by \code{+} operators.
    Default is \code{ ~1 }, i.e. the constant term is not under the
    constraint.  If this parameter is \code{NULL}, then all parameters
    are put under the constraint.
  }
  \item{x}{logical indicating if the model matrix should be returned in
    component \code{x}.}
  \item{y}{logical indicating if the response should be returned in
    component \code{y}.}
  \item{contrasts}{
    a list giving contrasts for some or all of the factors
    appearing in the model formula.
    The elements of the list should have the same name as the variable
    and should be either a contrast matrix (specifically, any full-rank
    matrix with as many rows as there are levels in the factor),
    or else a function to compute such a matrix given the number of levels.
  }
  \item{standardize}{
    logical flag: if \code{TRUE}, then the columns of the model matrix that
    correspond to parameters that are constrained are standardized to have
    emprical variance one.  The standardization is done after taking
    possible weights into account and after sweeping out variables whose
    parameters are not constrained.
  }
  \item{trace}{
    logical flag: if \code{TRUE}, then the status during each iteration of
    the fitting is reported.
  }
  \item{guess.constrained.coefficients}{
    initial guess for the parameters that are constrained.
  }
  \item{bound}{numeric, either a single number or a vector: the
    constraint(s) that is/are put onto the L1 norm of the parameters.}
  \item{absolute.t}{
    logical flag: if \code{TRUE}, then \code{bound} is an absolute bound and
    all entries in \code{bound} can be any positive number.  If
    \code{FALSE}, then \code{bound} is a relative bound and all entries
    must be between 0 and 1.}
}
\value{
  an object of class \code{l1ce} (if \code{bound} was a single value) or
  \code{l1celist} (if \code{bound} was a vector of values) is returned.
  See \code{\link{l1ce.object}} and \code{\link{l1celist.object}} for details.
}
\references{
  Osborne, M.R., Presnell, B. and Turlach, B.A. (2000)
  On the LASSO and its Dual,
  \emph{Journal of Computational and Graphical Statistics}
  \bold{9}(2), 319--337.

  Tibshirani, R. (1996)
  Regression shrinkage and selection via the lasso,
  \emph{Journal of the Royal Statistical Society, Series B}
  \bold{58}(1), 267--288.
}
\examples{
data(Iowa)
l1c.I <- l1ce(Yield ~ ., Iowa, bound = 10, absolute.t=TRUE)
l1c.I

## The same, printing information in each step:
l1ce(Yield ~ ., Iowa, bound = 10, trace = TRUE, absolute.t=TRUE)

data(Prostate)
l1c.P <- l1ce(lpsa ~ ., Prostate, bound=(1:30)/30)
length(l1c.P)# 30 l1ce models
l1c.P # -- MM: too large; should do this in summary(.)!
%% summary(l1c.P)
\testonly{
str(l1c.P, max.lev = 1)
}

plot(resid(l1c.I) ~ fitted(l1c.I))
abline(h = 0, lty = 3, lwd = .2)
}
\keyword{models}
\keyword{optimize}
\keyword{regression}
% Converted by Sd2Rd version 1.21.
