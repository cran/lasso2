\name{gl1ce}
\alias{gl1ce}
\alias{family.gl1ce}
\title{Generalized Regression With L1-constraint on the Parameters}
\description{
  Fit a generalized regression problem while imposing an L1 constraint on
  the parameters.  Returns an object of class \code{gl1ce}.
}
\usage{
gl1ce(formula, data = parent.frame(), weights, subset, na.action,
      family = gaussian, control = glm.control(\dots), sweep.out = ~ 1,
      x = FALSE, y = TRUE, contrasts = NULL, standardize = TRUE,
      guess.constrained.coefficients = double(p), bound = 0.5, \dots)
\method{family}{gl1ce}(object, \dots)
}
\arguments{
  \item{formula}{a \code{\link{formula}}, with the response on the left
    hand side of a \code{~} operator, and the terms, separated by a
    \code{+} operator, on the right hand side.}
  \item{data}{a \code{\link{data.frame}} in which to interpret the
    variables named in the formula, the \code{weights}, the
    \code{subset} and the \code{sweep.out} argument.  If this is
    missing, then the variables in the formula should be globally
    available.}
  \item{weights}{vector of observation weights.  The length of
    \code{weights} must be the same as the number of observations.  The
    weights must be strictly positive, since zero weights are ambiguous,
    compared to use of the \code{subset} argument.}
  \item{subset}{
    expression saying which subset of the rows of the data should be used in the
    fit.  This can be a logical vector (which is replicated to have length equal
    to the number of observations), or a numeric vector indicating which
    observation numbers are to be included, or a character vector of the
    row names to be included.  All observations are included by default.
  }
  \item{na.action}{
    a function to be applied to the \code{model.frame} after any
    \code{subset} argument has been used.  The default (with
    \code{na.fail}) is to create an error if any missing values are
    found.  A possible alternative is \code{na.omit}, which deletes
    observations that contain one or more missing values.
  }
  \item{family}{a \code{\link{family}} object - a list of functions and
    expressions for defining the link and variance functions,
    initialization and iterative weights.  Families supported are
    gaussian, binomial, poisson, Gamma, inverse.gaussian and quasi.
    Functions like binomial produce a family object, but can be given
    without the parentheses.  Family functions can take arguments, as in
    \code{binomial(link=probit)}.
  }
  \item{control}{
    a list of iteration and algorithmic  constants.  See glm.control for their
    names and default values.  These can also be set as arguments to gl1ce itself.
  }
  \item{sweep.out}{
    a formula object, variables whose parameters are not put under the
    constraint are swept out first.  The variables should appear on the
    right of a \code{~} operator and be separated by \code{+} operators.  Default is
    \code{~1}, i.e. the constant term is not under the constraint.  If this
    parameter is \code{NULL}, then all parameters are put under the constraint.
  }
  \item{x}{
    logical flag: if \code{TRUE}, the model matrix is returned in component \code{x}.
  }
  \item{y}{
    logical flag: if \code{TRUE}, the response is returned in component \code{y}.
  }
  \item{contrasts}{
    a list giving contrasts for some or all of the factors appearing in the model
    formula. The elements of the list should have the same name as the variable
    and should be either a contrast matrix (specifically, any full-rank
    matrix with as many rows as there are levels in the factor),
    or else a function to compute such a matrix given the number of levels.
  }
  \item{standardize}{
    logical flag: if \code{TRUE}, then the columns of the model matrix that
    correspond to parameters that are constrained are standardized to have
    empirical variance one.  The standardization is done after taking
    possible weights into account and after sweeping out variables whose
    parameters are not constrained.
  }
  \item{guess.constrained.coefficients}{
    initial guess for the parameters that are constrained.
  }
  \item{bound}{numeric, either a single number or a vector: the
    constraint(s) that is/are put onto the L1 norm of the parameters.}
  \item{\dots}{potential arguments for \code{\link[stats]{glm.control}},
    as default for the \code{control} argument above.}
  \item{object}{an \R object of class \code{"gl1ce"}.}
}% args
\value{
  an object of class \code{gl1ce} is returned by \code{gl1ce()}.
  See \code{\link{gl1ce.object}} for details.
}
\references{
  See the references in \code{\link{l1ce}}.

  Justin Lokhorst (1999).
  The LASSO and Generalised Linear Models,
  Honors Project, Nov. 1999, Dept.Statist., Univ. of Adelaide.
}
\seealso{\code{\link{glm}} for unconstrained generalized regression
  modeling.
}
\examples{
## example from base:
data(esoph)
summary(esoph)
## effects of alcohol, tobacco and interaction, age-adjusted
modEso <- formula(cbind(ncases, ncontrols) ~ agegp + tobgp * alcgp)
glm.E   <- glm(modEso, data = esoph, family = binomial())
gl1c.E <- gl1ce(modEso, data = esoph, family = binomial())
gl1c.E
plot(residuals(gl1c.E) ~ fitted(gl1c.E))

sg1c <- summary(gl1c.E)
sg1c

## Another comparison  glm() / gl1c.E:
plot(predict(glm.E,  type="link"), predict(glm.E,  type="response"),
     xlim = c(-3,0))
points(predict(gl1c.E, type="link"), predict(gl1c.E, type="response"),
       col = 2, cex = 1.5)

%%% mabye FIXME!!
labels(gl1c.E)#-- oops! empty!!
}
\keyword{models}
\keyword{optimize}
\keyword{regression}
% Converted by Sd2Rd version 1.21.
