##- Date: Mon, 7 Jul 2003 14:39:49 -0400
##- From: Jim_Garrett@bd.com
##- To: maechler@stat.math.ethz.ch
##- Subject: lasso2

## I started doing some work using the lasso2 package, and I've had some
## trouble using predict.gl1ce and predict.l1ce on new data.

###--- MM: heavily edited this and added  -------

set.seed(123)
d3 <- data.frame(x1 = rnorm(200), x2 = rnorm(200), x3 = rnorm(200))
d3$yb <- factor(rbinom(200, 1, prob = 1 - 1/(1 + exp(d3$x1))))
d3$y  <- rnorm(200, mean = 0.5 * d3$x1, sd = 0.5)

ndat <- d3[ 1:2, 1:3]

l1. <- l1ce (y  ~ x1 + x2 + x3, data = d3)
## The following  *FAILS* in R <= 1.3.1 (!) in model.matrix.default() ...
## but seems ok from 1.4.x {maybe I have "upgraded" it in that time}
gl1 <- gl1ce(yb ~ x1 + x2 + x3, data = d3, family = binomial())

## predict.l1ce ()

try( predict(l1., newdata = ndat) )

##
## This returns the error message

## >>   Error in eval(expr, envir, enclos) : couldn't find function "y"

##- I've gotten this result on a (Mandrake 9.1) Linux Intel platform, R 1.7.0,
##- and Windows2000, also R 1.7.0.

## predict.gl1ce() actually calls  predict.l1ce() eventually
## {using NextMethod()!}
predict(gl1, newdata = ndat)
traceback()


##- I see the description in the documentation for gl1ce is "Extracts the
##- fitted values from a gl1ce() object and returns a matrix of predictions"
##- which implies that perhaps predicting at new points is not yet implemented.
##- On the other hand there is a "newdata" argument.  Am I right in guessing
##- that predicting with new data simply isn't implemented yet?  At any rate,
##- are you aware of this behavior?

##- I've poked around a bit but not too deeply.  If I write a fix for this I
##- will of course send it.  Thanks for porting this to R!
