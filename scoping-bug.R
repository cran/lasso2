##- Date: Wed, 9 Jul 2003 14:47:13 -0400

##- From: Jim_Garrett@bd.com
##- To: Martin Maechler <maechler@stat.math.ethz.ch>
##- Subject: Re: lasso2


##- Since you sound so encouraging ;) here's another interesting feature I
##- discovered:


### ---> see also ./predict-new-bug.R
###		    ~~~~~~~~~~~~~~~~~

# same data as before
set.seed(123)
d3 <- data.frame(x1 = rnorm(200), x2 = rnorm(200), x3 = rnorm(200))
d3$yb <- factor(rbinom(200, 1, prob = 1 - 1/(1 + exp(d3$x1))))
d3$y  <- rnorm(200, mean = 0.5 * d3$x1, sd = 0.5)

## Ok:
gl1ce(yb ~ x1 + x2 + x3, data = d3, family = binomial())
# Now call it inside a function:
temp.fun <- function(InputData) {
  gl1ce(yb ~ x1 + x2 + x3, data = InputData, family = binomial())
}
temp.fun(d3)
## Error in model.frame.default(formula = y ~ x1 + x2 + x3 + 1, data= InputData) :
##    Object "InputData" not found

## --------- MM : -------------

## The same for simple glm()
glm(yb ~ x1 + x2 + x3, data = d3, family = binomial())
# Now call it inside a function:
myglm <- function(InputData) {
  glm(yb ~ x1 + x2 + x3, data = InputData, family = binomial())
}
myglm(d3)
## works fine

## Also with l1ce():
## Ok:
l1ce(y ~ x1 + x2 + x3, data = d3)
# Now call it inside a function:
myl1 <- function(InputData) l1ce(y ~ x1 + x2 + x3, data = InputData)
myl1(d3)
##-> Error

## The same for simple lm()
lm(y ~ x1 + x2 + x3, data = d3)
# Now call it inside a function:
mylm <- function(InputData) lm(y ~ x1 + x2 + x3, data = InputData)
mylm(d3)
## works fine


## Note, using
trace(model.frame)
##  I see that  l1ce() and gl1ce() call  model.frame()  as
model.frame(formula = y ~ x1 + x2 + x3 + 1, data = InputData)
##				      ====
## whereas in myglm() or mylm() it's
model.frame(formula = y ~ x1 + x2 + x3, data = InputData, drop.unused.levels = TRUE)
