####------ l1ce() tests ------------
library(lasso2)

###-------- +/- the same as  example(l1ce) :
data(Iowa)

l1c.I <- l1ce(Yield ~ ., Iowa, bound = 10, trace = TRUE, absolute.t=TRUE)

## the next ones give a l1ce LIST (one for each bound)
l1c.liI <- l1ce(Yield ~ ., Iowa, bound = seq(0,1, len= 17))
print(plot(l1c.liI))

data(Prostate)
l1c.P <- l1ce(lpsa ~ ., Prostate, bound= seq(0,1, len= 17))
print(plot(l1c.P))

###-------- Try a case where p > n :
n <- 100
p <- 120

RNGversion("1.6.0")
set.seed(n+p)
x <- matrix(runif(n*p), n,p, dimnames= list(NULL, paste("x",1:p,sep="")))
x <- as.data.frame(apply(x, 2, sort))
with(x,
     y <<- 5 + 4*x1 -3*x2 + 10*x3  + (eps <<- rnorm(n, sd = 1/200)))
d.ex <- cbind(y = y, x)
dim(d.ex)# 100 x 121
if(FALSE)
summary(lm(y ~ ., data = d.ex))# pretty nonsense

## gives something, but not at all the true model ...
l20 <- l1ce(y ~ ., data = d.ex, bound = 20, absolute.t = TRUE)
coef(l20)[coef(l20) > 0]

l15 <- l1ce(y ~ ., data = d.ex, bound = 15, absolute.t = TRUE)
coef(l15)[coef(l15) > 0]

sum(eps^2)
sum(resid(l20)^2) / sum(eps^2)
sum(resid(l15)^2) / sum(eps^2)

l1.lis <- l1ce(y ~ ., data = d.ex, bound = seq(0, 0.1, len=21))
pl1lis <- plot(l1.lis, ylim = c(-10,10))
round(pl1lis$mat,3)
bnds <- pl1lis$bounds[,"rel.bound"]
round(1000 * t(pl1lis$mat[0.01 < bnds & bnds < 0.06 ,]))
