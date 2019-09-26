# Test script

MAE <- function(x,y)
{
  mean(abs(x-y))/2
}

concordance <- function(x,y)
{
  MAE_max = 1/6
  1-MAE(x,y)/MAE_max
}

x <- runif(1000)
y <- runif(1000)

concordance(x,y)


a <- density(rnorm(10000))
plot(a)
polygon(a, col = "blue")











