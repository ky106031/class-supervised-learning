library(MASS)
library(GPfit)
library(lhs)
library(rgl)

### Gaussian process
n = 51
# 説明変数に対応するデータ
x = seq(0,1,length=n)

plot(0,0, type="n", xlim = c(0,1), ylim=c(-3,3), xlab="", ylab="")
for(k in 1:10){
  #カーネル関数からなる行列（グラム行列）
  K = matrix(nr=n, nc=n) 
  for(i in 1:n){
    for(j in 1:n){
      K[i,j] = exp(-(x[i]-x[j])^2/0.1) #ガウスカーネルを使用
    }
  }
  #グラム行列を分散共分散行列とした正規乱数を生成
#  y = mvrnorm(1, rep(0,n), K)
  y = mvrnorm(1, rep(0,n), K+diag(0.01, n))
  lines(x,y, col=k)
}



### GPfit
## 1D Example 2
n = 10  #時点数
d = 1  #次元
# 真の関数
computer_simulator <- function(x) {
  y <- log(x + 0.1) + sin(5 * pi * x)
  return(y)
}
# データ生成
x = maximinLHS(n, d)
y = computer_simulator(x) + rnorm(n, 0, 0.3)
## GP regression
GPmodel = GP_fit(x, y)　# 0<x<1
print(GPmodel, digits = 4)
plot(GPmodel)

# 新たな観測点(x)における予測
xnew <- seq(from = 0, to = 1, length.out = 21)
predres <- predict(GPmodel, xnew)



## 2D Example: GoldPrice Function
computer_simulator <- function(x) {
  x1 = 1 * x[, 1] - 2
  x2 = 1 * x[, 2] - 2
  t1 = 1 + (x1 + x2 + 1)^2 * (19 - 14 * x1 + 3 * x1^2 - 14 * x2 +
                                6 * x1 *x2 + 3 * x2^2)
  t2 = 30 + (2 * x1 - 3 * x2)^2 * (18 - 32 * x1 + 12 * x1^2 + 48 * x2 -
                                     36 * x1 * x2 + 27 * x2^2)
  y = t1 * t2
  return(y)
}
n = 30
d = 2
set.seed(100)
x = maximinLHS(n, d)
y = computer_simulator(x)
plot3d(x[,1],x[,2],y)

GPmodel = GP_fit(x, y)
print(GPmodel)
plot(GPmodel)

xnew1 <- seq(from = 0, to = 1, length.out = 21)
xnew = expand.grid(xnew1, xnew1)
ypred = predict(GPmodel, xnew)



### Bayesian Optimization
n = 8
d = 1
computer_simulator <- function(x) {
  y <- log(x + 0.1) + sin(5 * pi * x)
  return(y)
}
x = maximinLHS(n, d)
y = computer_simulator(x) 
## GP regression
GPmodel = GP_fit(x, y)
print(GPmodel, digits = 4)
plot(GPmodel)

#Calcuate acquisition function
x_new <- seq(0, 1, length.out = 100)
pred <- predict.GP(GPmodel, xnew = data.frame(x = x_new))
mu <- pred$Y_hat
sigma <- sqrt(pred$MSE)
acquisition = mu + sqrt(log(n)/n)*sigma
lines(x_new,acquisition, col=5)
t = (mu - max(GPmodel$Y))/sigma
acquisition2 = (mu - max(GPmodel$Y))*pnorm(t) + sigma*dnorm(t)
lines(x_new,acquisition2*20-3, col=5)

