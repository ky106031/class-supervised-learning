##### Lasso #####
library(glmnet) 

crime <- read.table("crime.txt")  # 犯罪データ読み込み
crime <- as.matrix(crime)
X <- crime[, 3:7]  # 説明変数
y <- crime[, 1]    # 目的変数

# 最小二乗推定
res_ls <- lm(y~X)
res_ls$coefficients #係数の推定値

# Ridge推定 (正則化パラメータを20に固定)
res_rid <- glmnet(x=X, y=y, alpha=0, lambda=20)
res_rid$beta #係数の推定値（切片除く，切片はres_rid$a0）

# Lasso推定
res_las <- glmnet(x=X, y=y)
# 解パス図描画
plot(res_las, xvar="lambda", label=TRUE, xlab="正則化パラメータの対数値", 
     ylab="回帰係数", lwd=2)
# 正則化パラメータの値を20と固定
res_las <- glmnet(x=X, y=y, lambda=20)  
res_las$beta  # 係数の推定値（切片除く，切片はres_las$a0）

# CVの計算
res.cv <- cv.glmnet(x=X, y=y)
# CV値の推移をプロット
plot(res.cv, xlab="正則化パラメータの対数値", ylab="２乗誤差")
# CV値が最小となる正則化パラメータ値を出力
res.cv$lambda.min
# 1標準誤差ルールにより選択された正則化パラメータの値を出力
res.cv$lambda.1se
# CVで選択された正則化パラメータの値でlasso実行
res2 <- glmnet(x=X, y=y, lambda=res.cv$lambda.min)  # CVmin
res2$beta  # 係数の推定値
res3 <- glmnet(x=X, y=y, lambda=res.cv$lambda.1se)  # CVmin
res3$beta  # 係数の推定値



###### Elastic net ######
# 人工データの発生
z1 <- runif(100, 0, 20)
z2 <- runif(100, 0, 20)
x <- scale(cbind(z1 + rnorm(100, 0, 0.1), -z1 + rnorm(100, 0, 0.1), z1 + rnorm(100, 0, 0.1), z2 + rnorm(100, 0, 0.1), -z2 + rnorm(100, 0, 0.1),z2 + rnorm(100, 0, 0.1)))
y <- 3*apply(x[ ,1:3], 1, sum) + apply(x[ ,4:6], 1, sum) + rnorm(100)
# エラスティックネット
elasticnet_fit <- glmnet(x, y, alpha = 0.5)
# エラスティックネットの解パス図の描画
plot(elasticnet_fit, xvar="lambda", label=TRUE, xlab="正則化パラメータの対数値", ylab="回帰係数", col="black", lwd=2.5)
# lasso
lasso_fit <- glmnet(x, y)
# lasso の解パス図の描画
plot(lasso_fit, xvar="lambda", label=TRUE, xlab="正則化パラメータの対数値", ylab="回帰係数", col="black", lwd=2.5)



##### fused lasso ######
### CGHデータ解析
library(genlasso)
y <- read.table("cgh.txt")
y <- unlist(y)            # 目的変数

# FLSA適用
res <- fusedlasso1d(y=y)
# FLSAによる当てはめ結果
plot(res, lambda=3, pch=20, col="gray", 
     xlab="遺伝子番号", ylab="コピー数の比の対数値",
     cex.lab=1.5, cex.axis=1.5)

## Fit Trend filtering
# trendfilterの引数ordの値を2,3,...とすることで
# 区分2,3,...次多項式関数を当てはめることができる
set.seed(1)
n = 100
y = 50 - abs(1:n-n/2) + rnorm(n, sd=5)  #データ生成
#トレンドフィルタリング実行
out = trendfilter(y, ord=1)  
plot(out, lambda=100) #結果plot


##### Group lasso ######
library(grpreg)
data(Birthwt)     # データセット

X <- Birthwt$X    # 説明変数
y <- Birthwt$bwt  # 目的変数
group <- Birthwt$group   # グループ変数のインデックス

# グループlassoの実行
res <- grpreg(X, y, group, penalty="grLasso")

# 解パス図の作成
par(mar=c(5, 5, 1, 1))
plot(0, 0, xlim=c(-0.01,0.15), ylim=range(res$beta[-1,]), 
     type="n", xlab="正則化パラメータの値", ylab="回帰係数")
for(i in c(1:3, 9, 14:16)){
  lines(res$lambda, res$beta[i+1, ], lty=as.numeric(group)[i],
        lwd=2)
  text(-0.01, res$beta[i+1,100], i, cex=1)
}
# 解パス図は次のプログラムでも作成可能
plot(res, lty=as.numeric(group), col=as.numeric(group), 
     xlim=c(0, 0.2), xlab="正則化パラメータの値", ylab="回帰係数")
res$beta[,10] # 回帰係数の推定値

##### Graphical model #####
dat.deca <- read.csv("decathlon.csv", header = T, row.names = 1)  # データの読み込み
cordat <- cor(dat.deca)

# グラフィカルlasso を当てはめる
library(glasso)  # グラフィカルlasso を実行するのに必要なパッケージの読み込み
library(igraph)  # グラフを描くのに必要なパッケージの読み込み
fit.glasso <- glasso(cordat, rho = 0.4)  # λ = 0.4 のときの結果
fit.glasso$wi  # 共分散逆行列の出力

# 隣接行列を作成する関数
makeadjmat <- function(mat, cordat) {
  ans <- mat
  diag(ans) <- 0
  ans <- ans != 0
  ans <- matrix(as.numeric(ans), nrow(ans), nrow(ans))
  colnames(ans) <- rownames(ans) <- colnames(cordat)
  ans
}
adjmat <- makeadjmat(fit.glasso$wi, cordat)  # 隣接行列を作成
g <- graph.adjacency(adjmat, mode = "undirected")  # グラフを描くためのオブジェクトの作成
plot(g, vertex.size = 40, vertex.shape = "rectangle", vertex.color = "#FFFF99")
