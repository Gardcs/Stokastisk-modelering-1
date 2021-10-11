#Defining the needed constants and lists
t <- 59
lambda <- 1.5
n <- 1000
gamma <- 10
colors <- c("green", "blue", "red", "magenta", "orange", "yellow", "indianred4", "purple", "pink", "cyan")
x0 <- 0:58
x1 <- 1:59


add0 <- function(vec){
  return(c(0, vec))
}

makeStep <- function(vec){
  sfun <- stepfun(1:59, add0(vec), f=0)
  return(sfun)
}


sumize <- function(vec){
  for(i in 1:(length(vec)-1)){
    vec[i+1] <- vec[i+1] + vec[i]
  }
  return(vec)
}

Xt <- function(t, lambda){
  vec <- rpois(t, lambda)
  return(sumize(vec))
}

checkAbove100 <- function(t, lambda, times){
  above100 <- 0
  for(i in 1:times){
    vec <- rpois(t, lambda)
    if(sum(vec) > 100){
      above100 <- above100 + 1
    }
  }
  return(above100/times)
}

XtMatrix <- function(t, lambda, times){
  vec <- matrix(rpois(t*times, lambda), t, times)
  for(i in 1:times){
    vec[, i] <- sumize(vec[, i])
  }
  return(vec)
}

#b



fexp <- function(c, gamma=10){
  return(gamma*exp(-c*gamma))
}

Zt <- function(Xt, gamma){
  Z = 1:length(Xt) * 0
  for(i in 1:length(Xt)){
    Z[i] <- sum(rexp(Xt[i], gamma))
  }
  for(i in 1:(length(Xt)-1)){
    Z[i+1] <- Z[i+1] + Z[i]
  }
  return(Z)
}

ZtMatrix <- function(t, lambda, gamma, times){
  zMatrix <- matrix(0, t, times)
  xMatrix <- matrix(rpois(t*times, lambda), t, times)
  for(i in 1:times){
    for(j in 1:t){
      zMatrix[j, i] <- sum(rexp(xMatrix[j, i], gamma))
    }
  }
  for(i in 1:times){
    zMatrix[, i] <- sumize(zMatrix[, i])
  }
  return(zMatrix)
}

plotMatrix <- function(mat, title, xlab, ylab){
  n = length(mat[1,])
  plot.new()
  plot.window(xlim=c(0, 59), ylim=c(0, max(mat)+1))
  axis(1)
  axis(2)
  box()
  title(main=title)
  title(xlab=xlab)
  title(ylab=ylab)
  for(i in 1:n){
    segments(x0, mat[, i], x1, mat[, i], lwd=2, col=colors[i])
  }
}

above8 <- function(t, lambda, gamma, n){
  above <- 0
  for(i in 1:n){
    Z <- Zt(rpois(t, lambda), gamma)
    if(Z[length(Z)] > 8){
      above <- above +1
    }
  }
  return(above/n)
}


above100 <- checkAbove100(t, lambda, 100000)
print(above100)

print(above8(t, lambda, gamma, 1000))

xMatrix <- XtMatrix(t, lambda, 10)
zMatrix <- ZtMatrix(t, lambda, gamma, 10)

plotMatrix(xMatrix, "X(t)", "Days", "Number of insurance claims")
plotMatrix(zMatrix, "Z(t)", "Days", "Value of insurance claims")








