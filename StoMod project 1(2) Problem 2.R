#Setting a seed to get the same pseudo random functiones each time
set.seed(13715)

#Setting given values
t <- 59
lambda <- 1.5
n <- 1000
gamma <- 10

#Creating lists that are used for plotting
x0 <- 0:58
x1 <- 1:59

#Creating a list of colors for plotting
colors <- c("green", "blue", "red", "magenta", "orange", "yellow", "indianred4", "purple", "pink", "cyan")

#___FUNCTIONS___

#sumize(vec)
#Takes a vector as input
#returns a vector where every entry is the sum of all entrys before in the input vector
sumize <- function(vec){
  for(i in 1:(length(vec)-1)){
    vec[i+1] <- vec[i+1] + vec[i]
  }
  return(vec)
}

#Xt(t, lambda)
#Takes an int t as length of vector and lambda as a rate in Poisson distribution
#Uses the rpois() function to generate a vector of length t, from a Poisson distribution. Then the function
#uses the sumize() function to create a cumulative poisson distributet vector and return it.
Xt <- function(t, lambda){
  vec <- rpois(t, lambda)
  return(sumize(vec))
}

#XtMatrix(t, lambda, n)
#takes an int t as length of vector, lambda as rate in Poisson distribution and n as number of vectors in
#the matrix.
#Uses the rpois() function to gnenerate a vector of length n*t from the Poisson distribution, then creates a matrix
#from this vector of n vectors of length t. This can be done this way becaouse poisson don't have a memory. For each
#of the vectors use sumize() to create cumulative poisson distributet vectors. Then return the matrix.
XtMatrix <- function(t, lambda, n){
  vec <- rpois(n*t, lambda)
  Xmatrix <- matrix(vec, t, n)
  for(i in 1:n){
    Xmatrix[, i] <- sumize(Xmatrix[, i])
  }
  return(Xmatrix)
}

#above100(t, lambda, n)
#takes an int t as length of vector, lambda as rate in Poisson distribution and n as number of vectors in
#the matrix.
#Creates n vectors with Xt() and counts how many of the have a last entry that is larger than 100. Then
#calculates the precentage of those above 100 and returns that value.
above100 <- function(t, lambda, n){
  above <- 0
  for(i in 1:n){
    vec <- Xt(t, lambda)
    if(vec[length(vec)] > 100){
      above <- above + 1
    }
  }
  return(above/n)
}

#Zt(t, lambda, gamma)
#takes an int t as length of vector, lambda as rate in Poisson distribution and gamma as rate in exponential
#distribution.
#First it uses rpois() to generate a poisson distributed vector. Then creates a vector Z of the same
#length but with zero at every entry. For each entry Z_i it calculates sum of all entries that are returned 
#from rexp(Xt_i, gamma), then uses the sumize() function and returns the vector.
Zt <- function(t, lambda, gamma){
  x <- rpois(t, lambda)
  Z <- 1:length(x) * 0
  for(i in 1:length(Z)){
    Z[i] <- sum(rexp(x[i], gamma))
  }
  return(sumize(Z))
}

#ZtMatrix(t, lambda, gamma, n)
#Takes an int t as lenth of vector, lambda as rate of Poisson distribution, gamma as rate of Exponential
#distribution and n as number of vectors in a matrix.
#First the function makes a n*t ling vector with rpois(), and make a matrix with n vectors of length t. This
#it can do because of the property that Poisson distribution has no memory. For each of the vectors in the 
#matrix the function calculates the value of claims as in Zt(). At the end the function uses sumize() on
#each vector and returns the matrix.
ZtMatrix <- function(t, lambda, gamma, n){
  vec <- rpois(n*t, lambda)
  Zmatrix <- matrix(vec, t, n)
  for(i in 1:n){
    for(j in 1:t){
      Zmatrix[j, i] <- sum(rexp(Zmatrix[j, i], gamma))
    }
    Zmatrix[, i] <- sumize(Zmatrix[, i])
  }
  return(Zmatrix)
}

#above8(t, lambda, gamma, n)
#Takes an int t as length of vector, lambda as rate of Poisson distribution, gamma as rate of Exponential
#distribution and n as number of iterations to calculate the average.
#The function creates n Zt vector and counts how many of them that have a last entry that is larger then 
#8. Then return the counter divided by n as the precentage of claims above 8 million.
above8 <- function(t, lambda, gamma, n){
  above <- 0
  for(i in 1:n){
    vec <- Zt(t, lambda, gamma)
    if(vec[length(vec)] > 8){
      above <- above +1
    }
  }
  return(above/n)
}

#plotMatrix(mat, title, xlab, ylab)
#Takes a matrix with values in y-direction to plot and three strings(title, label for x axis and label for y axis).
#The function creates a window and set parameters for size of window. Then it makes stair graphs for vectors in
#the matrix(without vertical lines).
plotMatrix <- function(mat, title, xlab, ylab){
  n = length(mat[1,])
  plot.new()
  plot.window(xlim=c(0, 59), ylim=c(0, max(mat)+(max(mat)/10)))
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


above100 <- checkAbove100(t, lambda, 10000)
print(above100)

above8 <- above8(t, lambda, gamma, 1000)
print(above8)

xMatrix <- XtMatrix(t, lambda, 10)
plotMatrix(xMatrix, "Insurance claims over a time period of 59 days", "Days", "Number of insurance claims")

zMatrix <- ZtMatrix(t, lambda, gamma, 10)
plotMatrix(zMatrix, "Total value of insurance claims over a time period of 59 days", "Days", "Value of insurance claims in millions")



#Non of the graphs in task 2 has a legend, this is because it's no point in a figure with ten graphs.
