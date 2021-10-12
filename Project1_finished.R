#Project 1

set.seed(7) #Set a seed, our favorite number. So we get the same pseudo random numbers each run

#Problem 1

#1c)

days = 7300 #20 years in days
P = matrix(c(0.99, 0.01, 0.0, 0, 0.9, 0.1, 0.005, 0, 0.995), byrow = TRUE, nrow = 3) #transition probability matrix

#MCSimulation is a function that simulates the Markov chain in problem 1.1. 
MCSimulation <- function(){
  #Declare a vector, that holds a state in [0, 1, 2] for each day the simulation runs
  #which in this case, is 7300 days = 20 years
  SIR_state = vector("numeric", length = days) 
  SIR_state[1] = 0 #Set the first day to state 0, which is the susceptible state
  
  #Do the simulation
  for(i in 1:(days-1)){
    previousDay = SIR_state[i] #state for day i
    todayProb = P[previousDay + 1, ] #the probabilities for which state you end in
    # in day i+1 is the row SIR_state[i] in the transition
    # probability matrix. The +1 is due to indexing in
    # R matrix. 
    # Sample state of next day from the vector [0, 1, 2] with probabilities todayProb
    SIR_state[i+1] = sample(c(0, 1, 2), size = 1, prob = todayProb)
  }
  
  #A vector of the states every year in the last 10 years, is the second half
  # of the previously declared vector SIR_state
  lastTenYears <-  SIR_state[(length(SIR_state)/2):days] 
  
  #For each run of this function, we estimate the number of days per year spent in each
  #state. This is just the days in one given state, divided by all days, times 365 (1 year)
  estimateSusceptibleDays = sum(lastTenYears == 0)/length(lastTenYears) * 365
  estimateInfectedDays = sum(lastTenYears == 1)/length(lastTenYears) * 365
  estimateImmuneDays = sum(lastTenYears == 2)/length(lastTenYears) * 365
  
  #returns a vector of the days per year estimates of each state
  return (c(estimateSusceptibleDays, estimateInfectedDays, estimateImmuneDays))
  
}

#CalculateCI is a function that calculates CI for the long-run mean numbers of 
#days spent in each state every year
calculateCI <- function(){
  #Run the function MCSimulation that runs the given Markov chain for 20 years, 30 times
  #The vector estimations is a vector of vectors which stores the estimates for days per year
  #spent in each state, for all 30 runs.
  estimations <- replicate(30, MCSimulation())
  
  #Stores number of susceptible days per year for each of the 30 simulations
  susceptibleDays = estimations[1,]
  #Stores number of infected days per year for each of the 30 simulations
  infectedDays = estimations[2,]
  #Stores number of immune days per year for each of the 30 simulations
  immuneDays = estimations[3,]
  
  #nSample is numbers of samples used in calculating CI. We have one sample for 
  #each run of the MC, so nSample = 30. All the three above defined vectors
  #have the same length. 
  nSample = length(susceptibleDays)
  
  #Calculate the mean of all three samples
  meanSusceptibleDays = mean(susceptibleDays)
  meanInfectedDays = mean(infectedDays)
  meanImmuneDays = mean(immuneDays)
  
  print(meanSusceptibleDays)
  print(meanInfectedDays)
  print(meanImmuneDays)
  
  #Calculate the standard deviation of all three samples. This is the corrected
  #sample standard deviation 
  sdSusceptibleDays = sd(susceptibleDays)
  sdInfectedDays = sd(infectedDays)
  sdImmuneDays = sd(immuneDays)
  
  #Critical value in t-distribution for 29 degrees of freedom and alpha = 0.25
  #t-distribution is symmetric, so the lower and upper quantile are equal
  quantile = 2.045
  
  #Calculate lower and upper bound of each of the CI. Describtion in report
  lowerBoundS = meanSusceptibleDays-(quantile*(sdSusceptibleDays/sqrt(nSample)))
  upperBoundS = meanSusceptibleDays+(quantile*(sdSusceptibleDays/sqrt(nSample)))
  lowerBoundI = meanInfectedDays-(quantile*(sdInfectedDays/sqrt(nSample)))
  upperBoundI = meanInfectedDays+(quantile*(sdInfectedDays/sqrt(nSample)))
  lowerBoundR = meanImmuneDays-(quantile*(sdImmuneDays/sqrt(nSample)))
  upperBoundR = meanImmuneDays+(quantile*(sdImmuneDays/sqrt(nSample)))
  
  #Last, print out the CI's
  
  print("Confidence interval for mean number of days spent in susceptible state: ")
  print(c('[', lowerBoundS, upperBoundS, ']'))
  
  print("Confidence interval for mean number of days spent in infected state: ")
  print(c('[', lowerBoundI, upperBoundI, ']'))
  
  print("Confidence interval for mean number of days spent in recovered state: ")
  print(c('[', lowerBoundR, upperBoundR, ']'))
}

calculateCI()

#1e)

n = 300 #steps
N = 1000 #number of people

#The function explosiveBehaviour has two input values. This is a boolean statement
#that tells us to return the matrix Y or not. The matrix Y has 300 (= steps) columns
#and 3 rows. For each column, the first value is the number of susceptible people this day,
#the second row is number of infected people this day, and the last row is number of
#immune people this day. The second input, is the number of vaccinated individuals.
#If the input value returnY is false, the function returns a vector with 2 elements. 
#
#For each step, we can model new infected, new recovered and new susceptible people as
#a binomial distribution. We have the number of people in each state at day i, and we also
#know the probability that any of these will change state. For example, for new infected, we
#have the numbers of susceptible people at day i as numbers of trials, and beta as probability
#of success, where "success" is that a person becomes infected. We do it similarly for the other states. 
explosiveBehaviour <- function(returnY, vaccinated){
  Y <- matrix(0, 3, n) #Create 3x300 matrix of zeros
  Y[,1] = c(950, 50, 0) #Set the first columns equal to Y_0
  Y[1, 1] = Y[1, 1] - vaccinated #If any vaccinated, the number of susceptible 
  #people is vaccinated people fewer than 950
  Y[3, 1] = Y[3, 1] + vaccinated #If any vaccinated, they are also immune
  gamma = 0.10    #probability of getting recovered if an individual is infected
  alpha = 0.005   #probability of becoming susceptible if you are immune
  
  #Simulate the first 300 days of infection
  for (i in 1:(n-1)){
    S_n <- Y[1, i] #number of susceptible at day i
    I_n <- Y[2, i] #number of infected on day i
    R_n <- Y[3, i] #number of immune on day i
    beta <- 0.5 * I_n / N #Calculate the probability of becoming infected
    newInfected <- rbinom(1, S_n, beta) #new infected on day i comes from a binomial
    #distr with S_n trials and beta prob of success
    newImmune <- rbinom(1, I_n, gamma) #new immune on day i comes from a binomial distr
    #with I_n trials and prob gamma of success
    newSusceptible <- rbinom(1, (R_n-vaccinated), alpha) #new susceptible on day i 
    #comes from binomial distr with R_n - vaccinated
    #trials and prob alpha of success. The vaccinated
    #cant become susceptible again.
    infected <- I_n + newInfected - newImmune #number of infected on day i+1 is 
    #number of infected on day i + new infected - new immune
    immune <- R_n + newImmune - newSusceptible #immune on day i+1 is immune on day i
    #+ new immune - new susceptible
    susceptible <- S_n + newSusceptible - newInfected #susceptible on day i+1 is susceptible on day i
    #+ new susceptible - new infected
    
    Y[,i+1] <- c(susceptible, infected, immune) #add a column to the matrix with
    #number of people in each state at day i+1
  }
  numbersOfInfectedPerDay = Y[2, ] #Vector of 300 elements, where the i'th element is
  #number of infected people on day i
  if(returnY){
    return(Y) #returns the matrix Y
  }
  #If returnY is false, return a vector of max infected a day, and which day that was
  return ( c(max(numbersOfInfectedPerDay), which.max(numbersOfInfectedPerDay)) )
}

#Run the above function, and plot the outlook of the measles outbreak (with 0 vaccinated)
Y = explosiveBehaviour(TRUE, 0)
plot(1:n, Y[1, ], col = "blue", type = "l", lwd = 2, main = "Simulation of a measles outbreak for 300 days", 
     xlab = "Day", ylab = "Number of people")
lines(1:n, Y[2, ], col = "red", lwd = 2)
lines(1:n, Y[3, ], col = "green", lwd = 2)
abline(v = 50)
legend(x = 200, y = 600, legend=c("Susceptible", "Infected", "Recovered"), 
       col = c("blue", "red", "green"), lty = 1)

#Run the simulation of the measles outbreak 1000 times with 0 vaccinated individuals, 
#store the number of max infected people for every simulation in the vector maxValsInfected0
#and the days that this happened in dayMaXInfected0. 
maxValsInfected0 <- replicate(1000, explosiveBehaviour(FALSE, 0))
maxInfected0 = maxValsInfected0[1, ]
dayMaxInfected0 = maxValsInfected0[2, ]

#This function calculates CI's for max number of infected per day in the 300 day
#simulation and at which day this happened. The function has 2 input values. 
#Both of the input values are from the simulation with returnY=FALSE. 
calculateCI_expectedValues <- function(maxInfected, dayMaxInfected){
  #Calculate the mean of bot input parameters
  meanMaxInfected = mean(maxInfected)
  meanDay = mean(dayMaxInfected)
  
  print("E[maxI]")
  print(meanMaxInfected)
  print("E[day]")
  print(meanDay)
  
  #Calculate corrected standard deviation of both input parameters
  sdMaxInfected = sd(maxInfected)
  sdDay = sd(dayMaxInfected)
  
  #nSample = 1000 = length of one of the input parameter
  nSample = length(maxInfected)
  quantile = 2.045 #Find this in blue book
  
  #Calculate the bounds for the CI's
  lowerBoundMaxInfected = meanMaxInfected-(quantile*(sdMaxInfected/sqrt(nSample)))
  upperBoundMaxInfected = meanMaxInfected+(quantile*(sdMaxInfected/sqrt(nSample)))
  lowerBoundDay = meanDay-(quantile*(sdDay/sqrt(nSample)))
  upperBoundDay = meanDay+(quantile*(sdDay/sqrt(nSample)))
  
  
  #Print out both of the CI's
  print("Confidence interval for maxInfected: ")
  print(c('[', lowerBoundMaxInfected, upperBoundMaxInfected, ']'))
  
  print("Confidence interval for dayMaxInfected: ")
  print(c('[', lowerBoundDay, upperBoundDay, ']'))
}

#Run the simulation of the measles outbreak again, with 100, 600 and 800
#vaccinated people respectively. I then plot the number of infected people
#per day, for all three cases and the case for 0 vaccinated
Y_100 = explosiveBehaviour(TRUE, 100)
Y_600 = explosiveBehaviour(TRUE, 600)
Y_800 = explosiveBehaviour(TRUE, 800)
plot(1:n, Y[2, ], col = "red", lwd = 2, type = "l", xlab = "Days", 
     ylab = "Infected people", main = "Infected people in a 300 day simulation of a measles outbreak")
lines(1:n, Y_100[2, ], col = "blue", lwd = 2)
lines(1:n, Y_600[2, ], col = "green", lwd = 2)
lines(1:n, Y_800[2, ], col = "orange", lwd = 2)
legend(x = 200, y = 500, legend=c("0 vaccinated", "100 vaccinated", "600 vaccinated", "800 vaccinated"), 
       col = c("red", "blue", "green", "orange"), lty = 1)

#Run the simulation of the outbreak 1000 times for each of the cases where
#the number of vaccinated people is 100, 600 and 800. Store the values in a
#matrix
maxValsInfected100 <- replicate(1000, explosiveBehaviour(FALSE, 100))
maxValsInfected600 <- replicate(1000, explosiveBehaviour(FALSE, 600))
maxValsInfected800 <- replicate(1000, explosiveBehaviour(FALSE, 800))

#Split up the above matrices to get vectors that have the max number of infected
#people for each of the cases, and vectors that store the day that this happens.
maxInfected100 = maxValsInfected100[1, ]
dayMaxInfected100 = maxValsInfected100[2, ]
maxInfected600 = maxValsInfected600[1, ]
dayMaxInfected600 = maxValsInfected600[2, ]
maxInfected800 = maxValsInfected800[1, ]
dayMaxInfected800 = maxValsInfected800[2, ]

#Calculate CI's for max number of infected per day in the 300 day simulation
#and at which day this happened for the cases with 0, 100, 600 and 800 vaccinated
#people at day 1.
calculateCI_expectedValues(maxInfected0, dayMaxInfected0)
calculateCI_expectedValues(maxInfected100, dayMaxInfected100)
calculateCI_expectedValues(maxInfected600, dayMaxInfected600)
calculateCI_expectedValues(maxInfected800, dayMaxInfected800)


#Problem 2

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
#uses the sumize() function to create a cumulative poisson distributed vector and return it.
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
#Creates n vectors with Xt() and counts how many of them that have a last entry that is larger than 100. Then
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
#First the function makes a n*t long vector with rpois(), and make a matrix with n vectors of length t. This
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


above100 <- above100(t, lambda, 10000)
print(above100)

above8 <- above8(t, lambda, gamma, 1000)
print(above8)

xMatrix <- XtMatrix(t, lambda, 10)
plotMatrix(xMatrix, "Insurance claims over a time period of 59 days", "Days", "Number of insurance claims")

zMatrix <- ZtMatrix(t, lambda, gamma, 10)
plotMatrix(zMatrix, "Total value of insurance claims over a time period of 59 days", "Days", "Value of insurance claims in millions")



#Non of the graphs in task 2 has a legend, this is because it's no point in a figure with ten graphs.