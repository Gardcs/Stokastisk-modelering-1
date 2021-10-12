#Project 1

set.seed(7) #Set a seed, our favorite number

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
  
  #A vector of the states every year in the last 10 years, is the seconc half
  # of the prevoiusly declared vector SIR_state
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

