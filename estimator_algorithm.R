## Algorithm to compute a log-concave density estimate for n observations in 
## R^d (d >= 2).

## We make the following assumptions:
## 1. The observations are independent copies of a random variable X which has 
##    mean zero and a log-concave density.
## 2. There exists an orthogonal matrix W such that WX has independent 
##    components.

## The idea underpinning the algorithm is that the density of WX factors into
## d 1D marginals by independence. This fact allows us to compute d 
## 1-dimensional log-concave MLEs instead of 1 d-dimensional log-concave MLE.
## This is much more computationally and statistically efficient.

## Suppose that X_1, ..., X_d are our observations. The steps of the algorithm
## are as follows:

## 1. For mathematical analysis reasons, and to ensure there is no "confounding"
##    between steps 2 and 3, first randomize the order of the samples, and then 
##    split them into two sets: Y_1, ..., Y_p and Z_1, ..., Z_q. The fraction
##    of samples going into the first set will be denoted by r, and is an 
##    adjustable parameter in the computation.
## 2. Compute an estimate W_hat of W by doing PCA on Y_1, ..., Y_p, and putting
##    the PCs into the rows of a matrix.
## 3. Compute the vectors W_hat*Z_1, ..., W_hat*Z_q.
## 4. For each 1 <= i <= d, take the ith components of the above vectors, and
##    use them to compute the 1D log-concave MLE f_i using logcondens package.
## 5. Let (W_hat)_i denote the ith row of W_hat. Return the density function
##    f: R^d -> R given by f(x) = f_1[(W_hat)_1 * x] ... f_d[(W_hat)_d * x] as 
##    the estimate of the true density (here * denotes the dot product).


## INSTRUCTIONS: after modifying get_data to fit the desired data/application, 
## run the following commands in the following order:

## 1. get_data()
## 2. randomize_and_split(r), where r (0 <= x <= 1) is the fraction of samples
##                            you want to use to estimate W_hat.
## 3. generate_estimator()
## 4. evaluate_estimator(x), where x is the vector in R^d at which you would 
##                           like to evaluate the density. Run this as many 
##                           times as needed after running the first 3 commands.


library(logcondens) 
library(LogConcDEAD)
library(cubature)

## Generates/gets pre-processed data to be worked with. This should be modified
## according to the desired application and data. Creates global variable 
## pre_Data, which is a matrix whose rows represent each data point.
get_data <- function() {
  library(MASS)
  num_samples = 2000
  covariance = diag(4)
  dimension = NROW(covariance)
  pre_data <<- mvrnorm(num_samples, numeric(dimension), covariance)
}


## Randomizes and splits the observations into two sets, with fraction r 
## (0 <= r <= 1) of them going into the first set. Must run 
## get_Data first or else error. Creates global variables data_1 and data_2.
randomize_and_split <- function(r) {
  dimension <<- NCOL(pre_data)
  num_samples = NROW(pre_data)
  num_samples_1 = floor(r * num_samples)
  randomized_data = pre_data[sample(num_samples),]
  data_1 <<- randomized_data[1:num_samples_1,]
  data_2 <<- randomized_data[(num_samples_1+1):num_samples,]
}


## Generates the density estimate using the algorithm above, 
## based on the two sets of data from randomize_and_split (must run 
## randomize_and_split first or else error).
generate_estimator <- function() {
  W_hat <<- prcomp(data_1)$rotation
  unmixed_obs = data_2 %*% W_hat
  marginals <<- list()
  
  for (i in 1:dimension) {
    marginals[i] <<- list(logConDens(unmixed_obs[,i], smoothed = FALSE))
  }
}


## Evaluates the estimator most recently generated in generate_estimator at the 
## vector x in R^d (must run generate_estimator beforehand or else error).
evaluate_estimator <- function(x) {
  result = 1
  for (i in 1:dimension) {
    result = result * evaluateLogConDens(W_hat[i,] %*% x, 
                                         marginals[i][[1]], 2)[3]
  }
  return(result)
}


## TEST FUNCTIONS (FOR A STANDARD GAUSSIAN):

## Generates n samples from a d-dimensional standard Gaussian and calculates
## the density estimate using the algorithm above, with sample splitting ratio
## r. If -1 is entered for err, then only the time taken to generate the 
## estimator is returned; if a nonnegative number is entered for err, then the
## squared Hellinger distance between the true Gaussian density and our 
## algorithm's density estimate is returned, with maximum tolerable absolute
## error err (if err = 0 then the best possible absolute error is obtained).
standard_gaussian_test <- function(d, n, r, err) {
  output <- matrix(data = NA, nrow = 1, ncol = 2, byrow = FALSE,
                   dimnames = list(c(""), 
                                   c("time (sec)        ", 
                                     "error (Hellinger^2)")))
  output[1,2] <- -1
  
  get_gaussian_test_data(d, n)
  
  start_time <- Sys.time()
  randomize_and_split(r)
  generate_estimator()
  end_time <- Sys.time()
  
  if (err >= 0) {
    integrand <- function(x) {
      return (0.5 * (sqrt(gaussian_test_data_pdf(x)) - 
                       sqrt(evaluate_estimator(x)))^2)
    }
    output[1,2] <- adaptIntegrate(
      integrand, rep(-Inf, dimension), rep(Inf, dimension), 
      absError = err)$integral
  }
  
  output[1, 1] <- end_time - start_time
  return(output)
} 

## Helper functions for test:

## Generates n samples from a d-dimensional standard Gaussian.
get_gaussian_test_data <- function(d, n) {
  library(MASS)
  covariance = diag(d)
  pre_data <<- mvrnorm(n, numeric(d), covariance)
}

## Evaluates the pdf of the Gaussian from get_gaussian_test_data at the point x.
gaussian_test_data_pdf <- function(x) {
  return ((2*pi)^(-0.5 * dimension) * exp(-0.5 * x %*% x))
}

















