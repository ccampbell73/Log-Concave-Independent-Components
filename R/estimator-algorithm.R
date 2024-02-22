## Algorithm to compute a log-concave density estimate for n observations in 
## R^d.

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
##    between steps 2 and 3,first randomize the order of the samples, and then 
##    split them into two halves: Y_1, ..., Y_p and Z_1, ..., Z_q.
## 2. Compute an estimate W_hat of W by doing PCA on Y_1, ..., Y_p, and putting
##    the PCs into the rows of a matrix.
## 3. Compute the vectors W_hat*Z_1, ..., W_hat*Z_q.
## 4. For each 1 <= i <= d, take the ith components of the above vectors, and
##    use them to compute the 1D log-concave MLE f_i using logcondens package.
## 5. Let (W_hat)_i denote the ith row of W_hat. Return the density function
##    f: R^d -> R given by f(x) = f_1[(W_hat)_1 * x] ... f_d[(W_hat)_d * x] as 
##    the estimate of the true density (here * denotes the dot product).



library(logcondens) 

## Generates/gets pre-processed data to be worked with. This should be modified
## according to the desired application and data. Creates global variable 
## pre_Data, which is a matrix whose rows represent each data point.
get_data <- function() {
  library(MASS)
  num_samples = 100
  covariance = matrix(c(2, -1, 0, 
                        -1, 2, -1, 
                        0, -1, 2), 3, 3)
  dimension = NROW(covariance)
  pre_data <<- mvrnorm(num_samples, numeric(dimension), covariance)
}

## Randomizes and splits the observations into two halves. If there is an odd
## number of observations, places the extra in the second half. Must run 
## get_Data first or else error. Creates global variables data_1 and data_2.
randomize_and_split <- function() {
  dimension <<- NCOL(pre_data)
  num_samples = NROW(pre_data)
  num_samples_1 = floor(num_samples / 2)
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





