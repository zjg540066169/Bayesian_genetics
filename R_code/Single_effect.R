# Compute useful posterior statistics in a simple "single effect"
# Bayesian logistic regression, in which y ~ sigmoid(b0 + x*b) and b ~
# N(mu0,s0), where b0 is the intercept, b is the regression
# coefficient, and mu0, s0 are the the prior mean and standard
# deviation of b. (The intercept is assigned a "flat" prior.)
#
# Input X is the n x p "data matrix", where n is the number of samples
# and p is the number of candidate predictors, and input y is the
# vector of n observed outcomes.
#
# Input argument p1 specifies the prior inclusion probabilities; it
# should be a vector of length p in which p1[i]/sum(p1) gives the
# prior probability that the jth candidate predictor has a nonzero
# effect on the outcome, Y.
#
# Note that these computations could probably be made faster and more
# accurate using a fast 2-d numerical integration (quadrature) method.
bayes.logistic <- function (X, y, p0, mu0, s0) {
  # These two variables define the 2-d grid used to compute the Monte
  # Carlo (importance sampling) estimates.
  b0  <- seq(-10,10,0.1) 
  b   <- seq(-10,10,0.1)
  # Create the 2-d grid.
  out <- expand.grid(b0 = b0,b = b)
  b0  <- out$b0
  b   <- out$b
  rm(out)
  # Get the number of candidate predictors (p) and importance weights (n).
  p <- ncol(X)
  n <- length(b)
  # Initialize storage for the marginal log-likelihoods (logZ) and the
  # posterior means (mu1) and standard deviations (s1) of the
  # coefficients.
  logZ <- rep(0,p)
  mu1  <- rep(0,p)
  s1   <- rep(0,p)
  # Repeat for each candidate predictor.
  for (j in 1:p)  {
    # Compute the log-importance weights, ignoring constant terms. The
    # important weight is a product of the (logistic) likelihood and
    # the (normal) prior. This is the step that requires the most
    # effort.
    logw <- rep(0,n)
    x    <- X[,j]
    for (i in 1:n) {
      u       <- b0[i] + x*b[i]
      logw[i] <- sum((y - 1)*u + logsigmoid(u)) - ((b[i] - mu0)/s0)^2/2
    }
    # Compute the importance sampling estimate of the marginal
    # log-likelihood (up to a constant of proportionality).
    u       <- max(logw)
    logZ[j] <- log(mean(exp(logw - u))) + u
    # Compute the normalized importance weights.
    w <- softmax(logw)
    # Compute the mean and standard deviation of the coefficient.
    mu1[j] <- sum(w*b)
    s1[j]  <- sqrt(sum(w*(b - mu1[j])^2))
  }
  # Compute the posterior inclusion probabilities.
  p1 <- softmax(logZ + log(p0))
  # Output the data frame containing the computed posterior
  # statistics.
  return(data.frame(p1 = p1,mu1 = mu1,s1 = s1))
}
# sigmoid(x) returns the sigmoid of the elements of x. The sigmoid
# function is also known as the logistic link function. It is the
# inverse of logit(x).
sigmoid <- function (x) {
  y   <- x
  y[] <- 0
  y[x > -500] <- 1/(1 + exp(-x))
  return(y)
}
# logpexp(x) returns log(1 + exp(x)). The computation is performed in a
# numerically stable manner. For large entries of x, log(1 + exp(x)) is
# effectively the same as x.
logpexp <- function (x) {
  y    <- x
  i    <- which(x < 16)
  y[i] <- log(1 + exp(x[i]))
  return(y)
}
# Use this instead of log(sigmoid(x)) to avoid loss of numerical precision.
logsigmoid <- function (x)
  -logpexp(-x)
# Compute the softmax of vector x in a more numerically prudent manner
# that avoids overflow or underflow.
softmax <- function (x) {
  y <- exp(x - max(x))
  return(y/sum(y))
}
###
###
###


n = 100
X = cbind(rnorm(n), rnorm(n, 1), rnorm(n, 1, 5))
b = c(2,5,-2)
Y = rbinom(n, 1, sigmoid(X %*% b + 2))



y_path = r"(../data/deletion.X.colnames_b30.simu_dele_30_0528.y)"

y = read_delim(y_path, "\n", col_names = "y", col_types = "i")
Y = y$y

X = read_delim("../data/block_439_442/block_439_442", "\t")

X = as.matrix(X - colMeans(X))

p0 = rep(1/length(b), length(b))


bayes.logistic(X, Y, p0, 0, 1)








X <- as.matrix(read.table(gzfile(${_input:r}), header = TRUE))
X <- scale(X, center = TRUE, scale = FALSE)
y <- as.matrix(read.table("${phenotype_file}"))
p <- dim(X)[2]
priors <- read.table("${hyperparam_file}")
p0 <- rep(priors[1,1], 1, p)
print ("Fitting Bayesian logistic regression model ...")
out <- bayes.logistic(X, y, p0, priors[2,1], priors[3,1])
print ("Model fitting completed!")
pip <- out$p1 * ${expected_effects}
names(pip) <- colnames(X)
write.table(t(t(pip)), ${_output:r}, sep='\t', col.names=FALSE, quote=FALSE)