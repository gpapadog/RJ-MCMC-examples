# Author: Georgia Papadogeorgou
# Date: 9.27.2016
# Description: Doing RJ MCMC where I'm using the marginal likelihood of the data
# to update the cutoffs, and then I'm sampling the heights and effect of X as linear
# term.

library(MASS)

set.seed(1234)
N <- 1000
sigma <- 1
minX <- 0
maxX <- 10
X <- rnorm(N, mean = 5, sd = 1.5)
minX <- min(X)
maxX <- max(X)
true_means <- c(1, 2, 3)
true_cut <- c(3, 6)
# mX <- true_means[1] * (X >= minX & X < true_cut[1]) +
#   true_means[2] * (X >= true_cut[1] & X < true_cut[2]) +
#   true_means[3] * (X >= true_cut[2])
mX <- X ^ 2 / 10
Y <- rnorm(N, mean = mX, sd = sigma)
plot(X, mX)
plot(X, Y)

max_cutoffs <- 200
lambda <- 30  # Poisson parameter for number of breaks.

# Calculating the c parameter of bk, dk.
bk <- rep(0, max_cutoffs + 1)
names(bk) <- paste0('s=', c(0:max_cutoffs))
dk <- bk
etak <- bk
pik <- bk

for (ii in 0:(max_cutoffs - 1)) {  # Number of current cutoffs.
  bk[ii + 1] <- dpois(ii + 1, lambda = lambda) / dpois(ii, lambda = lambda)
  dk[ii + 2] <- dpois(ii, lambda = lambda) / dpois(ii + 1, lambda = lambda)
}
dk[length(dk)] <- dpois(max_cutoffs - 1, lambda = lambda) /
  dpois(max_cutoffs, lambda = lambda)
bk <- sapply(bk, function(x) min(x, 1))
dk <- sapply(dk, function(x) min(x, 1))

c <- 0.9 / max(bk + dk)
bk <- c * bk
dk <- c * dk
for (ii in 2:(max_cutoffs + 1)) {
  etak[ii] <- (1 - bk[ii] - dk[ii]) / 2
  pik[ii] <- etak[ii]
}
etak[1] <- 1 - bk[1]



uninf <- 100 ^ 2
Nsims <- 3000
heights <- NULL
heights[[1]] <- rbind(c(4, 2, 2, 0), c(0, 0, 0, 0))
cutoffs <- NULL
cutoffs[[1]] <- c(2, 5, 7)
sigmas <- rep(NA, Nsims)
sigmas[1] <- 1
moves <- rep(NA, Nsims)
range_unif <- 1
nu <- 0.001
lambda_ssq <- 0.001

# moves correspond to 1 = 'H', 2 = 'P', 3 = 'B', 4 = 'D'.

D <- data.frame(X = X, Y = Y)
mu_prior <- mean(Y)

CalcLogMarginalLike <- function(nu, lambda, X, V, Y, mu) {
  if (class(X) == 'numeric') {
    X <- matrix(X, nrow = length(X), ncol = 1)
  }
  N <- length(Y)
  r <- lgamma((nu + N) / 2) + (nu / 2) * log(nu * lambda)
  r <- r - N / 2 * log(pi) - lgamma(nu / 2)
  mat <- diag(N) + X %*% V %*% t(X)
  r <- r - 1 / 2 * log(det(mat))
  r <- r - (nu + N) / 2 * log(lambda * nu + t(Y - X %*% mu) %*% solve(mat) %*%
                                (Y - X %*% mu))
  return(r)
}


for (ii in 2:Nsims) {
  if (ii %% 100 == 0) {
    print(ii)
  }
  curr_cut <- length(cutoffs[[ii - 1]])
  moves[ii] <- sample(2:4, 1, prob = c(pik[curr_cut + 1], bk[curr_cut + 1],
                                       dk[curr_cut + 1]))
  
  # MOVE 2: Change to a cutoff.
  if (moves[ii] == 2) {
    
    proposed_cutoffs <- cutoffs[[ii - 1]]
    
    # Which cutoff we will change.
    wh_cut <- sample(1:length(proposed_cutoffs), 1)
    cuts <- c(min(X), proposed_cutoffs, max(X))
    choose_from <- c(cuts[wh_cut], cuts[wh_cut + 2])
    proposed_cutoffs[wh_cut] <- runif(1, min = choose_from[1], max = choose_from[2])
    
    # Calculating the AR.
    sj <- cuts[wh_cut]
    sj1 <- cuts[wh_cut + 1]
    sj2 <- cuts[wh_cut + 2]
    sj1star <- proposed_cutoffs[wh_cut]
    wh1 <- which(X >= sj & X <= sj1)
    wh2 <- which(X > sj1 & X <= sj2)
    wh1star <- which(X >= sj & X <= sj1star)
    wh2star <- which(X > sj1star & X <= sj2)
    
    if (length(wh1star) > 0 & length(wh2star) > 0) {
      logAR <- CalcLogMarginalLike(nu = nu, lambda = lambda_ssq,
                                   X = cbind(rep(1, length(wh1star)), X[wh1star]),
                                   V = diag(rep(uninf, 2)), Y = Y[wh1star], mu = c(mean(Y), 0)) + 
        CalcLogMarginalLike(nu = nu, lambda = lambda_ssq,
                            X = cbind(rep(1, length(wh2star)), X[wh2star]),
                            V = diag(rep(uninf, 2)), Y = Y[wh2star], mu = c(mean(Y), 0)) -
        CalcLogMarginalLike(nu = nu, lambda = lambda_ssq,
                            X = cbind(rep(1, length(wh1)), X[wh1]),
                            V = diag(rep(uninf, 2)), Y = Y[wh1], mu = c(mean(Y), 0)) -
        CalcLogMarginalLike(nu = nu, lambda = lambda_ssq,
                            X = cbind(rep(1, length(wh2)), X[wh2]),
                            V = diag(rep(uninf, 2)), Y = Y[wh2], mu = c(mean(Y), 0))
      # I think the additional term is the same as in Poisson, but I need to check.
      logAR <- logAR + log(choose_from[2] - proposed_cutoffs[wh_cut])
      logAR <- logAR + log(proposed_cutoffs[wh_cut] - choose_from[1])
      logAR <- logAR - log(choose_from[2] - cutoffs[[ii - 1]][wh_cut])
      logAR <- logAR + log(cutoffs[[ii - 1]][wh_cut] - choose_from[1])
      
      cutoffs[[ii]] <- cutoffs[[ii - 1]]
      if (log(runif(1)) < logAR) {
        cutoffs[[ii]] <- proposed_cutoffs
      }
    } else {
      cutoffs[[ii]] <- cutoffs[[ii - 1]]
    }
    
    # MOVE 3: Birth to a new cutoff.
  } else if (moves[ii] == 3) {
    
    # Choose the new cutoff.
    sstar <- runif(1, minX, maxX)
    proposed_cutoffs <- sort(c(cutoffs[[ii - 1]], sstar))
    wh_cut <- which(proposed_cutoffs == sstar)
    
    sj <- ifelse(wh_cut == 1, minX, proposed_cutoffs[wh_cut - 1])
    sj1 <- ifelse(wh_cut == curr_cut + 1, maxX, proposed_cutoffs[wh_cut + 1])
    
    # Calculating the acceptance probability.
    ## Likelihood ratio:
    wh1 <- which(X >= sj & X <= sstar)
    wh2 <- which(X > sstar & X <= sj1 )
    
    if (length(wh1) > 0 & length(wh2) > 0) {
      logAR <- CalcLogMarginalLike(nu = nu, lambda = lambda_ssq,
                                   X = cbind(rep(1, length(wh1)), X[wh1]),
                                   V = diag(rep(uninf, 2)), Y = Y[wh1], mu = c(mean(Y), 0)) + 
        CalcLogMarginalLike(nu = nu, lambda = lambda_ssq,
                            X = cbind(rep(1, length(wh2)), X[wh2]),
                            V = diag(rep(uninf, 2)), Y = Y[wh2], mu = c(mean(Y), 0)) -
        CalcLogMarginalLike(nu = nu, lambda = lambda_ssq,
                            X = cbind(rep(1, length(wh1) + length(wh2)),
                                      X[c(wh1, wh2)]),
                            V = diag(rep(uninf, 2)), Y = Y[c(wh1, wh2)], mu = c(mean(Y), 0))
      
      ## Prior ratio:
      # For the number of cutoffs:
      k <- length(cutoffs[[ii - 1]])
      logAR <- logAR + dpois(k + 1, lambda, log = TRUE)
      logAR <- logAR - dpois(k, lambda, log = TRUE)
      # For the cutoffs:
      logAR <- logAR + log(2 * (k + 1) * (2 * k + 3)) - log((maxX - minX) ^ 2)
      logAR <- logAR + log((sstar - sj) * (sj1 - sstar) / (sj1 - sj))
      
      ## Proposal ratio:
      logAR <- logAR + log(dk[k + 2] * (maxX - minX) / (bk[k + 1] * (k + 1)))
      ## Jacobian is equal to 1.
      
      cutoffs[[ii]] <- cutoffs[[ii - 1]]
      if (log(runif(1)) < logAR) {
        cutoffs[[ii]] <- proposed_cutoffs
      }
    } else {
      cutoffs[[ii]] <- cutoffs[[ii - 1]]
    }
    
    # MOVE 4: DEATH OF A CUTOFF.
  } else {
    
    # Choose the cutoff to drop.
    wh_drop <- sample(1:curr_cut, 1)
    drop_cut <- cutoffs[[ii - 1]][wh_drop]
    proposed_cutoffs <- setdiff(cutoffs[[ii - 1]], drop_cut)
    
    # Dropping sj1
    sj <- ifelse(wh_drop == 1, minX, proposed_cutoffs[wh_drop - 1])
    sj1 <- drop_cut
    sj2 <- ifelse(wh_drop == curr_cut, maxX, proposed_cutoffs[wh_drop])
    
    # Calculate the AR.
    wh1 <- which(X > sj & X <= sj1)
    wh2 <- which(X > sj1 & X <= sj2)
    logAR <- CalcLogMarginalLike(nu = nu, lambda = lambda_ssq,
                                 X = cbind(rep(1, length(wh1) + length(wh2)),
                                           X[c(wh1, wh2)]),
                                 V = diag(rep(uninf, 2)), Y = Y[c(wh1, wh2)], mu = c(mean(Y), 0)) -
      CalcLogMarginalLike(nu = nu, lambda = lambda_ssq,
                          X = cbind(rep(1, length(wh1)), X[wh1]),
                          V = diag(rep(uninf, 2)), Y = Y[wh1], mu = c(mean(Y), 0)) -
      CalcLogMarginalLike(nu = nu, lambda = lambda_ssq,
                          X = cbind(rep(1, length(wh2)), X[wh2]),
                          V = diag(rep(uninf, 2)), Y = Y[wh2], mu = c(mean(Y), 0))
    
    # Add log prior ratio.
    # For the number of cutoffs:
    k <- length(cutoffs[[ii - 1]])
    logAR <- logAR + dpois(k - 1, lambda, log = TRUE) - dpois(k, lambda, log = TRUE)
    # For the cutoffs:
    logAR <- logAR + log((maxX - minX) ^ 2 / (2 * k * (2 * k - 1)))
    logAR <- logAR + log((sj2 - sj1) / ((sj1 - sj) * (sj2 - sj1)))
    
    # Add log proposal ratio.
    logAR <- logAR + log(bk[k] * k / (dk[k + 1] * (maxX - minX)))
    # Jacobian is 1.
    
    # Propose value and Accept/Reject.
    if (log(runif(1)) < logAR) {
      cutoffs[[ii]] <- proposed_cutoffs
    } else {
      cutoffs[[ii]] <- cutoffs[[ii - 1]]
    }
  }
}


for (ii in 2:Nsims) {
  cuts <- c(minX, cutoffs[[ii]], maxX)
  h <- matrix(nrow = 2, ncol = length(cuts) - 1)
  for (cc in 1:(length(cuts) - 1)) {
    wh_obs <- which(D$X > cuts[cc] & D$X <= cuts[cc + 1])
    design_mat <- cbind(rep(1, length(wh_obs)), X[wh_obs])
    Sigma_beta <- solve(t(design_mat) %*% design_mat / sigmas[ii - 1] +
                          diag(2) / uninf)
    mu_beta <- Sigma_beta %*% (t(design_mat) %*% Y[wh_obs] / sigmas[ii - 1] +
                                 c(mean(Y), 0) / uninf)
    h[, cc] <- mvrnorm(1, mu = mu_beta, Sigma = Sigma_beta)
  }
  heights[[ii]] <- h
  sigmas[ii] <- sigmas[ii - 1]
}



pred_x <- seq(minX, maxX, length.out = 100)
pred_x <- pred_x[- c(1, length(pred_x))]
pred_y <- matrix(NA, nrow = length(pred_x), ncol = length(heights))
for (ii in 1:length(heights)) {
  cuts <- rep(cutoffs[[ii]], each = 2)
  cuts <- c(minX, cuts, maxX)
  for (cc in 1:(length(cuts) / 2)) {
    wh_obs <- which(pred_x > cuts[2 * cc - 1] & pred_x <= cuts[2 * cc])
    pred_y[wh_obs, ii] <- heights[[ii]][1, cc] + heights[[ii]][2, cc] * pred_x[wh_obs]
  }
}
plot(1, xlim = c(0, 10), type = 'n', ylim = range(pred_y))
points(X, Y, pch = 16, cex = 0.3, col = 'red')
lines(pred_x, apply(pred_y, 1, mean))
lines(pred_x, apply(pred_y, 1, function(x) quantile(x, probs = 0.1)), col = 'green')
lines(pred_x, apply(pred_y, 1, function(x) quantile(x, probs = 0.9)), col = 'green')

pred_y <- pred_y[, - c(1:1000)]
plot(1, xlim = c(0, 10), type = 'n', ylim = range(pred_y), main = 'after burn in')
points(X, Y, pch = 16, cex = 0.3, col = 'red')
lines(pred_x, apply(pred_y, 1, mean))
lines(pred_x, apply(pred_y, 1, function(x) quantile(x, probs = 0.025)), col = 'green')
lines(pred_x, apply(pred_y, 1, function(x) quantile(x, probs = 0.975)), col = 'green')

table(sapply(cutoffs, length))






