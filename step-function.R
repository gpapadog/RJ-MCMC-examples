CalcLogLike <- function(Y, X, cutoffs, heights, sigma) {
  if (length(heights) != length(cutoffs) + 1) {
    stop('Error: heights should be one element longer than cutoffs.')
  }
  mean_function <- rep(heights[1], length(X))
  for (ii in 2:length(heights)) {
    wh <- which(X >= cutoffs[ii - 1])
    mean_function[wh] <- heights[ii]
  }
  log_like <- dnorm(Y, mean = mean_function, sd = sigma, log = TRUE)
  return(sum(log_like))
}


set.seed(1234)
N <- 1000
sigma <- 1
minX <- 0
maxX <- 10
#X <- runif(N, minX, maxX)
X <- rnorm(N, mean = 5, sd = 1.5)
minX <- min(X)
maxX <- max(X)
true_means <- c(1, 2, 3)
true_cut <- c(3, 6)
mX <- true_means[1] * (X >= minX & X < true_cut[1]) +
  true_means[2] * (X >= true_cut[1] & X < true_cut[2]) +
  true_means[3] * (X >= true_cut[2])
Y <- rnorm(N, mean = mX, sd = sigma)
#Y <- rnorm(N, mean = 0, sd = sigma)
plot(X, mX)
plot(X, Y)

max_cutoffs <- 10
lambda <- 5  # Poisson parameter for number of breaks.

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
Nsims <- 10000
heights <- NULL
heights[[1]] <- c(4, 2, 0)
cutoffs <- NULL
cutoffs[[1]] <- c(5, 9)
sigmas <- rep(NA, Nsims)
sigmas[1] <- 1
moves <- rep(NA, Nsims)
range_unif <- 1

# moves correspond to 1 = 'H', 2 = 'P', 3 = 'B', 4 = 'D'.



for (ii in 2:Nsims) {
  curr_cut <- length(cutoffs[[ii - 1]])
  moves[ii] <- sample(1:4, 1, prob = c(etak[curr_cut + 1], pik[curr_cut + 1],
                                       bk[curr_cut + 1], dk[curr_cut + 1]))
  
  # MOVE 1: Change in height.
  if (moves[ii] == 1) {
    
    # For change in heigth we can do Gibbs sampling, so we don't propose.
    cutoffs[[ii]] <- cutoffs[[ii - 1]]
    heights[[ii]] <- heights[[ii - 1]]
    
    # Choose which height to change randomly.
    wh_height <- sample(1:length(heights[[ii]]), 1)
    
    # Which observations correspond to that height.
    if (wh_height == 1) {
      wh_obs <- which(X < ifelse(length(cutoffs[[ii]]) > 0,
                                 cutoffs[[ii]][1], maxX + 1))
    } else if (wh_height == curr_cut + 1) {  # Last height.
      wh_obs <- which(X >= cutoffs[[ii]][curr_cut])
    } else {
      wh_obs <- which(X >= cutoffs[[ii]][wh_height - 1] &
                        X <= cutoffs[[ii]][wh_height])
    }
    # Update using Gibbs.
    ssq_h <- 1 / (length(wh_obs) / sigmas[ii - 1] + 1 / uninf)
    mu_h <- ssq_h * sum(Y[wh_obs] / sigmas[ii - 1])
    heights[[ii]][wh_height] <- rnorm(1, mean = mu_h, sd = sqrt(ssq_h))
    
    
    # MOVE 2: Change in position of cutoff.
  } else if (moves[ii] == 2) {
    
    proposed_cutoffs <- cutoffs[[ii - 1]]
    proposed_heights <- heights[[ii - 1]]
    
    # Which cutoff we will change.
    wh_cut <- sample(1:length(proposed_cutoffs), 1)
    cuts <- c(min(X), proposed_cutoffs, max(X))
    choose_from <- c(cuts[wh_cut], cuts[wh_cut + 2])
    proposed_cutoffs[wh_cut] <- runif(1, min = choose_from[1], max = choose_from[2])
    
    # Calculating the AR.
    logAR <- CalcLogLike(Y = Y, X = X, cutoffs = proposed_cutoffs,
                         heights = proposed_heights, sigma = sigmas[ii - 1])
    logAR <- logAR - CalcLogLike(Y = Y, X = X, cutoffs = cutoffs[[ii - 1]],
                                 heights = heights[[ii - 1]], sigma = sigmas[ii - 1])
    
    # I think the additional term is the same as in Poisson, but I need to check.
    logAR <- logAR + log(choose_from[2] - proposed_cutoffs[wh_cut])
    logAR <- logAR + log(proposed_cutoffs[wh_cut] - choose_from[1])
    logAR <- logAR - log(choose_from[2] - cutoffs[[ii - 1]][wh_cut])
    logAR <- logAR + log(cutoffs[[ii - 1]][wh_cut] - choose_from[1])
    
    
    cutoffs[[ii]] <- cutoffs[[ii - 1]]
    heights[[ii]] <- heights[[ii - 1]]
    if (log(runif(1)) < logAR) {
      cutoffs[[ii]] <- proposed_cutoffs
    }
    
    # MOVE 3: Birth to a new cutoff.
  } else if (moves[ii] == 3) {
    
    # Choose the new cutoff.
    sstar <- runif(1, minX, maxX)
    proposed_cutoffs <- sort(c(cutoffs[[ii - 1]], sstar))
    wh_cut <- which(proposed_cutoffs == sstar)
    
    sj <- ifelse(wh_cut == 1, minX, proposed_cutoffs[wh_cut - 1])
    sj1 <- ifelse(wh_cut == curr_cut + 1, maxX, proposed_cutoffs[wh_cut + 1])
    
    # Defining the proposed heights corresponding to the new cutoff.
    proposed_heights <- rep(NA, length(proposed_cutoffs) + 1)
    proposed_heights[- c(wh_cut, wh_cut + 1)] <- heights[[ii - 1]][- wh_cut]
    hj_prev <- heights[[ii - 1]][wh_cut]
    u <- runif(1, 0, range_unif)
    hj_new <- hj_prev - u * (sj1 - sstar) / (sj1 - sj)
    hj1_new <- hj_new + u * (sstar - sj) / (sj1 - sj)
    proposed_heights[c(wh_cut, wh_cut + 1)] <- c(hj_new, hj1_new)
    # NOTE: What if the change in heights is NOT uniform?
    
    # Calculating the acceptance probability.
    ## Likelihood ratio:
    logAR <- CalcLogLike(Y = Y, X = X, cutoffs = proposed_cutoffs,
                         heights = proposed_heights, sigma = sigmas[ii - 1])
    logAR <- logAR - CalcLogLike(Y = Y, X = X, cutoffs = cutoffs[[ii - 1]],
                                 heights = heights[[ii - 1]], sigma = sigmas[ii - 1])
    ## Prior ratio:
    # For the number of cutoffs:
    k <- length(cutoffs[[ii - 1]])
    logAR <- logAR + dpois(k + 1, lambda, log = TRUE)
    logAR <- logAR - dpois(k, lambda, log = TRUE)
    # For the cutoffs:
    logAR <- logAR + log(2 * (k + 1) * (2 * k + 3)) - log((maxX - minX) ^ 2)
    logAR <- logAR + log((sstar - sj) * (sj1 - sstar) / (sj1 - sj))
    # For the heights:
    logAR <- logAR + sum(dnorm(proposed_heights[c(wh_cut, wh_cut + 1)],
                               mean = 0, sd = sqrt(uninf), log = TRUE))
    logAR <- logAR - dnorm(heights[[ii - 1]][wh_cut], mean = 0, sd = sqrt(uninf),
                           log = TRUE)
    
    ## Proposal ratio:
    logAR <- logAR + log(dk[k + 2] * (maxX - minX) / (bk[k + 1] * (k + 1)))
    ## Jacobian is equal to 1.
    
    cutoffs[[ii]] <- cutoffs[[ii - 1]]
    heights[[ii]] <- heights[[ii - 1]]
    if (log(runif(1)) < logAR) {
      cutoffs[[ii]] <- proposed_cutoffs
      heights[[ii]] <- proposed_heights
    }
    
    # MOVE 4: DEATH OF A CUTOFF.
  } else {
    
    # Choose the cutoff to drop.
    wh_drop <- sample(1:curr_cut, 1)
    drop_cut <- cutoffs[[ii - 1]][wh_drop]
    
    # Proposing the new height.
    proposed_cutoffs <- setdiff(cutoffs[[ii - 1]], drop_cut)
    proposed_heights <- heights[[ii - 1]][- wh_drop]
    
    # Dropping sj1
    sj <- ifelse(wh_drop == 1, minX, proposed_cutoffs[wh_drop - 1])
    sj1 <- drop_cut
    sj2 <- ifelse(wh_drop == curr_cut, maxX, proposed_cutoffs[wh_drop])
    
    # Defining the new height as a weighted average of the two.
    hj_prev <- heights[[ii - 1]][wh_drop]
    hj1_prev <- heights[[ii - 1]][wh_drop + 1]
    h_new <- ((sj1 - sj) * hj_prev + (sj2 - sj1) * hj1_prev) / (sj2 - sj)
    proposed_heights[wh_drop] <- h_new
    
    # Calculate the AR.
    logAR <- CalcLogLike(Y, X, cutoffs = proposed_cutoffs, heights = proposed_heights,
                         sigma = sigmas[ii - 1])
    logAR <- logAR - CalcLogLike(Y, X, cutoffs = cutoffs[[ii - 1]],
                                 heights = heights[[ii - 1]], sigma = sigmas[ii - 1])
    # Add log prior ratio.
    # For the number of cutoffs:
    k <- length(cutoffs[[ii - 1]])
    logAR <- logAR + dpois(k - 1, lambda, log = TRUE) - dpois(k, lambda, log = TRUE)
    # For the cutoffs:
    logAR <- logAR + log((maxX - minX) ^ 2 / (2 * k * (2 * k - 1)))
    logAR <- logAR + log((sj2 - sj1) / ((sj1 - sj) * (sj2 - sj1)))
    # For the heights:
    logAR <- logAR +
      dnorm(h_new, mean = 0, sd = sqrt(uninf), log = TRUE) -
      sum(dnorm(c(hj_prev, hj1_prev), mean = 0, sd = sqrt(uninf), log = TRUE))
    
    # Add log proposal ratio.
    logAR <- logAR + log(bk[k] * k / (dk[k + 1] * (maxX - minX)))
    # Jacobian is 1.
    
    # Propose value and Accept/Reject.
    if (log(runif(1)) < logAR) {
      cutoffs[[ii]] <- proposed_cutoffs
      heights[[ii]] <- proposed_heights
    } else {
      cutoffs[[ii]] <- cutoffs[[ii - 1]]
      heights[[ii]] <- heights[[ii - 1]]
    }
  }
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
    pred_y[wh_obs, ii] <- heights[[ii]][cc]
  }
}
plot(1, xlim = c(0, 10), type = 'n', ylim = range(pred_y))
points(X, Y, pch = 16, cex = 0.3, col = 'red')
lines(pred_x, apply(pred_y, 1, mean))
lines(pred_x, apply(pred_y, 1, function(x) quantile(x, probs = 0.1)), col = 'green')
lines(pred_x, apply(pred_y, 1, function(x) quantile(x, probs = 0.9)), col = 'green')

pred_y <- pred_y[, - c(1:2000)]
plot(1, xlim = c(0, 10), type = 'n', ylim = range(pred_y), main = 'after burn in')
points(X, Y, pch = 16, cex = 0.3, col = 'red')
lines(pred_x, apply(pred_y, 1, mean))
lines(pred_x, apply(pred_y, 1, function(x) quantile(x, probs = 0.025)), col = 'green')
lines(pred_x, apply(pred_y, 1, function(x) quantile(x, probs = 0.975)), col = 'green')







