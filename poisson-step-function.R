CalcPoisLogLike <- function(Y, X, cutoffs, heights) {
  if (length(heights) != length(cutoffs) + 1) {
    stop('Error: heights should be one element longer than cutoffs.')
  }
  mean_function <- rep(heights[1], length(X))
  for (ii in 2:length(heights)) {
    wh <- which(X >= cutoffs[ii - 1])
    mean_function[wh] <- heights[ii]
  }
  log_like <- dpois(Y, lambda = mean_function, log = TRUE)
  return(sum(log_like))
}


set.seed(1234)
N <- 1000
minX <- 0
maxX <- 10
X <- runif(N, minX, maxX)
true_means <- c(2, 20, 40)
true_cut <- c(3, 6)
mX <- true_means[1] * (X < true_cut[1]) +
  true_means[2] * (X >= true_cut[1] & X < true_cut[2]) +
  true_means[3] * (X >= true_cut[2])
Y <- rpois(N, lambda = mX)
par(mar = c(10, 1, 1, 1))
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


alpha_prior <- 0.01
beta_prior <- 0.01
Nsims <- 20000
heights <- NULL
heights[[1]] <- true_means
cutoffs <- NULL
cutoffs[[1]] <- true_cut
moves <- rep(NA, Nsims)

# moves correspond to 1 = 'H', 2 = 'P', 3 = 'B', 4 = 'D'.

for (ii in 2:Nsims) {
  curr_cut <- length(cutoffs[[ii - 1]])
  moves[ii] <- sample(1:4, 1, prob = c(etak[curr_cut + 1], pik[curr_cut + 1],
                                       bk[curr_cut + 1], dk[curr_cut + 1]))

  if (moves[ii] == 1) {  # Change in height.
    
    proposed_cutoffs <- cutoffs[[ii - 1]]
    proposed_heights <- heights[[ii - 1]]

    # Choose which height to change randomly.
    wh_height <- sample(1:length(proposed_heights), 1)
    proposed_heights[wh_height] <- proposed_heights[wh_height] *
                                     exp(runif(1, - 1/ 2, 1 / 2))
    logAR <- CalcPoisLogLike(Y, X, cutoffs = proposed_cutoffs,
                             heights = proposed_heights)
    logAR <- logAR - CalcPoisLogLike(Y, X, cutoffs = cutoffs[[ii - 1]],
                                     heights = heights[[ii - 1]])
    # Including the part of the prior.
    logAR <- logAR + alpha_prior * (log(proposed_heights)[wh_height] -
                                      log(heights[[ii - 1]])[wh_height])
    logAR <- logAR - beta_prior * ((proposed_heights - heights[[ii - 1]])[wh_height])
    
    cutoffs[[ii]] <- cutoffs[[ii - 1]]
    heights[[ii]] <- heights[[ii - 1]]
    if (log(runif(1)) < logAR) {
    	heights[[ii]] <- proposed_heights
    }
    
  } else if (moves[ii] == 2) {  # Change position of cutoff.
    
    proposed_cutoffs <- cutoffs[[ii - 1]]
    proposed_heights <- heights[[ii - 1]]
    
    # Which cutoff we will change.
    wh_cut <- sample(1:length(proposed_cutoffs), 1)
    cuts <- c(minX, proposed_cutoffs, maxX)
    choose_from <- c(cuts[wh_cut], cuts[wh_cut + 2])
    proposed_cutoffs[wh_cut] <- runif(1, min = choose_from[1], max = choose_from[2])
	      
    # Calculating the AR.
    logAR <- CalcPoisLogLike(Y = Y, X = X, cutoffs = proposed_cutoffs,
                             heights = proposed_heights)
    logAR <- logAR - CalcPoisLogLike(Y = Y, X = X, cutoffs = cutoffs[[ii - 1]],
                                     heights = heights[[ii - 1]])
    logAR <- logAR + log(choose_from[2] - proposed_cutoffs[wh_cut])
    logAR <- logAR + log(proposed_cutoffs[wh_cut] - choose_from[1])
    logAR <- logAR - log(choose_from[2] - cutoffs[[ii - 1]][wh_cut])
    logAR <- logAR + log(cutoffs[[ii - 1]][wh_cut] - choose_from[1])

    cutoffs[[ii]] <- cutoffs[[ii - 1]]
    heights[[ii]] <- heights[[ii - 1]]
    if (log(runif(1)) < logAR) {
      cutoffs[[ii]] <- proposed_cutoffs
    }
    
  } else if (moves[ii] == 3) {  # Birth to a new cutoff.
    
    # Proposing a new cutoff.
    sstar <- runif(1, minX, maxX)
    proposed_cutoffs <- sort(c(cutoffs[[ii - 1]], sstar))
    wh_cut <- which(proposed_cutoffs == sstar)

    sj <- ifelse(wh_cut == 1, minX, proposed_cutoffs[wh_cut - 1])
    sj1 <- ifelse(wh_cut == curr_cut + 1, maxX, proposed_cutoffs[wh_cut + 1])
    
    # Defining the proposed heights corresponding to the new cutoff.
    proposed_heights <- rep(NA, length(proposed_cutoffs) + 1)
    proposed_heights[- c(wh_cut, wh_cut + 1)] <- heights[[ii - 1]][- wh_cut]
    hj_prev <- heights[[ii - 1]][wh_cut]
    u <- runif(1, 0, 1)
    hj_new <- hj_prev * exp((sj1 - sstar) / (sj1 - sj) * log(u / (1 - u)))
    proposed_heights[wh_cut] <- hj_new
    hj1_new <- hj_new * (1 - u) / u
    proposed_heights[wh_cut + 1] <- hj1_new
    
    # Calculating the acceptance probability.
    # The likelihood ratio
    logAR <- CalcPoisLogLike(Y, X, proposed_cutoffs, proposed_heights)
    logAR <- logAR - CalcPoisLogLike(Y, X, cutoffs[[ii - 1]], heights[[ii - 1]])
    
    # The prior ratio
    k <- length(cutoffs[[ii - 1]])
    logAR <- logAR +
      dpois(k + 1, lambda, log = TRUE) - dpois(k, lambda, log = TRUE) +
      log(2 * (k + 1) * (2 * k + 3) / (maxX - minX) ^ 2) +
      log((sstar - sj) * (sj1 - sstar) / (sj1 - sj)) +
      (alpha_prior - 1) * log(hj_new * hj1_new / hj_prev) -
      beta_prior * (hj_new + hj1_new - hj_prev)
    
    # The proposal ratio
    logAR <- logAR + log(dk[k + 2] * (maxX - minX) / (bk[k + 1] * (k + 1)))
    
    # The jacobian
    logAR <- logAR + log((hj_new + hj1_new) ^ 2 / hj_prev)
    
    heights[[ii]] <- heights[[ii - 1]]
    cutoffs[[ii]] <- cutoffs[[ii - 1]]
    if (log(runif(1)) < logAR) {
      heights[[ii]] <- proposed_heights
      cutoffs[[ii]] <- proposed_cutoffs
    }
    
  } else {  # Death of a cutoff.
    
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
    hj_prev <- heights[[ii - 1]][wh_drop]
    hj1_prev <- heights[[ii - 1]][wh_drop + 1]
    
    h_new <- exp((sj1 - sj) / (sj2 - sj) * log(hj_prev) +
                 (sj2 - sj1) / (sj2 - sj) * log(hj1_prev))
    proposed_heights[wh_drop] <- h_new

    # Calculate the AR.
    logAR <- CalcPoisLogLike(Y, X, proposed_cutoffs, proposed_heights)
    logAR <- logAR - CalcPoisLogLike(Y, X, cutoffs[[ii - 1]], heights[[ii - 1]])
    
    # Add log prior ratio.
    k <- length(cutoffs[[ii - 1]])
    logAR <- logAR + 
      dpois(k - 1, lambda, log = TRUE) - dpois(k, lambda, log = TRUE) +
      log((maxX - minX) ^ 2 / (2 * k * (2 * k - 1))) +  # This might be wrong.
      log((sj2 - sj1) / ((sj1 - sj) * (sj2 - sj1))) +
      (alpha_prior - 1) * log(h_new ^ 2 / (hj_prev * hj1_prev)) -
      beta_prior * (h_new - hj_prev - hj1_prev)
    
    # Add log proposal ratio.
    logAR <- logAR + log(bk[k] * k / (dk[k + 1] * (maxX - minX)))
    
    # Add log Jacobian.
    logAR <- logAR + log(h_new / (hj_prev + hj1_prev) ^ 2)
    
    # Propose value and Accept/Reject.
    heights[[ii]] <- heights[[ii - 1]]
    cutoffs[[ii]] <- cutoffs[[ii - 1]]
    if (log(runif(1)) < logAR) {
      heights[[ii]] <- proposed_heights
      cutoffs[[ii]] <- proposed_cutoffs
    }
  }
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
lines(pred_x, apply(pred_y, 1, function(x) quantile(x, probs = 0.025)), col = 'green')
lines(pred_x, apply(pred_y, 1, function(x) quantile(x, probs = 0.975)), col = 'green')

pred_y <- pred_y[, - c(1:2000)]
plot(1, xlim = c(0, 10), type = 'n', ylim = range(pred_y), main = 'after burn in')
points(X, Y, pch = 16, cex = 0.3, col = 'red')
lines(pred_x, apply(pred_y, 1, mean))
lines(pred_x, apply(pred_y, 1, function(x) quantile(x, probs = 0.025)), col = 'green')
lines(pred_x, apply(pred_y, 1, function(x) quantile(x, probs = 0.975)), col = 'green')






