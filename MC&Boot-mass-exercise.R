require(MASS)
?galaxies
data('galaxies')

galaxies_data

T_n <- function(x, p = 0.25) {
  n <- length(x)
  k <- round(n * p)  # Compute k
  sorted_x <- sort(x)  # Sort values
  trimmed_x <- sorted_x[(k + 1):(n - k)]  # Remove extreme values
  mean(trimmed_x)  # Compute the trimmed mean
}

T_n(galaxies)

#### MonteCarlo simulation to check , simulate from t student :

set.seed(123)

M <- 1000  # Number of simulations
n <- 21  # Sample size
p <- 0.25  # Trimming proportion
true_mean <- 0  # Theoretical mean of t-distribution (df=2)

mse_Xbar <- 0
mse_Tn <- 0

for (i in 1:M) {
  sample_data <- rt(n, df = 2)  # Generate sample from t-distribution
  
  Xbar <- mean(sample_data)  # Compute standard mean
  Tn_k <- T_n(sample_data, p)  # Compute trimmed mean
  
  mse_Xbar <- mse_Xbar + (Xbar - true_mean)^2
  mse_Tn <- mse_Tn + (Tn_k - true_mean)^2
}

mse_Xbar <- mse_Xbar / M
mse_Tn <- mse_Tn / M

cat("MSE of Sample Mean:", mse_Xbar, "\n")
cat("MSE of Trimmed Mean:", mse_Tn, "\n")

# the trimmed mean is more robust to outliers.


# Part 2: Bootstrap Confidence Interval for the Galaxies Dataset

set.seed(42)  # Ensure reproducibility

B <- 1000  # Number of bootstrap resamples
bootstrap_Tn <- numeric(B)

for (i in 1:B) {
  resample <- sample(galaxies_data, replace = TRUE)  # Resample with replacement
  bootstrap_Tn[i] <- T_n(resample, p)  # Compute trimmed mean
}

# Compute 95% confidence interval (percentiles)
CI <- quantile(bootstrap_Tn, c(0.025, 0.975))

cat("Bootstrap 95% Confidence Interval for Mean:", CI, "\n")

