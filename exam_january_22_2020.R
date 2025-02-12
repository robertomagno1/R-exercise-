# The data
xx <- c(9.47, 10.91, 41.35, 8.06, 15.15, 1.78, 21.61, 22.87, 9.96, 15.75)

# The MLE
beta_hat <- function(x) 1/mean(x)

# The (estimated) SE
se_hat <- function(x) abs(beta_hat(x))/sqrt(length(x))

# The Normal quantiles
alpha = 0.1
qq <- -qnorm(alpha/2)



qq

# Point estimate
round( beta_hat(xx), 3)
## [1] 0.064

# Asymptotic confidence interval (alternativa più chiara)
lower_bound <- beta_hat(xx) - qq * se_hat(xx)
upper_bound <- beta_hat(xx) + qq * se_hat(xx)

# Output finale in un formato più leggibile
round(c(lower_bound, upper_bound), 3)


# Asymptotic confidence interval
round( beta_hat(xx) + c(-1,+1)*qq*se_hat(xx), 3)



#### alternative version:.

# The data
xx <- c(9.47, 10.91, 41.35, 8.06, 15.15, 1.78, 21.61, 22.87, 9.96, 15.75)

# The MLE
beta_hat <- function(x) 1/mean(x)

# The (estimated) SE
se_hat <- function(x) abs(beta_hat(x))/sqrt(length(x))

# The Normal quantiles
alpha <- 0.1
qq <- qnorm(1 - alpha/2)  # Correzione del segno

# Point estimate
round(beta_hat(xx), 3)

# Asymptotic confidence interval (alternativa più chiara)
lower_bound <- beta_hat(xx) - qq * se_hat(xx)
upper_bound <- beta_hat(xx) + qq * se_hat(xx)

# Output finale in un formato più leggibile
round(c(lower_bound, upper_bound), 3)


require(MASS, quietly = TRUE)

??MASS

# Load the package
require(MASS, quietly = TRUE)


# For the Gamma model, the starting values for the quasi-Newton
# optimization method are not required.

gamma_fit <- fitdistr(xx, "gamma")
gamma_fit # maximum likelihood estimates with their (approximated) SE's


# Extract the estimates
shape_hat <- gamma_fit$estimate["shape"]
rate_hat <- gamma_fit$estimate["rate"]

# Extract MLE's full var/covar matrix
Sigma <- gamma_fit$vcov

# MLE for the population mean, by equivariace:
mean_hat <- shape_hat/rate_hat

# Function to evaluate the Jacobian of the transformation
grad <- function(shape, rate) c( 1/rate, -shape/(rate^2) ) # Jacobian of the transf.

# Standard error of <mean_hat> via delta-method
se_var <- sqrt( t(grad(shape_hat, rate_hat)) %*% Sigma %*% grad(shape_hat, rate_hat) )
round( se_var, 3 )


# Asymptotic confidence interval at 1% for the population mean (in minutes)
round( mean_hat + c(-1, +1) * qq * c(se_var), 3)

par(mfrow = c(1, 3))


### PDF
hist(xx, probability = TRUE, breaks = 5, col = "orchid", border = "white",
     xlab = "", ylab = "", main = "Densities")
curve( dgamma(x, shape = shape_hat, rate = rate_hat),
       lwd = 4, col = "purple", add = TRUE)
rug(xx)
grid()
legend("topright", c("Data", "Gamma fit"), lwd = 4, col = c("orchid", "purple"),
       bty = "n", cex = .8)


### CDF
curve( pgamma(x, shape = shape_hat, rate = rate_hat),
       from = min(xx), to = max(xx),
       lwd = 4, col = "purple", xlab = "", ylab = "",
       main = "Cumulative Distributions" )
plot( ecdf( xx ), col = "orchid", cex = .7, add = TRUE)
grid()
legend("topleft", c("Data", "Gamma fit"), lwd = 4, col = c("orchid", "purple"),
       bty = "n", cex = .8)


### QQ-Plot / Compare empirical vs theoretical quantiles at various levels
qqplot(qgamma(ppoints(length(xx)), shape = shape_hat, rate = rate_hat), xx,
       pch = 19, cex = .8, col = "orchid",
       main = "Gamma Q-Q Plot", xlab = "")
qqline(xx, distribution = function(p) qgamma(p, shape = shape_hat, rate = rate_hat),
       prob = c(0.1, 0.6), col = rgb(0,0,0,.3), lty = 3, lwd = 3)

