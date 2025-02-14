#####fixed ::: 

x1 <- c(0.10, 0.86, 0.64, 0.38, 0.57, 0.60, 0.50, 0.08, 0.10, 0.07)

# MLE Estimator function
psi.ml <- function(x) -length(x) / sum(log(1 - x))
round(psi.ml(x1), 2)

# MLE for τ
mle_tau <- function(x) 1 / (1 + psi.ml(x))
tau_v1 <- mle_tau(x1)
round(tau_v1, 2)

# SE for τ (Professor's version)
se.tau <- function(x) (1 / ((1 + psi.ml(x))^2)) * (abs(psi.ml(x))) / sqrt(length(x))
round(se.tau(x1), 2)

# Confidence Interval setup
z1 <- qnorm(1 - 0.01/ 2)

# Computing Confidence Interval
round(c(tau_v1 - z1 * se.tau(x1), tau_v1 + z1 * se.tau(x1)), 3)


###execise 1 : 2020/january/22::: we need to find the value of the mle ::

x1 <- c(10.2, 10.0, 8.0, 1.8, 6.5, 1.1, 9.3, 7.6, 10.1, 7.8)

psi_mle <- function(x) sqrt(sum(x^2)/(2*length(x)))
round(psi_mle(x1), 3)

tau_mle <- function(x) (2-pi/2)*psi_mle(x)^2
round(tau_mle(x1),3)



#### implement a non-parametric bootstrap ::


n <- length(x1)
B <- 1000
tau_boot = rep(NA, B)
set.seed(123)
for (b in 1:B) {
  x.boot <- rbeta(n,1, psi_mle(x1))
  tau_boot[b] <- tau_mle(psi_mle(x.boot))
}

round(mean(tau_boot) - tau_mle(psi_mle(x1)), 3)
round(sd(tau_boot), 3)




### corrected and short version ::;  still same result as before :...


n <- length(x1)  # Sample size
B <- 1000        # Bootstrap size
tau.pboot <- rep(NA, B)  # Init
set.seed(123)

for (b in 1:B) {
  x.boot <- rbeta(n, 1, psi_mle(x1))  # Generate synthetic data
  tau.pboot[b] <- tau_mle(psi_mle(x.boot))
}

# Bias & SE via Parametric Bootstrap
cat('--- PBoot / Bias ----\n', round(mean(tau.pboot) - tau_mle(psi_mle(x1)), 3))
cat('\n--- PBoot / SE ------\n', round(sd(tau.pboot), 3))

tau_ml1 <- tau_mle(x1)


### exercise find wald test for tau :: 
tau0 = 15

se_tau <- function(x) (abs(2*(2-pi/2)*psi_mle(x)))*abs(psi_mle(x))/(2*(sqrt(length(x))))

se_tau_modified <- function(x) abs(2 * (2 - pi / 2) * psi_mle(x) * psi_mle(x)) / (2 * sqrt(length(x)))


round(abs((se_tau(x1) - tau0)/se_tau(psi_mle(x1))),2)

#prof's version :: 
abs((tau_ml1 - tau0)/se_tau(x1))

abs((tau_ml1 - tau0)/se_tau_modified(x1))

se.tau1 = function(x) ( abs(2*(2 - pi/2) * psi_mle(x)) ) * (abs(psi_mle(x))) / (2 * sqrt(length(x)))

abs((tau_ml1 - tau0)/se.tau1(x1))


##since 0.39 < zaplha/2 ::: 1.69 for aplha = 0.05 we accept H0 so reject null Hp.....






