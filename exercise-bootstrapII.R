###########

# Data

xx <- c(0.55, 0.45, 0.10, 0.15, 0.80, 0.65, 0.38, 0.63, 0.10, 0.03)

### CHECK MLE for TAU !!!! 

psi_ml <- function(x) - length(x)/sum(log(1- x))
round(psi_ml(xx),2)


tau_ml <- function(x) 1/(1+psi_ml(x))
round(tau_ml(xx),2)

#### check c.i. for tau with wald test ::: with tau0 = 10 , and  aplha = 0.05 ::

tau0 = 10
se_tau <- function(x) (1/(1+psi_ml(x))^2)*(abs(psi_ml(x)))/(sqrt(length(x)))
abs(se_tau(xx)-tau0)/se_tau(xx)
qnorm(1-0.05/2)

### cause 134 > 1.95 (Za/2) we reject the null at aplha = 0.05


################################################################

## 1. Plot the likelihood function L(ψ | x10) highlithing the mle.

# Likelihood function
L <- function(psi, x){
  n <- length(x)
  s <- prod((1 - x))
  out <- (psi^n)*(s^(psi - 1))
  return(out)
}

# ---- simplest /shortest way ::
L1 <- function(psi, x) (psi^length(x)) * (prod(1 - x)^(psi - 1))

# MLE
mle = psi_ml(xx)

# Plot
curve(L(x, xx), from = 0, to = 4, lwd = 4, col = "orchid",
      xlab = expression(psi), ylab = "Likelihood function")
segments(mle, 0, mle, L(mle, xx), lty = 3, col = "green")
text(mle, 0.05, bquote(hat(psi) == .(round(mle,2))),
     pos = 2, cex = .8)
grid()


# Plot the likelihood function for L1 
curve(L1(x, xx), from = 0, to = 4, lwd = 2, col = "green",
      xlab = expression(psi), ylab = "Likelihood1")

################################################################



# Notice that the model we are dealing with is just a Beta(1, ψ) so use the function rbeta() to
# implement a parametric bootstrap to evaluate the bias and standard error of the mle for τ . Then
# build a 99% confidence interval again for the population mean τ .

### for the bootstrap we need a 10000 B sample as follow ::

## PARAMETRIC BOOTSTRAP ::

n <- length(xx)
B <- 1000
tau.boot <- rep(NA, B)
set.seed(123)
for (b in 1:B) {
  x.boot <- rbeta(n, 1, mle)
  tau.boot[b] <- 1/(1+ psi_ml(x.boot))
}
round(mean(tau.boot)- tau_ml(xx),3)
round(sd(tau.boot),3)

round(
  c(lower = tau_ml(xx) - qnorm(1-0.01/2)*sd(tau.boot),
    upper = tau_ml(xx) + qnorm(1-0.01/2)*sd(tau.boot)),3)



##### non parametric::


n <- length(xx)
B <- 1000
tau.boot <- rep(NA, B)
set.seed(123)
for (b in 1:B) {
  idx <- sample(1:n , replace = T)
  x.boot <- xx[idx]
  tau.nboot[b] <- 1/(1 + psi_ml(x.boot))
}
round(mean(tau.nboot) - tau_ml(xx),3)
round(sd(tau.nboot),3)

