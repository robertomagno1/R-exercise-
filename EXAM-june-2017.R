
### this exercise required to find the mle for tau ....

# Data
xx <- c(0.65, 0.48, 0.08, 0.11, 0.80, 0.67, 0.38, 0.63, 0.10, 0.03)
# MLE for psi

ml_psi <- function(x) -length(x)/sum(log(1- x))
round(ml_psi(xx), 2)

ml_tau <- function(x) 1/(1+ml_psi(x))
round(ml_tau(xx), 3)







tau0 <- 10


se_tau <- function(x) (1/((1+ml_psi(x))^2) )*(abs(ml_psi(x)))/(sqrt(length(x)))
abs((ml_tau(xx)-tau0)/se_tau(xx))







### <C.I. for tau :::

tau0 = 10
# SE for tau
se.tau = function(x) ( 1/((1 + psi.ml(x))^2) )*(abs(psi.ml(x)))/(sqrt(length(x)))
# Wald statistics
abs((tau.ml(xx) - tau0)/se.tau(xx))

## cause 8.823 is bigger than t at alpha = 0.05 --->> 1.96 

### see this 
round(qnorm(1-0.05/2), 2)

########################################################


# Likelihood function
L <- function(psi, x){
  n <- length(x)
  s <- prod((1 - x))
  out <- (psi^n)*(s^(psi - 1))
  return(out)
}
# MLE
mle = psi.ml(xx)
# Plot
curve(L(x, xx), from = 0, to = 4, lwd = 4, col = "orchid",
      xlab = expression(psi), ylab = "Likelihood function")
segments(mle, 0, mle, L(mle, xx), lty = 3, col = "green")
text(mle, 0.05, bquote(hat(psi) == .(round(mle,2))),
     pos = 2, cex = .8)
grid()




# Data
xx <- c(0.55, 0.45, 0.10, 0.15, 0.80, 0.65, 0.38, 0.63, 0.10, 0.03)

mle_psi <- function(x) - length(x)/sum(log(1-x))
round(mle_psi(xx),2)

mle_tau <- function(x) 1/(1+mle_psi(x))
round(mle_tau(xx),2)

tau0 = 10

# SE for tau
se.tau = function(x) ( 1/((1 + psi.ml(x))^2) )*(abs(psi.ml(x)))/(sqrt(length(x)))

# Wald statistics
abs((tau.ml(xx) - tau0)/se.tau(xx))

## since this value is way larger than 1.96 so we definitely reject the null at α = 0.05.
qnorm(1-0.05/2)

