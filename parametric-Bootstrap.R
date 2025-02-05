
# Data (v1 and v2)
x1 <- c(0.65, 0.48, 0.08, 0.11, 0.80, 0.67, 0.38, 0.63, 0.10, 0.03)
x2 <- c(0.10, 0.86, 0.64, 0.38, 0.57, 0.60, 0.50, 0.08, 0.10, 0.07)

# MLE for psi / function
psi.ml <- function(x) - length(x)/sum(log(1-x))

## MLE for v1 and v2  ::
mle_v1 <- psi.ml(x1) ; mle_v2 <- psi.ml(x2)
c(mle_v1 = round(mle_v1, 2))
c(mle_v2 = round(mle_v2, 2))


### mle for tau:::, v1 only in this case ....

tau.ml <- function(x) 1/(1 + psi.ml(x))
tau_v1 = tau.ml(x1) ; tau_v2 <- tau.ml(x2)
c(tau_v1 = round(tau_v1, 2))
c(tau_v2 = round(tau_v2, 2))

## SE for tau ::
se.tau <- function(x) (1/((1+psi.ml(x))^2)*(abs(psi.ml(x)))/(sqrt(length(x))))

#alpha(s)
alpha1 <- 0.1 ; alpha2 <- 0.01 

#normal quantiles ::
z1 <- qnorm(1- alpha1/2) ; z2 =qnorm(1 - alpha2/2)

# intervals :: 
round(c(tau_v1 -z1*se.tau(x1), tau_v1 + z1*se.tau(x1)), 3)
round(c(tau_v2 -z2*se.tau(x2), tau_v2 + z2*se.tau(x2), 3)

      
      
####### PARAMETRIC BOOTSTRAP !!! :::
 
n <- length(x1)   # sample size 
B <- 1000         # bootstrap size 

tau.boot <- rep(NA, B)
set.seed(123)

for (b in 1:B) {
  x.boot <- rbeta(n, 1, mle_v1)
  tau.boot[b] <- 1/(1 + psi.ml(x.boot))
}

cat('PBoot / BIAS :')
round(mean(tau.boot) - tau.ml(x1),3)

cat('PBoot / SE :')
round(sd(tau.boot), 3)

cat('Asymptotic - SE :')
round(se.tau(x1), 3)

cat( ' Pboot CI(0.90) --/n')
round(
  c(lower = tau.ml(x1) -qnorm(1-0.01/2)*sd(tau.boot),
    upper = tau.ml(x1) + qnorm(1-0.01/2)*sd(tau.boot)), 3)

    