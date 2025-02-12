##Consider μb0, that is the estimator obtained by setting k = 0, and the sample mean X -n 


#In R, plot their mse’s as functions of the parameter μ: does X ̄n  uniformly dominate μb0?

# Definiamo le funzioni di rischio (MSE)
risk_mu0 <- function(mu, n) {
  1/((n+1)^2) * (mu^2) + n/((n+1)^2)
}

risk_xbar <- function(mu, n) {
  rep(1/n, length(mu))  # Restituisce un valore costante per tutti i mu
}

# Valori di n che vogliamo testare
n_seq <- c(10, 50, 100, 1000)

# Impostiamo un layout per più grafici
par(mfrow = c(2,2), mar = c(4,3,3,1))

# Creiamo i grafici per ciascun valore di n
for (nn in n_seq) {
  curve(risk_mu0(x, nn), -2, 2, 
        main = paste("n =", nn), 
        xlab = expression(mu), ylab = "MSE", 
        col = "purple", lwd = 5)
  
  curve(risk_xbar(x, nn), col = "pink", lwd = 5, add = TRUE)
  
  legend("top", c("new estimator", "sample average"),
         col = c("purple", "pink"), lwd = 5, bty = "n", cex = 0.8)
}
