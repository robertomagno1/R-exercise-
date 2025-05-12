# ───────────────────────────────────────────────────────────────────────────────
# Mood Prediction via Custom Kernel Ridge + Conformal + Marginal Effects
# ───────────────────────────────────────────────────────────────────────────────

# 1) LOAD LIBRARIES -------------------------------------------------------------
pkgs <- c("caret","dplyr","ggplot2")
for(p in pkgs) if(!p %in% installed.packages()[,"Package"]) install.packages(p)
lapply(pkgs, library, character.only=TRUE)

# 2) LOAD & PREPROCESS ----------------------------------------------------------
data_path <- "/Users/roberto/Desktop/UNI_DATASCIENCE/Stat_learning/final_project/risposte_umore.csv"
df        <- read.csv(data_path, stringsAsFactors=FALSE)

# Clean up column names
colnames(df) <- c("Timestamp","Energy","SleepQuality","Hunger",
                  "PhysicalActivity","StatsLiking","Informed",
                  "BogosBinted","Anger","Fear","Happiness",
                  "Sadness","Surprise","Disgust")

# Turn PhysicalActivity into factor No/Low/High
df <- df %>%
  mutate(
    PhysAct = case_when(
      grepl("No",      PhysicalActivity, ignore.case=TRUE) ~ "No",
      grepl("leggera", PhysicalActivity, ignore.case=TRUE) ~ "Low",
      grepl("intensa", PhysicalActivity, ignore.case=TRUE) ~ "High",
      TRUE                                                  ~ "No"
    ),
    PhysAct = factor(PhysAct, levels=c("No","Low","High"))
  ) %>%
  select(-PhysicalActivity)

# Features & target
features <- df %>% select(Energy, SleepQuality, Hunger,
                          StatsLiking, Informed, BogosBinted, PhysAct)
y_all    <- df$Happiness

# 3) SPLIT 60/20/20 --------------------------------------------------------------
set.seed(2025)
n     <- nrow(df)
i_tr  <- sample(n, size = 0.6*n)
rest  <- setdiff(seq_len(n), i_tr)
i_ca  <- sample(rest, size = 0.5*length(rest))
i_te  <- setdiff(rest, i_ca)

X_tr_raw <- features[i_tr,]; y_tr  <- y_all[i_tr]
X_ca_raw <- features[i_ca,]; y_ca  <- y_all[i_ca]
X_te_raw <- features[i_te,]; y_te  <- y_all[i_te]

# 4) DUMMY‐ENCODE & SCALE -------------------------------------------------------
dumm   <- caret::dummyVars(~ ., data=X_tr_raw)
X_tr_n <- predict(dumm, X_tr_raw)
X_ca_n <- predict(dumm, X_ca_raw)
X_te_n <- predict(dumm, X_te_raw)

pp     <- caret::preProcess(X_tr_n, method=c("center","scale"))
X_tr   <- predict(pp, X_tr_n)
X_ca   <- predict(pp, X_ca_n)
X_te   <- predict(pp, X_te_n)

# Convert to numeric matrices
X_tr_mat <- as.matrix(X_tr)
X_ca_mat <- as.matrix(X_ca)
X_te_mat <- as.matrix(X_te)

# 5) DEFINE RBF KERNEL & GRID ---------------------------------------------------
rbf_kernel <- function(X, Y, sigma) {
  # returns nrow(X)×nrow(Y) kernel matrix
  X2 <- rowSums(X^2)
  Y2 <- rowSums(Y^2)
  # squared distances
  D2 <- outer(X2, Y2, "+") - 2 * (X %*% t(Y))
  exp(-D2 / (2 * sigma^2))
}

# Heuristic for sigma: median pairwise distance on train
paird <- as.vector(dist(X_tr_mat))
sigma0 <- median(paird)

sigmas  <- sigma0 * c(0.5, 1, 2)
lambdas <- 10^seq(-3, 1, length.out=5)

best <- list(mse=Inf)

# 6) GRID SEARCH ON CALIBRATION -------------------------------------------------
for(s in sigmas) {
  K_tr <- rbf_kernel(X_tr_mat, X_tr_mat, s)
  for(lambda in lambdas) {
    # solve α = (K + λ I)^{-1} y_tr
    A     <- K_tr + lambda * diag(nrow(K_tr))
    alpha <- solve(A, y_tr)
    # predict on calibration
    K_ca  <- rbf_kernel(X_ca_mat, X_tr_mat, s)  # n_ca × n_tr
    pred_ca <- K_ca %*% alpha
    mse_ca  <- mean((pred_ca - y_ca)^2)
    if(mse_ca < best$mse) {
      best <- list(sigma=s, lambda=lambda, alpha=alpha)
    }
  }
}

# unpack best
sigma_hat <- best$sigma
lambda_hat<- best$lambda
alpha_hat <- best$alpha

# 7) CONFORMAL CALIBRATION ------------------------------------------------------
alpha   <- 0.1
# recompute cal predictions
K_ca    <- rbf_kernel(X_ca_mat, X_tr_mat, sigma_hat)
pred_ca <- K_ca %*% alpha_hat
nc      <- abs(y_ca - pred_ca)
qlev    <- (1 - alpha) * (1 + 1/length(nc))
q_krr   <- quantile(nc, probs = qlev, type = 1)

# 8) PREDICTION + INTERVAL FUNCTION ---------------------------------------------
predict_krr <- function(Xnew) {
  K_new <- rbf_kernel(as.matrix(Xnew), X_tr_mat, sigma_hat)
  mu    <- K_new %*% alpha_hat
  data.frame(pred=as.numeric(mu),
             lo  = as.numeric(mu) - q_krr,
             hi  = as.numeric(mu) + q_krr)
}

# 9) EVALUATE ON TEST -----------------------------------------------------------
res_krr  <- predict_krr(X_te)
mse_krr  <- mean((res_krr$pred - y_te)^2)
cov_krr  <- mean(y_te >= res_krr$lo & y_te <= res_krr$hi)
wid_krr  <- mean(res_krr$hi - res_krr$lo)

# 10) BOGOSBINTED SHUFFLE TEST --------------------------------------------------
shuf_feat             <- features
shuf_feat$BogosBinted <- sample(shuf_feat$BogosBinted)

# reapply pipeline
Xsh_tr   <- predict(pp, predict(dumm, shuf_feat[i_tr,]))
Xsh_te   <- predict(pp, predict(dumm, shuf_feat[i_te,]))
Xsh_tr_m <- as.matrix(Xsh_tr)
Xsh_te_m <- as.matrix(Xsh_te)

# recompute α with same σ, λ
K_tr_sh   <- rbf_kernel(Xsh_tr_m, Xsh_tr_m, sigma_hat)
alpha_sh  <- solve(K_tr_sh + lambda_hat*diag(nrow(K_tr_sh)), y_tr)
pred_shuf <- (rbf_kernel(Xsh_te_m, Xsh_tr_m, sigma_hat) %*% alpha_sh)[,1]
mse_shuf  <- mean((pred_shuf - y_te)^2)

# 11) REPORT RESULTS -----------------------------------------------------------
cat("\n===== KRR (custom) =====\n")
cat(sprintf(" σ = %.4g   λ = %.4g\n", sigma_hat, lambda_hat))
cat(sprintf(" Test MSE      = %.3f\n",      mse_krr))
cat(sprintf(" Coverage(90%%) = %.1f%%\n",    100*cov_krr))
cat(sprintf(" Avg. Width     = %.3f\n\n",    wid_krr))

cat("===== BOGOSBINTED SHUFFLE =====\n")
cat(sprintf(" Test MSE after shuffle = %.3f\n\n", mse_shuf))

# 12) PLOT --------------------------------------------------------------------
plot_df <- cbind(res_krr, True=y_te)
ggplot(plot_df, aes(x=pred, y=True)) +
  geom_point(alpha=0.6) +
  geom_errorbar(aes(ymin=lo,ymax=hi), alpha=0.3) +
  labs(
    title    = "KRR + 90% Conformal on Test",
    subtitle = sprintf("MSE=%.3f | Cov=%.1f%% | W=%.2f",
                       mse_krr, 100*cov_krr, wid_krr),
    x="Predicted Happiness", y="True Happiness"
  ) +
  theme_minimal()

# 13) MARGINAL EFFECTS ---------------------------------------------------------
step   <- 1e-3
n_tr   <- nrow(X_tr_mat)
p      <- ncol(X_tr_mat)
meffs  <- matrix(0, n_tr, p)
colnames(meffs) <- colnames(X_tr_mat)

# We already have K_tr = rbf_kernel(X_tr_mat, X_tr_mat, sigma_hat)
K_tr <- rbf_kernel(X_tr_mat, X_tr_mat, sigma_hat)
for(j in seq_len(n_tr)) {
  # difference matrix: X_i - x_j
  diff_mat <- sweep(X_tr_mat, 2, X_tr_mat[j,], FUN="-")
  # for each i: ∂k(x_i,x_j)/∂x_j = k_ij * (x_j - x_i) / sigma^2
  # so ∂f/∂x_j = sum_i α_i * ∂k/∂x_j = -1/sigma^2 * sum_i α_i * K[i,j] * diff_mat[i,]
  # we take the *negative* of diff_mat because derivative wrt x_j
  meffs[j, ] <- colSums(alpha_hat * K_tr[,j] * (-diff_mat)) / sigma_hat^2
}

avg_me    <- colMeans(meffs)
quart_me  <- apply(meffs, 2, quantile, probs=c(0.25,0.5,0.75))

cat("\nAverage Marginal Effects:\n")
print(round(avg_me, 8))
cat("\nQuartiles of Marginal Effects:\n")
print(round(quart_me, 8))

# ───────────────────────────────────────────────────────────────────────────────
# End of Script

