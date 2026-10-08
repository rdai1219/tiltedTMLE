library(tiltedTMLE)

## --- Internal-tilted smoke test (mirrors run_internal_tilted_sec4_1.R DGP) ---
true_g_fn <- function(X1, X2, X3) plogis(1.0 * X1 - 0.8 * X2 + 0.8 * X3)
true_Q_fn <- function(X1, X2, X3, A) {
  gX <- true_g_fn(X1, X2, X3)
  plogis(-0.4 + 0.7 * X1 - 0.5 * X2 + 0.6 * X3 + A * (-0.2 + 1.0 * gX))
}
set.seed(2027)
n <- 400
X1 <- rnorm(n); X2 <- rnorm(n); X3 <- rbinom(n, 1, 0.5)
gX <- true_g_fn(X1, X2, X3)
A <- rbinom(n, 1, gX)
Q <- true_Q_fn(X1, X2, X3, A)
Y <- rbinom(n, 1, Q)
W <- data.frame(X1 = X1, X2 = X2, X3 = X3)

cat("=== Internal-tilted (ATE/ATT/ATC/ATO) ===\n")
for (est in c("ATE", "ATT", "ATC", "ATO")) {
  spec <- tilt_spec(est)
  out <- estimate_tau_h_internal_tilted(Y, A, W, eta = spec$eta, eta_dot = spec$eta_dot,
                                         Q_method = "glm", Qform = Y ~ A * (X1 + X2 + X3),
                                         g_method = "glm", gform = A ~ X1 + X2 + X3,
                                         truncate = FALSE)
  cat(sprintf("  %-4s TMLE=%.6f SE=%.6f iter=%d\n", est, out$TMLE, out$TMLE_se, out$TMLE_iter))
}

## --- External-tilted smoke test (mirrors run_external_tilted_sec4_2.R DGP) ---
true_s_fn <- function(X1, X2, X3) plogis(1.5 * X1 - 1.5 * X2 + X3)
true_g_fn2 <- function(X1, X2, X3) plogis(-0.2 + 0.8 * X1 - 0.6 * X2)
true_Q_fn2 <- function(X1, X2, X3, A) {
  sX <- true_s_fn(X1, X2, X3)
  plogis(-0.4 + 0.7 * X1 - 0.4 * X2 + 0.5 * X3 + A * (-0.20 + 2.00 * (1 - sX)))
}
simulate_dgp <- function(n = 400, m = 800, N_super = 5000) {
  repeat {
    X1 <- rnorm(N_super); X2 <- rnorm(N_super); X3 <- rbinom(N_super, 1, 0.5)
    sX <- true_s_fn(X1, X2, X3)
    D_pool <- rbinom(N_super, 1, sX)
    if (sum(D_pool == 1) >= n && sum(D_pool == 0) >= m) break
  }
  idx1 <- sample(which(D_pool == 1), n)
  idx0 <- sample(which(D_pool == 0), m)
  X1_int <- X1[idx1]; X2_int <- X2[idx1]; X3_int <- X3[idx1]
  gX <- true_g_fn2(X1_int, X2_int, X3_int)
  A <- rbinom(n, 1, gX)
  Q <- true_Q_fn2(X1_int, X2_int, X3_int, A)
  Y <- rbinom(n, 1, Q)
  Ws <- data.frame(X1 = c(X1[idx1], X1[idx0]), X2 = c(X2[idx1], X2[idx0]), X3 = c(X3[idx1], X3[idx0]))
  D <- c(rep(1L, n), rep(0L, m))
  list(Y = Y, A = A, D = D, Ws = Ws)
}

set.seed(4242)
d <- simulate_dgp(400, 800)

cat("\n=== External-tilted (ATE-I/ATE-E/ATE-O) ===\n")
for (est in c("ATE-I", "ATE-E", "ATE-O")) {
  spec <- tilt_spec_external(est)
  out <- estimate_tau_h_external_tilted(
    d$Y, d$A, d$D, d$Ws, eta = spec$eta, eta_dot = spec$eta_dot,
    Q_method = "glm", Qform = Y ~ A * (X1 + X2 + X3),
    g_method = "glm", gform = A ~ X1 + X2 + X3,
    s_method = "glm", sform = A ~ X1 + X2 + X3,
    truncate_g = TRUE, gbounds = c(0.025, 0.975),
    truncate_s = TRUE, sbounds = c(0.025, 0.975)
  )
  cat(sprintf("  %-6s TMLE=%.6f SE=%.6f iter=%d\n", est, out$TMLE, out$TMLE_se, out$TMLE_iter))
}

cat("\n=== estimate_tau_h_fixed (fixed/precomputed h, ATE) ===\n")
out_fixed <- estimate_tau_h_fixed(Y, A, W, h = rep(1, length(Y)),
                                   Q_method = "glm", Qform = Y ~ A * (X1 + X2 + X3),
                                   g_method = "glm", gform = A ~ X1 + X2 + X3,
                                   truncate = FALSE, n_boot = 50)
cat(sprintf("  TMLE=%.6f SE=%.6f AIPW=%.6f\n", out_fixed$TMLE, out_fixed$TMLE_se, out_fixed$AIPW))

cat("\nAll smoke tests completed without error.\n")
