## Shared data-generating-process helpers for tests, mirroring the DGPs
## in codes/run_internal_tilted_sec4_1.R (Sec 4.1 of the reference
## paper) and codes/run_external_tilted_sec4_2.R (Sec 4.2).

sim_internal_dgp <- function(n = 400) {
  true_g_fn <- function(X1, X2, X3) stats::plogis(1.0 * X1 - 0.8 * X2 + 0.8 * X3)
  true_Q_fn <- function(X1, X2, X3, A) {
    gX <- true_g_fn(X1, X2, X3)
    stats::plogis(-0.4 + 0.7 * X1 - 0.5 * X2 + 0.6 * X3 + A * (-0.2 + 1.0 * gX))
  }
  X1 <- stats::rnorm(n)
  X2 <- stats::rnorm(n)
  X3 <- stats::rbinom(n, 1, 0.5)
  gX <- true_g_fn(X1, X2, X3)
  A <- stats::rbinom(n, 1, gX)
  Q <- true_Q_fn(X1, X2, X3, A)
  Y <- stats::rbinom(n, 1, Q)
  list(Y = Y, A = A, W = data.frame(X1 = X1, X2 = X2, X3 = X3))
}

approximate_true_tau_internal <- function(eta, n_large = 2e5, seed = 20270101) {
  true_g_fn <- function(X1, X2, X3) stats::plogis(1.0 * X1 - 0.8 * X2 + 0.8 * X3)
  true_Q_fn <- function(X1, X2, X3, A) {
    gX <- true_g_fn(X1, X2, X3)
    stats::plogis(-0.4 + 0.7 * X1 - 0.5 * X2 + 0.6 * X3 + A * (-0.2 + 1.0 * gX))
  }
  set.seed(seed)
  X1 <- stats::rnorm(n_large)
  X2 <- stats::rnorm(n_large)
  X3 <- stats::rbinom(n_large, 1, 0.5)
  gX <- true_g_fn(X1, X2, X3)
  Q1 <- true_Q_fn(X1, X2, X3, 1)
  Q0 <- true_Q_fn(X1, X2, X3, 0)
  eta_g <- eta(gX)
  sum(eta_g * (Q1 - Q0)) / sum(eta_g)
}

sim_external_dgp <- function(n = 400, m = 800, N_super = 5000) {
  true_s_fn <- function(X1, X2, X3) stats::plogis(1.5 * X1 - 1.5 * X2 + X3)
  true_g_fn <- function(X1, X2, X3) stats::plogis(-0.2 + 0.8 * X1 - 0.6 * X2)
  true_Q_fn <- function(X1, X2, X3, A) {
    sX <- true_s_fn(X1, X2, X3)
    stats::plogis(-0.4 + 0.7 * X1 - 0.4 * X2 + 0.5 * X3 + A * (-0.20 + 2.00 * (1 - sX)))
  }
  repeat {
    X1 <- stats::rnorm(N_super)
    X2 <- stats::rnorm(N_super)
    X3 <- stats::rbinom(N_super, 1, 0.5)
    sX <- true_s_fn(X1, X2, X3)
    D_pool <- stats::rbinom(N_super, 1, sX)
    if (sum(D_pool == 1) >= n && sum(D_pool == 0) >= m) break
  }
  idx1 <- sample(which(D_pool == 1), n)
  idx0 <- sample(which(D_pool == 0), m)
  X1_int <- X1[idx1]
  X2_int <- X2[idx1]
  X3_int <- X3[idx1]
  gX <- true_g_fn(X1_int, X2_int, X3_int)
  A <- stats::rbinom(n, 1, gX)
  Q <- true_Q_fn(X1_int, X2_int, X3_int, A)
  Y <- stats::rbinom(n, 1, Q)
  Ws <- data.frame(X1 = c(X1[idx1], X1[idx0]), X2 = c(X2[idx1], X2[idx0]), X3 = c(X3[idx1], X3[idx0]))
  D <- c(rep(1L, n), rep(0L, m))
  list(Y = Y, A = A, D = D, Ws = Ws)
}

approximate_true_tau_external <- function(eta, n_large = 2e5, seed = 20270102) {
  true_s_fn <- function(X1, X2, X3) stats::plogis(1.5 * X1 - 1.5 * X2 + X3)
  true_Q_fn <- function(X1, X2, X3, A) {
    sX <- true_s_fn(X1, X2, X3)
    stats::plogis(-0.4 + 0.7 * X1 - 0.4 * X2 + 0.5 * X3 + A * (-0.20 + 2.00 * (1 - sX)))
  }
  set.seed(seed)
  X1 <- stats::rnorm(n_large)
  X2 <- stats::rnorm(n_large)
  X3 <- stats::rbinom(n_large, 1, 0.5)
  sX <- true_s_fn(X1, X2, X3)
  Q1 <- true_Q_fn(X1, X2, X3, 1)
  Q0 <- true_Q_fn(X1, X2, X3, 0)
  eta_s <- eta(sX)
  sum(sX * eta_s * (Q1 - Q0)) / sum(sX * eta_s)
}
