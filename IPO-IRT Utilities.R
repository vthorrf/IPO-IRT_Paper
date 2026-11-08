### IPO-IRT: Utilities

### SECTION 1 - Dependencies====================================================
.required_packages <- c("statmod", "numDeriv", "imv")
for (.pkg in .required_packages) {
  if (!requireNamespace(.pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required. Install it with install.packages('%s').",
                 .pkg, .pkg))
  }
}
lapply(.required_packages, library, character.only = TRUE)

### SECTION 2 - Basic general utilities=========================================
`%||%` <- function(x, y) if (is.null(x)) y else x

.ipo_eps <- 1e-10

.ipo_normalize_mode <- function(mode) {
  mode <- match.arg(mode,
                    c("uniform", "relative_k", "k_relative",
                      "k_item", "k_free", "hw_free"))
  switch(mode, hw_free = "uniform",
         k_relative = "relative_k",
         k_item = "relative_k",
         k_free = "relative_k",
         mode)
}

.ipo_validate_Y <- function(Y, require_item_observed = TRUE) {
  Y <- as.matrix(Y)
  storage.mode(Y) <- "double"
  bad <- !is.na(Y) & !(Y %in% c(0, 1))
  if (any(bad)) stop("Y must contain only 0, 1, and NA values.")
  if (nrow(Y) < 1L || ncol(Y) < 1L) stop("Y must have at least one row and one column.")
  if (require_item_observed && any(colSums(!is.na(Y)) == 0L)) {
    stop("Every item must have at least one observed response.")
  }
  Y
}

.ipo_observed_matrices <- function(Y) {
  obs <- !is.na(Y)
  Y0 <- Y
  Y0[!obs] <- 0
  O <- matrix(as.numeric(obs), nrow(Y), ncol(Y))
  list(Y0 = Y0, O = O)
}

.ipo_logsumexp_rows <- function(x) {
  m <- apply(x, 1L, max)
  m + log(rowSums(exp(sweep(x, 1L, m, "-"))))
}

.ipo_safe_vcov <- function(H, tol = 1e-8) {
  if (is.null(H) || any(!is.finite(H))) return(NULL)
  H <- (H + t(H)) / 2
  ee <- tryCatch(eigen(H, symmetric = TRUE), error = function(e) NULL)
  if (is.null(ee)) return(NULL)
  cutoff <- tol * max(1, max(abs(ee$values)))
  if (any(ee$values <= cutoff)) return(NULL)
  ee$vectors %*% diag(1 / ee$values, nrow = length(ee$values)) %*% t(ee$vectors)
}

.ipo_quad <- function(n_quad) {
  statmod::gauss.quad.prob(n_quad, dist = "normal", mu = 0, sigma = 1)
}

### SECTION 3 - IPO state probabilities and ICC=================================
ipo_state_prob <- function(theta, b, hw, k) {
  #' IPO-IRT state probabilities
  #' @param theta Ability value(s).
  #' @param b Item difficulty, scalar or broadcastable object.
  #' @param hw Hesitation-band half-width, scalar or broadcastable object, > 0.
  #' @param k Boundary sharpness, scalar or broadcastable object, > 0.
  #' @return A list containing pi_ignorance, pi_overlap, and pi_mastery.

  D <- theta - b
  upper_cdf <- pnorm(k * (D + hw))
  lower_cdf <- pnorm(k * (D - hw))

  pi_ignorance <- 1 - upper_cdf
  pi_overlap   <- upper_cdf - lower_cdf
  pi_mastery   <- lower_cdf

  # Numerical guards only; the ordered-CDF construction is analytically valid.
  pi_ignorance <- pmax(pi_ignorance, 0)
  pi_overlap   <- pmax(pi_overlap, 0)
  pi_mastery   <- pmax(pi_mastery, 0)
  total <- pi_ignorance + pi_overlap + pi_mastery

  Results <- list(pi_ignorance = pi_ignorance / total,
                  pi_overlap = pi_overlap / total,
                  pi_mastery = pi_mastery / total)
  return(Results)
}

ipo_overlap_correct <- function(theta, b, alpha = 1, K = 2) {
  #' Correct-response probability conditional on the overlap state
  #' The use of a probit function here keeps the entire model smooth while avoiding
  #' the exact cancellation produced by combining normal-CDF state boundaries with
  #' a logistic overlap curve at particular half-width values.
  if (K <= 1) stop("K must be greater than 1.")
  if (!is.finite(alpha) || alpha <= 0) stop("alpha must be positive and finite.")
  c_guess <- 1 / K
  D <- theta - b
  c_guess + (1 - c_guess) * pnorm(alpha * D)
}

ipo_p <- function(theta, b, hw, k, alpha = 1, K = 2) {
  #' IPO-IRT item response probability
  #' @param theta Ability value(s).
  #' @param b Item difficulty.
  #' @param hw Hesitation-band half-width.
  #' @param k Boundary sharpness.
  #' @param alpha Fixed/global slope of the within-overlap response process.
  #' @param K Number of response alternatives.
  #' @return P(Y=1 | theta).
  
  st <- ipo_state_prob(theta, b, hw, k)
  c_guess <- 1 / K
  rho <- ipo_overlap_correct(theta, b, alpha = alpha, K = K)

  p <- c_guess * st$pi_ignorance + rho * st$pi_overlap + st$pi_mastery

  pmin(pmax(p, .ipo_eps), 1 - .ipo_eps)
}

.ipo_item_probability_matrix <- function(nodes, b, hw, k, alpha, K) {
  J <- length(b)
  Q <- length(nodes)
  theta_mat <- matrix(nodes, nrow = J, ncol = Q, byrow = TRUE)
  b_mat <- matrix(b, nrow = J, ncol = Q)
  hw_vec <- rep(hw, length.out = J)
  k_vec <- rep(k, length.out = J)
  hw_mat <- matrix(hw_vec, nrow = J, ncol = Q)
  k_mat <- matrix(k_vec, nrow = J, ncol = Q)
  ipo_p(theta_mat, b_mat, hw_mat, k_mat, alpha = alpha, K = K)
}

.ipo_loglik_iq <- function(Y, p_item) {
  mats <- .ipo_observed_matrices(Y)
  log_p <- log(p_item)
  log_1mp <- log1p(-p_item)
  mats$Y0 %*% log_p + (mats$O - mats$Y0) %*% log_1mp
}

### SECTION 4 - Parameterization, MML objective, and direct MML estimator======= 
.ipo_decode_relative_k <- function(v, J, k_fixed = 1) {
  if (J < 1L) stop("J must be positive.")
  if (!is.finite(k_fixed) || k_fixed <= 0) stop("k_fixed must be positive.")
  if (J == 1L) {
    u <- 0
  } else {
    if (length(v) != J - 1L) stop("v must have length J - 1.")
    u_raw <- c(v, 0)
    u <- u_raw - mean(u_raw)
  }
  list(u = u, k = k_fixed * exp(u))
}

.ipo_encode_relative_k <- function(k, k_fixed = 1) {
  k <- as.numeric(k)
  J <- length(k)
  if (J < 1L || any(!is.finite(k)) || any(k <= 0)) {
    stop("k must contain positive finite values.")
  }
  u <- log(k / k_fixed)
  u <- u - mean(u)
  if (J == 1L) numeric(0) else u[seq_len(J - 1L)] - u[J]
}

.ipo_k_boundary_penalty <- function(k, k_bounds, multiplier = 1e6) {
  log_k <- log(k)
  lower_violation <- pmax(log(k_bounds[1]) - log_k, 0)
  upper_violation <- pmax(log_k - log(k_bounds[2]), 0)
  multiplier * sum(lower_violation^2 + upper_violation^2)
}

.ipo_build_objective <- function(Y, mode, n_quad = 21, K = 2, k_fixed = 1,
                                 alpha_fixed = 1, start = NULL,
                                 hw_bounds = c(0.01, 4), k_bounds = c(0.05, 20),
                                 relative_logk_bounds = c(-4, 4)) {

  Y <- .ipo_validate_Y(Y)
  mode <- .ipo_normalize_mode(mode)
  J <- ncol(Y)
  gh <- .ipo_quad(n_quad)
  nodes <- gh$nodes
  weights <- gh$weights
  start <- start %||% list()

  if (length(hw_bounds) != 2L || hw_bounds[1] <= 0 || hw_bounds[2] <= hw_bounds[1]) {
    stop("hw_bounds must be an increasing positive vector of length two.")
  }
  if (length(k_bounds) != 2L || k_bounds[1] <= 0 || k_bounds[2] <= k_bounds[1]) {
    stop("k_bounds must be an increasing positive vector of length two.")
  }
  if (length(relative_logk_bounds) != 2L ||
      relative_logk_bounds[2] <= relative_logk_bounds[1]) {
    stop("relative_logk_bounds must be an increasing vector of length two.")
  }
  if (!is.finite(k_fixed) || k_fixed <= 0) stop("k_fixed must be positive.")
  if (!is.finite(alpha_fixed) || alpha_fixed <= 0) stop("alpha_fixed must be positive.")

  b0  <- rep(start$b %||% 0, length.out = J)
  hw0 <- rep(start$hw %||% 0.6, length.out = 1L)

  if (mode == "uniform") {
    par0 <- c(b0, log(hw0))
    lower <- c(rep(-Inf, J), log(hw_bounds[1]))
    upper <- c(rep( Inf, J), log(hw_bounds[2]))

    decode <- function(par) {
      list(b = par[seq_len(J)], hw = exp(par[J + 1L]),
           k = k_fixed, u = rep(0, J), alpha = alpha_fixed)
    }
    natural_vector <- function(par) {
      prm <- decode(par)
      c(prm$b, prm$hw)
    }
    natural_names <- c(paste0("b", seq_len(J)), "hw")
    raw_names     <- c(paste0("b", seq_len(J)), "log_hw")

  } else {
    k0    <- rep(start$k %||% k_fixed, length.out = J)
    v0    <- .ipo_encode_relative_k(k0, k_fixed = k_fixed)
    par0  <- c(b0, log(hw0), v0)
    lower <- c(rep(-Inf, J), log(hw_bounds[1]),
               rep(relative_logk_bounds[1], max(0L, J - 1L)))
    upper <- c(rep( Inf, J), log(hw_bounds[2]),
               rep(relative_logk_bounds[2], max(0L, J - 1L)))

    decode <- function(par) {
      v <- if (J == 1L) numeric(0) else par[J + 1L + seq_len(J - 1L)]
      rel <- .ipo_decode_relative_k(v, J = J, k_fixed = k_fixed)
      list(b = par[seq_len(J)], hw = exp(par[J + 1L]),
           k = rel$k, u = rel$u, alpha = alpha_fixed)
    }
    natural_vector <- function(par) {
      prm <- decode(par)
      c(prm$b, prm$hw, prm$k)
    }
    natural_names <- c(paste0("b", seq_len(J)), "hw", paste0("k", seq_len(J)))
    raw_names <- c(paste0("b", seq_len(J)), "log_hw",
                   if (J > 1L) paste0("relative_logk", seq_len(J - 1L)) else character(0))
  }

  if (length(par0) != length(lower) || any(!is.finite(par0))) {
    stop("Starting values are incomplete or non-finite for the selected mode.")
  }

  pure_neg_loglik <- function(par) {
    if (any(!is.finite(par))) return(1e100)
    prm <- decode(par)
    p_item <- .ipo_item_probability_matrix(nodes, prm$b, prm$hw,
                                           prm$k, prm$alpha, K)
    if (any(!is.finite(p_item))) return(1e100)
    loglik_iq <- .ipo_loglik_iq(Y, p_item)
    row_ll <- .ipo_logsumexp_rows(sweep(loglik_iq, 2L, log(weights), "+"))
    val <- -sum(row_ll)
    if (is.finite(val)) val else 1e100
  }

  neg_loglik <- function(par) {
    val <- pure_neg_loglik(par)
    if (mode == "relative_k" && is.finite(val)) {
      prm <- decode(par)
      val <- val + .ipo_k_boundary_penalty(prm$k, k_bounds)
    }
    if (is.finite(val)) val else 1e100
  }

  list(Y = Y,
       J = J,
       nodes = nodes,
       weights = weights,
       par0 = par0,
       lower = lower,
       upper = upper,
       decode = decode,
       natural_vector = natural_vector,
       natural_names = natural_names,
       raw_names = raw_names,
       neg_loglik = neg_loglik,
       pure_neg_loglik = pure_neg_loglik,
       mode = mode,
       K = K,
       k_fixed = k_fixed,
       alpha_fixed = alpha_fixed,
       hw_bounds = hw_bounds,
       k_bounds = k_bounds,
       relative_logk_bounds = relative_logk_bounds)
}

.ipo_attach_se <- function(result, objective, par, verbose = FALSE) {
  if (verbose) cat("Computing standard errors via numerical Hessian ...\n")
  H <- tryCatch(numDeriv::hessian(objective$pure_neg_loglik, par),
                error = function(e) NULL)
  vcov_raw <- .ipo_safe_vcov(H)
  if (is.null(vcov_raw)) {
    result$se <- NULL
    result$vcov <- NULL
    result$vcov_raw <- NULL
    result$hessian <- H
    result$hessian_ok <- FALSE
    result$hessian_condition <- Inf
    return(result)
  }

  G <- tryCatch(numDeriv::jacobian(objective$natural_vector, par),
                error = function(e) NULL)
  if (is.null(G) || any(!is.finite(G))) {
    result$se <- NULL
    result$vcov <- NULL
    result$vcov_raw <- vcov_raw
    result$hessian <- H
    result$hessian_ok <- FALSE
    result$hessian_condition <- Inf
    return(result)
  }

  vcov_nat <- G %*% vcov_raw %*% t(G)
  dimnames(vcov_nat) <- list(objective$natural_names, objective$natural_names)
  se_nat <- sqrt(pmax(diag(vcov_nat), 0))
  J <- objective$J

  if (objective$mode == "uniform") {
    result$se <- list(b = se_nat[seq_len(J)], hw = se_nat[J + 1L])
  } else {
    result$se <- list(b = se_nat[seq_len(J)],
                      hw = se_nat[J + 1L],
                      k = se_nat[J + 1L + seq_len(J)])
  }

  ee <- eigen((H + t(H)) / 2, symmetric = TRUE, only.values = TRUE)$values
  result$vcov <- vcov_nat
  result$vcov_raw <- vcov_raw
  result$hessian <- H
  result$hessian_eigenvalues <- ee
  result$hessian_condition <- max(ee) / min(ee)
  result$hessian_ok <- TRUE
  result
}

#' Fit an IPO-IRT model by direct marginal maximum likelihood
ipo_mml <- function(Y, mode = c("uniform", "relative_k"), n_quad = 21, K = 2,
                    k_fixed = 1, alpha_fixed = 1, start = NULL, maxit = 1000,
                    se = TRUE, verbose = TRUE, hw_bounds = c(0.01, 4),
                    k_bounds = c(0.05, 20), relative_logk_bounds = c(-4, 4)) {

  mode <- .ipo_normalize_mode(mode[1])
  obj <- .ipo_build_objective(Y = Y, mode = mode,
                              n_quad = n_quad,
                              K = K, k_fixed = k_fixed,
                              alpha_fixed = alpha_fixed,
                              start = start,
                              hw_bounds = hw_bounds,
                              k_bounds = k_bounds,
                              relative_logk_bounds = relative_logk_bounds)

  if (verbose) {
    cat(sprintf("Starting direct MML optimization (mode = %s) ...\n", mode))
  }

  fit <- optim(par = obj$par0,
               fn = obj$neg_loglik,
               method = "L-BFGS-B",
               lower = obj$lower,
               upper = obj$upper,
               control = list(maxit = maxit, trace = as.integer(verbose), REPORT = 50))

  prm <- obj$decode(fit$par)
  result <- list(b = prm$b,
                 hw = prm$hw,
                 k = prm$k,
                 u = prm$u,
                 k_fixed = k_fixed,
                 k_geometric_mean = exp(mean(log(prm$k))),
                 alpha = prm$alpha,
                 K = K,
                 c = 1 / K,
                 mu = 0,
                 sigma = 1,
                 mode = mode,
                 loglik = -obj$pure_neg_loglik(fit$par),
                 n_par = length(fit$par),
                 convergence = fit$convergence,
                 converged = identical(fit$convergence, 0L),
                 message = fit$message,
                 counts = fit$counts,
                 raw_par = fit$par,
                 fit = fit,
                 call = match.call())

  if (se) result <- .ipo_attach_se(result, obj, fit$par, verbose = verbose)
  class(result) <- "ipo_MML"
  result
}

### SECTION 5 - EM / generalized-EM estimator=================================== 
.ipo_e_step <- function(Y, nodes, weights, b, hw, k, alpha, K) {
  p_item <- .ipo_item_probability_matrix(nodes, b, hw, k, alpha, K)
  loglik_iq <- .ipo_loglik_iq(Y, p_item)
  logpost <- sweep(loglik_iq, 2L, log(weights), "+")
  row_ll <- .ipo_logsumexp_rows(logpost)
  m <- apply(logpost, 1L, max)
  r <- exp(sweep(logpost, 1L, m, "-"))
  r <- r / rowSums(r)
  list(r = r, loglik = sum(row_ll), person_loglik = row_ll, p_item = p_item)
}

ipo_em <- function(Y, mode = c("uniform", "relative_k"), n_quad = 21, K = 2,
                   k_fixed = 1, alpha_fixed = 1, start = NULL, max_iter = 300,
                   mstep_maxit = 500, tol = 1e-6, se = TRUE, verbose = TRUE,
                   hw_bounds = c(0.01, 4), k_bounds = c(0.05, 20),
                   relative_logk_bounds = c(-4, 4)) {
  #' Fit an IPO-IRT model by EM / generalized EM
  #' Model U has a pooled common-width M-step. Model K uses a joint M-step because
  #' the sum-to-zero relative-log-sharpness constraint couples the items.

  Y <- .ipo_validate_Y(Y)
  mode <- .ipo_normalize_mode(mode[1])
  mats <- .ipo_observed_matrices(Y)
  obj <- .ipo_build_objective(Y = Y, mode = mode, n_quad = n_quad, K = K,
                              k_fixed = k_fixed, alpha_fixed = alpha_fixed, start = start,
                              hw_bounds = hw_bounds, k_bounds = k_bounds,
                              relative_logk_bounds = relative_logk_bounds)

  par_current <- obj$par0
  prm <- obj$decode(par_current)
  current <- .ipo_e_step(Y, obj$nodes, obj$weights, prm$b,
                         prm$hw, prm$k, prm$alpha, K)
  ll <- current$loglik
  converged <- FALSE
  last_mstep <- NULL

  for (iter in seq_len(max_iter)) {
    r <- current$r
    R <- t(r) %*% mats$Y0
    Nobs <- t(r) %*% mats$O
    W <- Nobs - R

    neg_Q <- function(par) {
      prm_c <- obj$decode(par)
      p <- .ipo_item_probability_matrix(
        obj$nodes, prm_c$b, prm_c$hw, prm_c$k, prm_c$alpha, K
      )
      val <- -sum(R * t(log(p)) + W * t(log1p(-p)))
      if (mode == "relative_k" && is.finite(val)) {
        val <- val + .ipo_k_boundary_penalty(prm_c$k, k_bounds)
      }
      if (is.finite(val)) val else 1e100
    }

    last_mstep <- optim(par = par_current,
                        fn = neg_Q,
                        method = "L-BFGS-B",
                        lower = obj$lower,
                        upper = obj$upper,
                        control = list(maxit = mstep_maxit))
    candidate_par <- last_mstep$par
    candidate_prm <- obj$decode(candidate_par)
    updated <- .ipo_e_step(Y, obj$nodes, obj$weights,
                           candidate_prm$b, candidate_prm$hw, candidate_prm$k,
                           candidate_prm$alpha, K)
    ll_new <- updated$loglik
    change <- ll_new - ll

    if (verbose) {
      cat(sprintf(
        "EM iter %3d | log-lik = %.6f | change = %+0.6e | M-step code = %d\n",
        iter, ll_new, change, last_mstep$convergence
      ))
    }

    # Numerical optimizers may occasionally return a tiny decrease. Retain the
    # previous iterate rather than violating the EM ascent property.
    if (!is.finite(ll_new) || change < -1e-7) {
      if (verbose) cat("M-step did not improve the observed likelihood; stopping.\n")
      break
    }

    par_current <- candidate_par
    prm <- candidate_prm
    current <- updated
    ll <- ll_new

    if (abs(change) < tol) {
      converged <- TRUE
      break
    }
  }

  result <- list(b = prm$b,
                 hw = prm$hw,
                 k = prm$k,
                 u = prm$u,
                 k_fixed = k_fixed,
                 k_geometric_mean = exp(mean(log(prm$k))),
                 alpha = prm$alpha,
                 K = K,
                 c = 1 / K,
                 mu = 0,
                 sigma = 1,
                 mode = mode,
                 loglik = ll,
                 n_par = length(par_current),
                 iter = iter,
                 converged = converged,
                 convergence = if (converged) 0L else 1L,
                 message = if (converged) "EM convergence criterion satisfied" else
                   "EM stopped before satisfying the convergence criterion",
                 mstep = last_mstep,
                 raw_par = par_current,
                 call = match.call())

  if (se) result <- .ipo_attach_se(result, obj, par_current, verbose = verbose)
  class(result) <- "ipo_EM"
  result
}

### SECTION 6 - EAP scoring, prediction, and marginal held-out likelihood======= 
#' Expected-a-posteriori ability estimates
ipo_eap <- function(Y, fit, n_quad = 31, K = fit$K %||% 2) {
  if (!inherits(fit, c("ipo_MML", "ipo_EM"))) {
    stop("fit must be an object returned by ipo_mml() or ipo_em().")
  }
  Y <- .ipo_validate_Y(Y, require_item_observed = FALSE)
  if (ncol(Y) != length(fit$b)) stop("Y and fit contain different numbers of items.")

  gh <- .ipo_quad(n_quad)
  estep <- .ipo_e_step(Y, gh$nodes, gh$weights,
                       b = fit$b, hw = fit$hw, k = fit$k,
                       alpha = fit$alpha, K = K)
  post <- estep$r
  theta_eap <- as.numeric(post %*% gh$nodes)
  theta_var <- as.numeric(post %*% (gh$nodes^2)) - theta_eap^2

  data.frame(theta_eap = theta_eap,
             theta_sd = sqrt(pmax(theta_var, 0)),
             n_observed = rowSums(!is.na(Y)))
}

#' Predicted probabilities using EAP scores based on all supplied responses
ipo_predict <- function(Y, fit, n_quad = 31, K = fit$K %||% 2) {
  scores <- ipo_eap(Y, fit, n_quad = n_quad, K = K)
  N <- nrow(Y)
  J <- ncol(Y)
  P <- matrix(NA_real_, N, J)
  for (j in seq_len(J)) {
    k_j <- if (length(fit$k) == 1L) fit$k else fit$k[j]
    P[, j] <- ipo_p(scores$theta_eap, fit$b[j], fit$hw, k_j,
                    alpha = fit$alpha, K = K)
  }
  P
}

ipo_predict_loo <- function(Y, fit, n_quad = 31, K = fit$K %||% 2) {
  #' Leave-one-item-out cell predictions
  #' For item j, the response to j is masked before the respondent's EAP score is
  #' computed. This prevents the target response from contributing to its own
  #' prediction and is the recommended input for pairwise cross-validated IMV.
  Y <- .ipo_validate_Y(Y, require_item_observed = FALSE)
  N <- nrow(Y)
  J <- ncol(Y)
  P <- matrix(NA_real_, N, J)
  for (j in seq_len(J)) {
    Y_minus_j <- Y
    Y_minus_j[, j] <- NA_real_
    scores <- ipo_eap(Y_minus_j, fit, n_quad = n_quad, K = K)
    k_j <- if (length(fit$k) == 1L) fit$k else fit$k[j]
    P[, j] <- ipo_p(scores$theta_eap, fit$b[j], fit$hw, k_j,
                    alpha = fit$alpha, K = K)
  }
  P
}

#' Marginal person log likelihood under a fitted IPO model
ipo_marginal_person_loglik <- function(Y, fit, n_quad = 31, K = fit$K %||% 2) {
  Y <- .ipo_validate_Y(Y, require_item_observed = FALSE)
  gh <- .ipo_quad(n_quad)
  p_item <- .ipo_item_probability_matrix(gh$nodes, fit$b, fit$hw,
                                         fit$k, fit$alpha, K)
  loglik_iq <- .ipo_loglik_iq(Y, p_item)
  .ipo_logsumexp_rows(sweep(loglik_iq, 2L, log(gh$weights), "+"))
}

ipo_marginal_loglik <- function(Y, fit, n_quad = 31, K = fit$K %||% 2) {
  sum(ipo_marginal_person_loglik(Y, fit, n_quad = n_quad, K = K))
}

### SECTION 7 - IPO data simulation============================================= 
ipo_balanced_difficulties <- function(J, lower = -2, upper = 2) {
  if (J < 1L) stop("J must be positive.")
  if (J == 1L) return(0)
  seq(lower, upper, length.out = J)
}

#' Simulate responses from Model U or the relative-sharpness Model K
ipo_simulate <- function(N = 500, J = 20, mode = c("uniform", "relative_k"),
                         b = NULL, hw = 0.6, k = NULL, k_fixed = 1,
                         alpha_fixed = 1, k_sdlog = 0.35, K = 2,
                         missing_rate = 0, balanced_b = FALSE, seed = NULL) {

  mode <- .ipo_normalize_mode(mode[1])
  if (!is.null(seed)) set.seed(seed)
  if (is.null(b)) {
    b <- if (balanced_b) ipo_balanced_difficulties(J) else rnorm(J, 0, 1)
  }
  if (length(b) != J) stop("b must have length J.")
  if (length(hw) != 1L || !is.finite(hw) || hw <= 0) {
    stop("hw must be one positive finite common half-width.")
  }

  theta <- rnorm(N)

  if (mode == "uniform") {
    k <- k_fixed
    u <- rep(0, J)
  } else {
    if (is.null(k)) {
      u <- rnorm(J, mean = 0, sd = k_sdlog)
      u <- u - mean(u)
      k <- k_fixed * exp(u)
    } else {
      k <- rep(k, length.out = J)
      if (any(!is.finite(k)) || any(k <= 0)) stop("k must be positive and finite.")
      u <- log(k / k_fixed)
      u <- u - mean(u)
      k <- k_fixed * exp(u)
    }
  }

  Y <- matrix(NA_real_, N, J)
  P <- matrix(NA_real_, N, J)
  for (j in seq_len(J)) {
    k_j <- if (length(k) == 1L) k else k[j]
    P[, j] <- ipo_p(theta, b[j], hw, k_j, alpha_fixed, K)
    Y[, j] <- rbinom(N, size = 1, prob = P[, j])
  }

  if (missing_rate < 0 || missing_rate >= 1) {
    stop("missing_rate must be in [0, 1).")
  }
  if (missing_rate > 0) {
    miss <- matrix(runif(N * J) < missing_rate, N, J)
    Y[miss] <- NA_real_
  }

  list(Y = Y,
       P = P,
      theta = theta,
      b = b,
      hw = hw,
      k = k,
      u = u,
      k_fixed = k_fixed,
      k_geometric_mean = exp(mean(log(rep(k, length.out = J)))),
      alpha = alpha_fixed,
      K = K,
      mode = mode,
      missing_rate = missing_rate)
}

### SECTION 8 - Model comparison and pairwise InterModel Vigorish=============== 
compare_irt_ic <- function(..., criterion = c("AIC", "BIC", "EBIC", "HQ"),
                           nobs = NULL, ebic_gamma = 0.5, ebic_p = NULL,
                           sort = TRUE, check_same_nobs = TRUE) {
  
  criterion <- match.arg(toupper(criterion[1L]),
                         choices = c("AIC", "BIC", "EBIC", "HQ"))
  
  models <- list(...)
  
  if (length(models) < 2L) {
    stop("Supply at least two fitted models.")
  }
  
  # ---------------------------------------------------------------------------
  # 1. Recover informative model names
  # ---------------------------------------------------------------------------
  model_expressions <- as.list(substitute(list(...)))[-1L]
  model_names <- names(models)
  
  if (is.null(model_names)) {
    model_names <- rep("", length(models))
  }
  
  for (i in seq_along(models)) {
    if (is.na(model_names[i]) || model_names[i] == "") {
      model_names[i] <- paste(deparse(model_expressions[[i]], width.cutoff = 80L),
                              collapse = "")
    }
  }
  
  model_names <- make.unique(model_names)
  
  # ---------------------------------------------------------------------------
  # 2. Resolve scalar, positional, or named arguments
  # ---------------------------------------------------------------------------
  resolve_argument <- function(x, model_name, index, total, argument_name) {
    if (is.null(x)) {
      return(NULL)
    }
    
    if (length(x) == 1L) {
      return(as.numeric(x))
    }
    
    if (!is.null(names(x)) && model_name %in% names(x)) {
      return(as.numeric(x[[model_name]]))
    }
    
    if (length(x) == total) {
      return(as.numeric(x[index]))
    }
    
    stop(argument_name,
         " must be NULL, a scalar, a vector with one value per model, ",
         "or a named vector matching the model names.")
  }
  
  # ---------------------------------------------------------------------------
  # 3. Identify mirt objects
  # ---------------------------------------------------------------------------
  is_mirt_model <- function(object) {
    inherits(object, c("SingleGroupClass",
                       "MultipleGroupClass",
                       "MixedClass",
                       "DiscreteGroupClass",
                       "MixtureClass"))
  }
  
  # ---------------------------------------------------------------------------
  # 4. Extract the sample size from a mirt object
  # ---------------------------------------------------------------------------
  extract_mirt_nobs <- function(object) {
    raw_data <- tryCatch(mirt::extract.mirt(object, "data"),
                         error = function(e) NULL)
    
    if (!is.null(raw_data) && length(dim(raw_data)) == 2L) {
      return(nrow(raw_data))
    }
    
    frequencies <- tryCatch(mirt::extract.mirt(object, "freq"),
                            error = function(e) NULL)
    
    if (!is.null(frequencies) &&
        length(frequencies) > 0L &&
        all(is.finite(frequencies))) {
      return(sum(frequencies))
    }
    
    row_ids <- tryCatch(mirt::extract.mirt(object, "rowID"),
                        error = function(e) NULL)
    
    if (!is.null(row_ids)) {
      return(length(row_ids))
    }
    
    NA_real_
  }
  
  # ---------------------------------------------------------------------------
  # 5. Extract logLik, number of parameters, and N from one model
  # ---------------------------------------------------------------------------
  extract_model_information <- function(object, model_name, index,
                                        total, nobs_argument) {
    
    supplied_nobs <- resolve_argument(x = nobs_argument, model_name = model_name,
                                      index = index, total = total,
                                      argument_name = "nobs")
    
    # -------------------------------------------------------------------------
    # IPO-IRT models from ipo_mml() or ipo_em()
    # -------------------------------------------------------------------------
    if (inherits(object, c("ipo_MML", "ipo_EM"))) {
      
      log_likelihood <- object$loglik %||% object$logLik
      number_parameters <- object$n_par %||%
        object$npar %||%
        object$df
      
      model_nobs <- supplied_nobs %||%
        object$nobs %||%
        object$N %||%
        object$n_obs
      
      if (is.null(log_likelihood) ||
          length(log_likelihood) != 1L ||
          !is.finite(log_likelihood)) {
        stop("Could not extract a finite log-likelihood from IPO model '",
             model_name, "'.")
      }
      
      if (is.null(number_parameters) ||
          length(number_parameters) != 1L ||
          !is.finite(number_parameters)) {
        stop("Could not extract the number of parameters from IPO model '",
             model_name, "'.")
      }
      
      if (is.null(model_nobs)) {
        model_nobs <- NA_real_
      }
      
      return(list(model_class = class(object)[1L],
                  loglik = as.numeric(log_likelihood),
                  npar = as.integer(number_parameters),
                  nobs = as.numeric(model_nobs)))
    }
    
    # -------------------------------------------------------------------------
    # Logistic models fitted by logistic_mml() in the IPO code
    # -------------------------------------------------------------------------
    if (inherits(object, "logistic_MML")) {
      
      model_nobs <- supplied_nobs %||%
        object$nobs %||%
        object$N %||%
        object$n_obs
      
      if (is.null(model_nobs)) {
        model_nobs <- NA_real_
      }
      
      return(list(model_class = "logistic_MML",
                  loglik = as.numeric(object$loglik),
                  npar = as.integer(object$n_par),
                  nobs = as.numeric(model_nobs)))
    }
    
    # -------------------------------------------------------------------------
    # mirt models
    # -------------------------------------------------------------------------
    if (is_mirt_model(object)) {
      if (!requireNamespace("mirt", quietly = TRUE)) {
        stop("The mirt package is required to inspect mirt models.")
      }
      
      log_likelihood <- mirt::extract.mirt(object, "logLik")
      number_parameters <- mirt::extract.mirt(object, "nestpars")
      
      model_nobs <- supplied_nobs
      
      if (is.null(model_nobs)) {
        model_nobs <- extract_mirt_nobs(object)
      }
      
      return(list(model_class = class(object)[1L],
                  loglik = as.numeric(log_likelihood),
                  npar = as.integer(number_parameters),
                  nobs = as.numeric(model_nobs)))
    }
    
    # -------------------------------------------------------------------------
    # Generic fallback for objects with a logLik() method
    # -------------------------------------------------------------------------
    loglik_object <- tryCatch(stats::logLik(object),
                              error = function(e) NULL)
    
    if (!is.null(loglik_object)) {
      log_likelihood <- as.numeric(loglik_object)
      number_parameters <- attr(loglik_object, "df")
      
      model_nobs <- supplied_nobs
      
      if (is.null(model_nobs)) {
        model_nobs <- attr(loglik_object, "nobs")
      }
      
      if (is.null(model_nobs)) {
        model_nobs <- tryCatch(stats::nobs(object),
                               error = function(e) NA_real_)
      }
      
      if (is.null(number_parameters) || !is.finite(number_parameters)) {
        stop("The logLik() method for model '", model_name,
             "' did not report the number of free parameters.")
      }
      
      return(list(model_class = class(object)[1L],
                  loglik = log_likelihood,
                  npar = as.integer(number_parameters),
                  nobs = as.numeric(model_nobs)))
    }
    
    stop("Unsupported model '", model_name,
         "' with class: ", paste(class(object), collapse = ", "), ".")
  }
  
  # Local null-coalescing operator
  `%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0L) y else x
  }
  
  # ---------------------------------------------------------------------------
  # 6. Extract information from all supplied models
  # ---------------------------------------------------------------------------
  extracted <- lapply(seq_along(models), function(i) {
      extract_model_information(object = models[[i]], model_name = model_names[i],
                                index = i, total = length(models),
                                nobs_argument = nobs)})
  
  comparison <- data.frame(model = model_names,
                           class = vapply(extracted, function(x) x$model_class,
                                          character(1)),
                           logLik = vapply(extracted, function(x) x$loglik,
                                           numeric(1)),
                           n_parameters = vapply(extracted, function(x) x$npar,
                                                 integer(1)),
                           nobs = vapply(extracted, function(x) x$nobs,
                                         numeric(1)),
                           stringsAsFactors = FALSE)
  
  # ---------------------------------------------------------------------------
  # 7. Validate comparability
  # ---------------------------------------------------------------------------
  if (any(!is.finite(comparison$logLik))) {
    stop("At least one model has a non-finite log-likelihood.")
  }
  
  if (any(comparison$n_parameters < 0L)) {
    stop("The number of model parameters cannot be negative.")
  }
  
  criterion_requires_n <- criterion %in% c("BIC", "EBIC", "HQ")
  
  if (criterion_requires_n &&
      any(!is.finite(comparison$nobs))) {
    stop(criterion,
         " requires the sample size. Supply nobs explicitly for any IPO ",
         "fit that does not store it, for example nobs = nrow(training_data).")
  }
  
  if (criterion_requires_n && any(comparison$nobs <= 1)) {
    stop("Sample sizes must be greater than one for ", criterion, ".")
  }
  
  if (check_same_nobs) {
    finite_nobs <- comparison$nobs[is.finite(comparison$nobs)]
    
    if (length(unique(finite_nobs)) > 1L) {
      warning("The models report different sample sizes. Information criteria ",
              "should ordinarily be compared only when models were fitted to ",
              "the same observations.")
    }
  }
  
  # ---------------------------------------------------------------------------
  # 8. Compute the requested information criterion
  # ---------------------------------------------------------------------------
  k <- comparison$n_parameters
  ll <- comparison$logLik
  n <- comparison$nobs
  
  IC <- switch(criterion,
               AIC = {-2 * ll + 2 * k},
               BIC = {-2 * ll + log(n) * k},
               HQ = {-2 * ll + 2 * log(log(n)) * k},
               EBIC = {if (is.null(ebic_p)) {
                 stop("For EBIC, supply ebic_p: the total number of candidate ",
                      "parameters or features in the common model universe.")
                 }
                 if (length(ebic_p) != 1L ||
                     !is.finite(ebic_p) ||
                     ebic_p <= 0 ||
                     ebic_p != round(ebic_p)) {
                 stop("ebic_p must be one positive integer.")
                 }
                 if (ebic_p < max(k)) {
                 stop("ebic_p must be at least as large as the greatest number ",
                      "of free parameters among the candidate models.")
                 }
                 if (length(ebic_gamma) != 1L ||
                     !is.finite(ebic_gamma) ||
                     ebic_gamma < 0) {
                 stop("ebic_gamma must be one nonnegative finite number.")
                 }
                 # Chen and Chen form:
                 # EBIC = BIC + 2 * gamma * log(choose(P, k))
                 -2 * ll + log(n) * k + 2 * ebic_gamma * lchoose(ebic_p, k)
                 })
  
  # ---------------------------------------------------------------------------
  # 9. Delta IC, relative likelihoods, and model weights
  # ---------------------------------------------------------------------------
  delta_IC <- IC - min(IC)
  
  # IC-based likelihood ratio relative to the best candidate.
  # The best model therefore has a value of 1.
  likelihood_ratio <- exp(-0.5 * delta_IC)
  
  model_weight <- likelihood_ratio / sum(likelihood_ratio)
  
  comparison$criterion <- criterion
  comparison$IC <- IC
  comparison$delta_IC <- delta_IC
  comparison$likelihood_ratio <- likelihood_ratio
  comparison$model_weight <- model_weight
  
  if (sort) {
    comparison <- comparison[order(comparison$IC),,drop = FALSE]
  }
  rownames(comparison) <- NULL
  attr(comparison, "criterion") <- criterion
  attr(comparison, "ebic_gamma") <- if (criterion == "EBIC") {
    ebic_gamma
  } else {
    NULL
  }
  attr(comparison, "ebic_p") <- if (criterion == "EBIC") {
    ebic_p
  } else {
    NULL
  }
  class(comparison) <- c("irt_ic_comparison", class(comparison))
  comparison
}

imv_compare_irt <- function(train_data, test_data,
                            fit_baseline, fit_enhanced,
                            baseline_name = "Baseline",
                            enhanced_name = "Enhanced",
                            prediction = c("loo", "eap"),
                            n_quad = 31, sigma = 1e-4) {
  
  prediction <- match.arg(prediction)
  
  # ---------------------------------------------------------------------------
  # 1. Validate and align the datasets
  # ---------------------------------------------------------------------------
  as_binary_matrix <- function(x, object_name) {
    x <- as.matrix(x)
    
    if (!is.numeric(x) && !is.integer(x) && !is.logical(x)) {
      stop(object_name, " must be a numeric, integer, or logical matrix.")
    }
    
    storage.mode(x) <- "double"
    
    bad <- !is.na(x) & !(x %in% c(0, 1))
    if (any(bad)) {
      stop(object_name, " must contain only 0, 1, and NA.")
    }
    
    x
  }
  
  train_data <- as_binary_matrix(train_data, "train_data")
  test_data  <- as_binary_matrix(test_data,  "test_data")
  
  if (ncol(train_data) != ncol(test_data)) {
    stop("train_data and test_data must contain the same number of items.")
  }
  
  # Reorder test columns to match training columns when names are available.
  if (!is.null(colnames(train_data)) && !is.null(colnames(test_data))) {
    if (!setequal(colnames(train_data), colnames(test_data))) {
      stop("train_data and test_data do not contain the same named items.")
    }
    
    test_data <- test_data[, colnames(train_data), drop = FALSE]
  }
  
  J <- ncol(test_data)
  
  # ---------------------------------------------------------------------------
  # 2. Identify the model type
  # ---------------------------------------------------------------------------
  model_type <- function(fit) {
    if (inherits(fit, c("ipo_MML", "ipo_EM"))) {
      return("ipo")
    }
    
    if (inherits(fit, "SingleGroupClass")) {
      return("mirt")
    }
    
    stop("Unsupported model class: ",
         paste(class(fit), collapse = ", "),
         ". Fits must be IPO models or mirt SingleGroupClass models.")
  }
  
  baseline_type <- model_type(fit_baseline)
  enhanced_type <- model_type(fit_enhanced)
  
  # ---------------------------------------------------------------------------
  # 3. IPO predictions
  # ---------------------------------------------------------------------------
  predict_ipo <- function(Y, fit, prediction, n_quad) {
    if (length(fit$b) != ncol(Y)) {
      stop("The IPO fit and test_data contain different numbers of items.")
    }
    
    if (prediction == "loo") {
      if (!exists("ipo_predict_loo", mode = "function")) {
        stop("ipo_predict_loo() was not found. ",
             "Source the complete IPO-IRT code before using this function.")
      }
      
      P <- ipo_predict_loo(Y = Y, fit = fit, n_quad = n_quad)
      
    } else {
      if (!exists("ipo_predict", mode = "function")) {
        stop("ipo_predict() was not found. ",
             "Source the complete IPO-IRT code before using this function.")
      }
      
      P <- ipo_predict(Y = Y, fit = fit, n_quad = n_quad)
    }
    
    as.matrix(P)
  }
  
  # ---------------------------------------------------------------------------
  # 4. mirt predictions
  # ---------------------------------------------------------------------------
  predict_mirt <- function(Y, fit, prediction, n_quad) {
    if (!requireNamespace("mirt", quietly = TRUE)) {
      stop("The mirt package is required for predicting from mirt models.")
    }
    
    n_items_fit <- mirt::extract.mirt(fit, "nitems")
    
    if (n_items_fit != ncol(Y)) {
      stop("The mirt fit and test_data contain different numbers of items.")
    }
    
    categories <- mirt::extract.mirt(fit, "K")
    
    if (any(categories != 2L)) {
      stop("The current IMV function supports only dichotomous mirt items. ",
           "At least one item in the supplied mirt model is polytomous.")
    }
    
    n_factors <- mirt::extract.mirt(fit, "nfact")
    fit_item_names <- mirt::extract.mirt(fit, "itemnames")
    
    # Put the response matrix in the same item order used by the mirt fit.
    if (!is.null(colnames(Y)) && !is.null(fit_item_names)) {
      if (!setequal(colnames(Y), fit_item_names)) {
        stop("The item names in test_data do not match the item names ",
             "in the mirt model.")
      }
      
      Y_fit_order <- Y[, fit_item_names, drop = FALSE]
    } else {
      Y_fit_order <- Y
    }
    
    N <- nrow(Y_fit_order)
    P_fit_order <- matrix(NA_real_, nrow = N, ncol = n_items_fit,
                          dimnames = list(NULL, fit_item_names))
    
    # Score a response matrix and return the factor-score matrix.
    estimate_theta <- function(response_matrix) {
      Theta <- matrix(0, nrow = nrow(response_matrix), ncol = n_factors)
      
      # If all responses are missing, the EAP equals the prior mean, zero.
      scorable <- rowSums(!is.na(response_matrix)) > 0L
      
      if (any(scorable)) {
        scores <- mirt::fscores(fit, method = "EAP",
                                response.pattern = response_matrix[scorable, ,
                                                                   drop = FALSE],
                                full.scores = TRUE, full.scores.SE = FALSE,
                                quadpts = n_quad, verbose = FALSE)
        
        scores <- as.matrix(scores)
        
        Theta[scorable, ] <- scores[, seq_len(n_factors), drop = FALSE]
      }
      
      Theta
    }
    
    if (prediction == "eap") {
      # Descriptive plug-in predictions. Each response contributes to its
      # respondent's EAP score.
      Theta <- estimate_theta(Y_fit_order)
      
      for (j in seq_len(n_items_fit)) {
        item_object <- mirt::extract.item(fit, j)
        trace <- as.matrix(mirt::probtrace(item_object, Theta))
        
        if (ncol(trace) != 2L) {
          stop("Item ", j, " did not produce two probability traces.")
        }
        
        # The final column corresponds to the highest category, coded as 1.
        P_fit_order[, j] <- trace[, ncol(trace)]
      }
      
    } else {
      # Leave-one-item-out predictions. Response j is masked before EAP
      # scoring and then predicted from the remaining items.
      for (j in seq_len(n_items_fit)) {
        Y_minus_j <- Y_fit_order
        Y_minus_j[, j] <- NA_real_
        
        Theta_minus_j <- estimate_theta(Y_minus_j)
        
        item_object <- mirt::extract.item(fit, j)
        trace <- as.matrix(mirt::probtrace(item_object, Theta_minus_j))
        
        if (ncol(trace) != 2L) {
          stop("Item ", j, " did not produce two probability traces.")
        }
        
        P_fit_order[, j] <- trace[, ncol(trace)]
      }
    }
    
    # Return probabilities in the same order as test_data.
    if (!is.null(colnames(Y)) && !is.null(fit_item_names)) {
      P_fit_order[, colnames(Y), drop = FALSE]
    } else {
      P_fit_order
    }
  }
  
  # ---------------------------------------------------------------------------
  # 5. Obtain predictions from both models
  # ---------------------------------------------------------------------------
  
  get_predictions <- function(Y, fit, type) {
    switch(type, ipo = predict_ipo(Y = Y, fit = fit, prediction = prediction,
                                   n_quad = n_quad),
      mirt = predict_mirt(Y = Y, fit = fit, prediction = prediction,
                          n_quad = n_quad))
  }
  
  P_baseline <- get_predictions(test_data, fit_baseline, baseline_type)
  
  P_enhanced <- get_predictions(test_data, fit_enhanced, enhanced_type)
  
  if (!identical(dim(P_baseline), dim(test_data)) ||
      !identical(dim(P_enhanced), dim(test_data))) {
    stop("The predicted-probability matrices have incorrect dimensions.")
  }
  
  # ---------------------------------------------------------------------------
  # 6. Construct inputs for imv.binary()
  # ---------------------------------------------------------------------------
  observed <- as.vector(test_data)
  p_baseline <- as.vector(P_baseline)
  p_enhanced <- as.vector(P_enhanced)
  
  keep <- (!is.na(observed) & is.finite(p_baseline) & is.finite(p_enhanced))
  
  m0 <- observed[keep]
  m1 <- p_baseline[keep]
  p2 <- p_enhanced[keep]
  
  if (length(m0) == 0L) {
    stop("No jointly valid observed responses and predictions were found.")
  }
  
  # imv.binary() performs its own probability clipping through sigma.
  if (!requireNamespace("imv", quietly = TRUE)) {
    stop("The imv package is required. Install it with ",
         "install.packages('imv').")
  }
  
  imv_value <- imv::imv.binary(y = m0, p0 = m1, p1 = p2, sigma = sigma)
  
  favours <- if (imv_value > 0) {
    enhanced_name
  } else if (imv_value < 0) {
    baseline_name
  } else {
    "Neither model"
  }
  
  # ---------------------------------------------------------------------------
  # 7. Return the IMV and the exact vectors supplied to imv.binary()
  # ---------------------------------------------------------------------------
  result <- list(IMV = as.numeric(imv_value),
                 favours = favours,
                 baseline_model = baseline_name,
                 enhanced_model = enhanced_name,
                 prediction = prediction,
                 n_test_responses = length(m0),
                 # Exact inputs passed to imv.binary()
                 y  = m0,
                 p0 = m1,
                 p1 = p2,
                 # Matrices can be useful for item-level diagnostics.
                 baseline_probability_matrix = P_baseline,
                 enhanced_probability_matrix = P_enhanced)
  result
}

### SECTION 9 - Print methods and example usage================================= 
print.ipo_MML <- function(x, ...) {
  cat(sprintf("IPO-IRT MML fit: mode = %s\n", x$mode))
  cat(sprintf("logLik = %.3f | n_par = %d | converged = %s | alpha = %.3f\n",
              x$loglik, x$n_par, x$converged, x$alpha))
  if (x$mode == "uniform") {
    cat(sprintf("common hw = %.4f | fixed k = %.4f\n", x$hw, x$k))
    print(data.frame(item = seq_along(x$b), b = x$b), row.names = FALSE)
  } else {
    cat(sprintf("common hw = %.4f | geometric mean k = %.4f (fixed anchor %.4f)\n",
                x$hw, x$k_geometric_mean, x$k_fixed))
    print(data.frame(item = seq_along(x$b), b = x$b, k = x$k, u = x$u),
          row.names = FALSE)
  }
  invisible(x)
}

print.ipo_EM <- function(x, ...) {
  cat(sprintf("IPO-IRT EM fit: mode = %s\n", x$mode))
  cat(sprintf("logLik = %.3f | n_par = %d | converged = %s | iterations = %d | alpha = %.3f\n",
              x$loglik, x$n_par, x$converged, x$iter, x$alpha))
  if (x$mode == "uniform") {
    cat(sprintf("common hw = %.4f | fixed k = %.4f\n", x$hw, x$k))
    print(data.frame(item = seq_along(x$b), b = x$b), row.names = FALSE)
  } else {
    cat(sprintf("common hw = %.4f | geometric mean k = %.4f (fixed anchor %.4f)\n",
                x$hw, x$k_geometric_mean, x$k_fixed))
    print(data.frame(item = seq_along(x$b), b = x$b, k = x$k, u = x$u),
          row.names = FALSE)
  }
  invisible(x)
}

###########################====---- THE END ----====############################