### IPO-IRT Simulation Study

### SECTION 1 - Dependencies and source file ===================================
.required_packages <- c("mirt", "imv", "future", "future.apply", "ggplot2")
for (.pkg in .required_packages) {
  if (!requireNamespace(.pkg, quietly = TRUE)) {
    stop(sprintf(
      "Package '%s' is required. Install it with install.packages('%s').",
      .pkg, .pkg
    ))
  }
}

UTILITIES_FILE <- "IPO-IRT Utilities.R"
if (!exists("ipo_mml", mode = "function")) {
  if (!file.exists(UTILITIES_FILE)) {
    stop(
      "Could not find '", UTILITIES_FILE, "'. Set UTILITIES_FILE to the path ",
      "of the current IPO-IRT Utilities.R file."
    )
  }
  source(UTILITIES_FILE)
}

### SECTION 2 - Simulation configuration =======================================
# Core crossed design: N is the model-fitting sample size. Each replication also
# generates an independent test sample for theta recovery and predictive IMV.
N_VALUES <- c(250L, 500L, 1000L, 2000L)
J_VALUES <- c(10L, 20L, 50L, 100L)
DGP_IPO  <- c("U", "K")

N_REPS       <- 1000L
N_TEST       <- 1000L
BASE_SEED    <- 20260807L

# Base IPO values. These are NOT crossed in the final design.
HW_TRUE      <- 0.60
K_FIXED      <- 1.00
ALPHA_FIXED  <- 1.00
K_SDLOG      <- 0.35      # fixed heterogeneity of relative log-k under IPO-K
K_ALT        <- 2L        # number of answer alternatives; outcomes remain 0/1

# Item difficulty generation. Set TRUE to reduce nuisance variation across
# replications; FALSE reproduces the earlier N(0,1) random-difficulty design.
BALANCED_DIFFICULTIES <- FALSE

# Logistic DGP parameters used in the large-sample misspecification arm.
LOGISTIC_A_SDLOG <- 0.25
THREE_PL_C_CONCENTRATION <- 20

# Numerical settings.
IPO_N_QUAD    <- 21L
PRED_GRID     <- seq(-6, 6, length.out = 61L)
HW_BOUNDS     <- c(0.01, 4)
HW_BOUND_TOL  <- 1e-6
IMV_SIGMA     <- 1e-4
IC_CRITERIA   <- c("AIC", "BIC", "HQ")
IMV_PROTOCOLS <- c("pp", "eap")  # posterior predictive is the primary protocol

# Models fitted under IPO truth. 3PL is reserved for the large logistic arm to
# avoid unnecessary computation and unstable small-N 3PL fits.
MODELS_MAIN <- c("IPO-U", "IPO-K", "Rasch", "2PL", "2.5PL")
MODELS_LOGISTIC_ARM <- c("IPO-U", "IPO-K", "Rasch", "2PL", "2.5PL", "3PL")

# Optional supplement-style IMV protocol sweep. It deliberately includes very
# short forms because the EAP-plugin vs posterior-predictive discrepancy is a
# short-test phenomenon.
PROTOCOL_J_VALUES <- c(3L, 5L, 8L, 12L, 20L, 40L)
PROTOCOL_N_TRAIN  <- 4000L
PROTOCOL_N_TEST   <- 4000L
PROTOCOL_N_REPS   <- 50L

# Output.
RESULTS_DIR <- "IPO_IRT_Simulation_Results"

### SECTION 3 - General helpers ================================================
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

.empty_df <- function() data.frame(stringsAsFactors = FALSE)

.is_fit_error <- function(x) inherits(x, "simulation_fit_error")

.safe_fit <- function(expr) {
  tryCatch(
    expr,
    error = function(e) {
      structure(
        list(message = conditionMessage(e)),
        class = "simulation_fit_error"
      )
    }
  )
}

.safe_numeric <- function(expr, default = NA_real_) {
  tryCatch(as.numeric(expr), error = function(e) default)
}

.make_item_names <- function(J) paste0("X", seq_len(J))

.make_difficulties <- function(J) {
  if (isTRUE(BALANCED_DIFFICULTIES)) {
    ipo_balanced_difficulties(J)
  } else {
    rnorm(J, 0, 1)
  }
}

.metric_summary <- function(estimate, truth) {
  estimate <- as.numeric(estimate)
  truth <- as.numeric(truth)
  keep <- is.finite(estimate) & is.finite(truth)
  if (!any(keep)) {
    return(c(bias = NA_real_, rmse = NA_real_, cor = NA_real_, n = 0))
  }
  estimate <- estimate[keep]
  truth <- truth[keep]
  d <- estimate - truth
  rr <- if (length(estimate) >= 2L && stats::sd(estimate) > 0 && stats::sd(truth) > 0) {
    suppressWarnings(stats::cor(estimate, truth))
  } else {
    NA_real_
  }
  c(bias = mean(d), rmse = sqrt(mean(d^2)),
    cor = rr, n = length(d))
}

.fit_is_converged <- function(fit) {
  if (.is_fit_error(fit)) return(FALSE)
  if (inherits(fit, c("ipo_MML", "ipo_EM"))) {
    return(isTRUE(fit$converged))
  }
  if (inherits(fit, "SingleGroupClass")) {
    return(tryCatch(isTRUE(fit@OptimInfo$converged), error = function(e) NA))
  }
  NA
}

.ipo_hw_at_bound <- function(fit) {
  if (.is_fit_error(fit) || !inherits(fit, c("ipo_MML", "ipo_EM"))) return(FALSE)
  any(abs(fit$hw - HW_BOUNDS) <= HW_BOUND_TOL)
}

.model_info <- function(fit) {
  if (.is_fit_error(fit)) {
    return(list(loglik = NA_real_, npar = NA_integer_))
  }
  if (inherits(fit, c("ipo_MML", "ipo_EM"))) {
    return(list(loglik = as.numeric(fit$loglik), npar = as.integer(fit$n_par)))
  }
  if (inherits(fit, "SingleGroupClass")) {
    return(list(
      loglik = .safe_numeric(mirt::extract.mirt(fit, "logLik")),
      npar = as.integer(.safe_numeric(mirt::extract.mirt(fit, "nestpars")))
    ))
  }
  list(loglik = NA_real_, npar = NA_integer_)
}

### SECTION 4 - Data-generating processes ======================================
.simulate_ipo_pair <- function(N_train, N_test, J, dgp, seed) {
  dgp <- match.arg(dgp, c("U", "K"))
  set.seed(seed)

  b <- .make_difficulties(J)
  if (dgp == "U") {
    u <- rep(0, J)
    k <- rep(K_FIXED, J)
  } else {
    u <- rnorm(J, 0, K_SDLOG)
    u <- u - mean(u)
    k <- K_FIXED * exp(u)
  }

  theta_train <- rnorm(N_train)
  theta_test  <- rnorm(N_test)

  generate_Y <- function(theta) {
    N <- length(theta)
    P <- matrix(NA_real_, N, J)
    Y <- matrix(NA_real_, N, J)
    for (j in seq_len(J)) {
      P[, j] <- ipo_p(
        theta = theta,
        b = b[j],
        hw = HW_TRUE,
        k = k[j],
        alpha = ALPHA_FIXED,
        K = K_ALT
      )
      Y[, j] <- rbinom(N, 1L, P[, j])
    }
    colnames(Y) <- colnames(P) <- .make_item_names(J)
    list(Y = Y, P = P)
  }

  tr <- generate_Y(theta_train)
  te <- generate_Y(theta_test)

  list(
    train = tr$Y,
    test = te$Y,
    theta_train = theta_train,
    theta_test = theta_test,
    b = b,
    hw = HW_TRUE,
    k = k,
    u = u,
    dgp = paste0("IPO-", dgp),
    K = K_ALT
  )
}

.simulate_logistic_pair <- function(N_train, N_test, J, dgp, seed) {
  dgp <- match.arg(dgp, c("Rasch", "2PL", "3PL"))
  set.seed(seed)

  b <- .make_difficulties(J)
  a <- if (dgp == "Rasch") rep(1, J) else exp(rnorm(J, 0, LOGISTIC_A_SDLOG))

  if (dgp == "3PL") {
    c_mean <- 1 / K_ALT
    phi <- THREE_PL_C_CONCENTRATION
    c_item <- stats::rbeta(J, shape1 = c_mean * phi, shape2 = (1 - c_mean) * phi)
    c_item <- pmin(pmax(c_item, 1e-4), 0.49)  # stays inside common 3PL guessing bounds
  } else {
    c_item <- rep(0, J)
  }

  theta_train <- rnorm(N_train)
  theta_test  <- rnorm(N_test)

  generate_Y <- function(theta) {
    N <- length(theta)
    eta <- outer(theta, a) - matrix(a * b, N, J, byrow = TRUE)
    base <- plogis(eta)
    P <- matrix(c_item, N, J, byrow = TRUE) +
      (1 - matrix(c_item, N, J, byrow = TRUE)) * base
    Y <- matrix(rbinom(N * J, 1L, as.vector(P)), N, J)
    colnames(Y) <- colnames(P) <- .make_item_names(J)
    list(Y = Y, P = P)
  }

  tr <- generate_Y(theta_train)
  te <- generate_Y(theta_test)

  list(
    train = tr$Y,
    test = te$Y,
    theta_train = theta_train,
    theta_test = theta_test,
    b = b,
    a = a,
    c = c_item,
    dgp = dgp,
    K = K_ALT
  )
}

### SECTION 5 - Model fitting ===================================================
.fit_models <- function(Y, model_names) {
  Y <- as.matrix(Y)
  storage.mode(Y) <- "double"
  J <- ncol(Y)
  if (is.null(colnames(Y))) colnames(Y) <- .make_item_names(J)

  fits <- list()

  if ("IPO-U" %in% model_names) {
    fits[["IPO-U"]] <- .safe_fit(ipo_mml(Y, mode = "uniform", n_quad = IPO_N_QUAD,
                                         K = K_ALT, k_fixed = K_FIXED,
                                         alpha_fixed = ALPHA_FIXED,
                                         se = FALSE, verbose = FALSE,
                                         hw_bounds = HW_BOUNDS))
  }

  if ("IPO-K" %in% model_names) {
    fits[["IPO-K"]] <- .safe_fit(ipo_mml(Y, mode = "relative_k", n_quad = IPO_N_QUAD,
                                         K = K_ALT, k_fixed = K_FIXED,
                                         alpha_fixed = ALPHA_FIXED,
                                         se = FALSE, verbose = FALSE,
                                         hw_bounds = HW_BOUNDS))
  }

  if ("Rasch" %in% model_names) {
    fits[["Rasch"]] <- .safe_fit(suppressMessages(mirt::mirt(data.frame(Y), 1, itemtype = "Rasch",
                                                             SE = FALSE, verbose = FALSE)))
  }

  if ("2PL" %in% model_names) {
    fits[["2PL"]] <- .safe_fit(suppressMessages(mirt::mirt(data.frame(Y), 1, itemtype = "2PL",
                                                           SE = FALSE, verbose = FALSE)))
  }

  if ("2.5PL" %in% model_names) {
    fits[["2.5PL"]] <- .safe_fit(suppressMessages(mirt::mirt(data.frame(Y), 1, itemtype = "2PL",
                                                             guess = 1 / K_ALT,
                                                             SE = FALSE, verbose = FALSE)))
  }

  if ("3PL" %in% model_names) {
    fits[["3PL"]] <- .safe_fit(suppressMessages(mirt::mirt(data.frame(Y), 1, itemtype = "3PL",
                                                           SE = FALSE, verbose = FALSE)))
  }

  # Return in requested order.
  fits[intersect(model_names, names(fits))]
}

### SECTION 6 - Unified latent-posterior prediction machinery ==================
# Following the coauthor protocol, every fitted model is reduced to an ICC on a
# common theta grid and a model-specific latent-distribution weight vector.
# This ensures that pp vs eap differences are prediction-protocol differences,
# not differences in how IPO and mirt models are scored.
.model_on_grid <- function(fit, J, theta_grid = PRED_GRID) {
  if (.is_fit_error(fit)) stop("Cannot construct predictions from a failed fit.")

  if (inherits(fit, c("ipo_MML", "ipo_EM"))) {
    icc <- sapply(seq_len(J), function(j) {
      k_j <- if (length(fit$k) == 1L) fit$k else fit$k[j]
      ipo_p(theta_grid, fit$b[j], fit$hw, k_j, fit$alpha, fit$K)
    })
    w <- dnorm(theta_grid)  # IPO fixes theta ~ N(0,1)

  } else if (inherits(fit, "SingleGroupClass")) {
    icc <- sapply(seq_len(J), function(j) {
      tr <- mirt::probtrace(mirt::extract.item(fit, j),
                            matrix(theta_grid, ncol = 1L))
      tr[, ncol(tr)]
    })

    gp <- mirt::coef(fit, simplify = TRUE)
    m <- tryCatch(as.numeric(gp$means[1]), error = function(e) 0)
    s <- tryCatch(sqrt(as.numeric(gp$cov[1, 1])), error = function(e) 1)
    if (!is.finite(m)) m <- 0
    if (!is.finite(s) || s <= 0) s <- 1
    w <- dnorm(theta_grid, mean = m, sd = s)

  } else {
    stop("Unsupported fitted-model class for prediction.")
  }

  icc <- as.matrix(icc)
  if (nrow(icc) != length(theta_grid)) icc <- matrix(icc, nrow = length(theta_grid))
  if (ncol(icc) != J) stop("ICC grid has the wrong number of items.")
  icc <- pmin(pmax(icc, 1e-10), 1 - 1e-10)
  w <- w / sum(w)

  list(icc = icc, w = w, theta_grid = theta_grid)
}

.full_posterior_scores <- function(Y, mg) {
  Y <- as.matrix(Y)
  O <- (!is.na(Y)) * 1
  Y0 <- Y
  Y0[is.na(Y0)] <- 0

  log_icc <- log(mg$icc)
  log_1m  <- log1p(-mg$icc)
  liq <- Y0 %*% t(log_icc) + (O - Y0) %*% t(log_1m)
  lp <- sweep(liq, 2L, log(mg$w), "+")
  lp <- lp - apply(lp, 1L, max)
  post <- exp(lp)
  post <- post / rowSums(post)

  mu <- as.numeric(post %*% mg$theta_grid)
  v <- as.numeric(post %*% (mg$theta_grid^2)) - mu^2

  data.frame(theta_eap = mu, theta_sd = sqrt(pmax(v, 0)),
             n_observed = rowSums(!is.na(Y)))
}

.loo_predict_both <- function(Y, mg) {
  Y <- as.matrix(Y)
  N <- nrow(Y)
  J <- ncol(Y)

  icc <- mg$icc
  lw <- log(mg$w)
  O <- (!is.na(Y)) * 1
  Y0 <- Y
  Y0[is.na(Y0)] <- 0

  log_icc <- log(icc)
  log_1m <- log1p(-icc)
  liq <- Y0 %*% t(log_icc) + (O - Y0) %*% t(log_1m)

  P_pp <- P_eap <- matrix(NA_real_, N, J)
  post_sd_item <- numeric(J)

  for (j in seq_len(J)) {
    contrib <- outer(Y0[, j], log_icc[, j]) +
      outer(O[, j] - Y0[, j], log_1m[, j])

    lp <- liq - contrib + matrix(lw, N, length(mg$theta_grid), byrow = TRUE)
    lp <- lp - apply(lp, 1L, max)
    post <- exp(lp)
    post <- post / rowSums(post)

    P_pp[, j] <- as.numeric(post %*% icc[, j])
    theta_eap <- as.numeric(post %*% mg$theta_grid)
    second <- as.numeric(post %*% (mg$theta_grid^2))
    post_sd_item[j] <- mean(sqrt(pmax(second - theta_eap^2, 0)))
    P_eap[, j] <- approx(mg$theta_grid, icc[, j], xout = theta_eap,
                         rule = 2)$y
  }

  list(pp = P_pp, eap = P_eap,
       mean_loo_post_sd = mean(post_sd_item),
       mean_signed_protocol_gap = mean(P_pp - P_eap, na.rm = TRUE),
       mean_abs_protocol_gap = mean(abs(P_pp - P_eap), na.rm = TRUE))
}

# Positional arguments avoid the p0/p1 versus p1/p2 argument-name difference
# across imv package versions.
.imv_binary_compat <- function(y, p_baseline, p_enhanced, sigma = IMV_SIGMA) {
  y <- as.numeric(y)
  p_baseline <- as.numeric(p_baseline)
  p_enhanced <- as.numeric(p_enhanced)
  keep <- is.finite(y) & is.finite(p_baseline) & is.finite(p_enhanced)
  if (!any(keep)) return(NA_real_)
  tryCatch(as.numeric(imv::imv.binary(y[keep], p_baseline[keep], p_enhanced[keep],
                                      sigma = sigma)),
           error = function(e) NA_real_)
}

### SECTION 7 - Recovery and model-comparison extractors =======================
.fit_health_table <- function(fits, meta) {
  out <- lapply(names(fits), function(nm) {
    ft <- fits[[nm]]
    info <- .model_info(ft)
    data.frame(arm = meta$arm,
               dgp = meta$dgp,
               N = meta$N,
               J = meta$J,
               rep = meta$rep,
               seed = meta$seed,
               model = nm,
               fit_ok = !.is_fit_error(ft),
               converged = .fit_is_converged(ft),
               logLik = info$loglik,
               n_parameters = info$npar,
               hw = if (inherits(ft, c("ipo_MML", "ipo_EM"))) ft$hw else NA_real_,
               hw_at_bound = .ipo_hw_at_bound(ft),
               error = if (.is_fit_error(ft)) ft$message else NA_character_,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

.parameter_recovery_table <- function(fits, sim, meta) {
  rows <- list()
  add_row <- function(model, parameter, estimate, truth) {
    met <- .metric_summary(estimate, truth)
    rows[[length(rows) + 1L]] <<- data.frame(arm = meta$arm,
                                             dgp = meta$dgp,
                                             N = meta$N,
                                             J = meta$J,
                                             rep = meta$rep,
                                             seed = meta$seed,
                                             model = model,
                                             correctly_specified = (meta$dgp == "IPO-U" && model == "IPO-U") ||
                                                                   (meta$dgp == "IPO-K" && model == "IPO-K"),
                                             parameter = parameter,
                                             true_mean = mean(as.numeric(truth), na.rm = TRUE),
                                             estimate_mean = mean(as.numeric(estimate), na.rm = TRUE),
                                             bias = met["bias"],
                                             rmse = met["rmse"],
                                             correlation = met["cor"],
                                             n_parameter_values = as.integer(met["n"]),
                                             hw_at_bound = .ipo_hw_at_bound(fits[[model]]),
                                             stringsAsFactors = FALSE)
  }

  for (model in intersect(c("IPO-U", "IPO-K"), names(fits))) {
    ft <- fits[[model]]
    if (.is_fit_error(ft)) next
    add_row(model, "b", ft$b, sim$b)
    add_row(model, "hw", ft$hw, sim$hw)

    if (model == "IPO-K") {
      true_k <- sim$k %||% rep(K_FIXED, meta$J)
      true_u <- sim$u %||% rep(0, meta$J)
      add_row(model, "k", ft$k, true_k)
      add_row(model, "u_log_relative_k", ft$u, true_u)
    }
  }

  if (!length(rows)) .empty_df() else do.call(rbind, rows)
}

.theta_recovery_table <- function(fits, Y_test, theta_true, meta) {
  rows <- list()
  for (nm in names(fits)) {
    ft <- fits[[nm]]
    if (.is_fit_error(ft)) next

    mg <- tryCatch(.model_on_grid(ft, ncol(Y_test)), error = function(e) NULL)
    if (is.null(mg)) next
    sc <- tryCatch(.full_posterior_scores(Y_test, mg), error = function(e) NULL)
    if (is.null(sc)) next

    met <- .metric_summary(sc$theta_eap, theta_true)
    rows[[length(rows) + 1L]] <- data.frame(arm = meta$arm,
                                            dgp = meta$dgp,
                                            N = meta$N,
                                            J = meta$J,
                                            rep = meta$rep,
                                            seed = meta$seed,
                                            model = nm,
                                            theta_bias = met["bias"],
                                            theta_rmse = met["rmse"],
                                            theta_correlation = met["cor"],
                                            mean_theta_estimate = mean(sc$theta_eap),
                                            mean_theta_true = mean(theta_true),
                                            mean_posterior_sd = mean(sc$theta_sd),
                                            hw_at_bound = .ipo_hw_at_bound(ft),
                                            stringsAsFactors = FALSE)
  }
  if (!length(rows)) .empty_df() else do.call(rbind, rows)
}

.ic_table <- function(fits, nobs, meta) {
  valid <- fits[!vapply(fits, .is_fit_error, logical(1))]
  if (length(valid) < 2L) return(.empty_df())

  base <- do.call(rbind, lapply(names(valid), function(nm) {
    info <- .model_info(valid[[nm]])
    data.frame(model = nm,
               logLik = info$loglik,
               n_parameters = info$npar,
               hw_at_bound = .ipo_hw_at_bound(valid[[nm]]),
               stringsAsFactors = FALSE)
  }))
  base <- base[is.finite(base$logLik) & is.finite(base$n_parameters), , drop = FALSE]
  if (nrow(base) < 2L) return(.empty_df())

  rows <- list()
  for (criterion in IC_CRITERIA) {
    IC <- switch(criterion,
                 AIC = -2 * base$logLik + 2 * base$n_parameters,
                 BIC = -2 * base$logLik + log(nobs) * base$n_parameters,
                 HQ  = -2 * base$logLik + 2 * log(log(nobs)) * base$n_parameters)
    delta <- IC - min(IC)
    rel <- exp(-0.5 * delta)
    weight <- rel / sum(rel)

    rows[[length(rows) + 1L]] <- data.frame(arm = meta$arm,
                                            dgp = meta$dgp,
                                            N = meta$N,
                                            J = meta$J,
                                            rep = meta$rep,
                                            seed = meta$seed,
                                            model = base$model,
                                            criterion = criterion,
                                            IC = IC,
                                            deltaIC = delta,
                                            relative_likelihood = rel,
                                            model_weight = weight,
                                            selected = delta == 0,
                                            hw_at_bound = base$hw_at_bound,
                                            stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

.lr_u_k_table <- function(fits, meta) {
  if (!all(c("IPO-U", "IPO-K") %in% names(fits))) return(.empty_df())
  U <- fits[["IPO-U"]]
  K <- fits[["IPO-K"]]
  if (.is_fit_error(U) || .is_fit_error(K)) return(.empty_df())

  iu <- .model_info(U)
  ik <- .model_info(K)
  if (!all(is.finite(c(iu$loglik, ik$loglik, iu$npar, ik$npar)))) return(.empty_df())

  stat <- 2 * (ik$loglik - iu$loglik)
  df <- ik$npar - iu$npar
  data.frame(arm = meta$arm,
             dgp = meta$dgp,
             N = meta$N,
             J = meta$J,
             rep = meta$rep,
             seed = meta$seed,
             comparison = "IPO-U_vs_IPO-K",
             LR = stat,
             df = df,
             p_value = if (df > 0 && is.finite(stat) && stat >= 0) {
               stats::pchisq(stat, df = df, lower.tail = FALSE)
             } else {
               NA_real_
             },
             U_hw_at_bound = .ipo_hw_at_bound(U),
             K_hw_at_bound = .ipo_hw_at_bound(K),
             stringsAsFactors = FALSE)
}

.pairwise_imv_tables <- function(fits, Y_test, meta) {
  valid <- fits[!vapply(fits, .is_fit_error, logical(1))]
  if (length(valid) < 2L) {
    return(list(imv = .empty_df(), prediction = .empty_df()))
  }

  pred <- list()
  pred_rows <- list()
  for (nm in names(valid)) {
    mg <- tryCatch(.model_on_grid(valid[[nm]], ncol(Y_test)), error = function(e) NULL)
    if (is.null(mg)) next
    pp <- tryCatch(.loo_predict_both(Y_test, mg), error = function(e) NULL)
    if (is.null(pp)) next
    pred[[nm]] <- pp
    pred_rows[[length(pred_rows) + 1L]] <- data.frame(arm = meta$arm,
                                                      dgp = meta$dgp,
                                                      N = meta$N,
                                                      J = meta$J,
                                                      rep = meta$rep,
                                                      seed = meta$seed,
                                                      model = nm,
                                                      mean_loo_posterior_sd = pp$mean_loo_post_sd,
                                                      mean_signed_pp_minus_eap = pp$mean_signed_protocol_gap,
                                                      mean_abs_pp_minus_eap = pp$mean_abs_protocol_gap,
                                                      hw_at_bound = .ipo_hw_at_bound(valid[[nm]]),
                                                      stringsAsFactors = FALSE)
  }

  if (length(pred) < 2L) {
    return(list(
      imv = .empty_df(),
      prediction = if (length(pred_rows)) do.call(rbind, pred_rows) else .empty_df()
    ))
  }

  pairs <- combn(names(pred), 2L, simplify = FALSE)
  y <- as.vector(Y_test)
  rows <- list()

  for (p in pairs) {
    baseline <- p[1]
    enhanced <- p[2]
    for (protocol in IMV_PROTOCOLS) {
      p0 <- as.vector(pred[[baseline]][[protocol]])
      p1 <- as.vector(pred[[enhanced]][[protocol]])
      val <- .imv_binary_compat(y, p0, p1)

      rows[[length(rows) + 1L]] <- data.frame(arm = meta$arm,
                                              dgp = meta$dgp,
                                              N = meta$N,
                                              J = meta$J,
                                              rep = meta$rep,
                                              seed = meta$seed,
                                              baseline = baseline,
                                              enhanced = enhanced,
                                              protocol = protocol,
                                              IMV = val,
                                              favours = if (is.na(val)) NA_character_ else if (val > 0) enhanced else if (val < 0) baseline else "Tie",
                                              baseline_hw_at_bound = .ipo_hw_at_bound(valid[[baseline]]),
                                              enhanced_hw_at_bound = .ipo_hw_at_bound(valid[[enhanced]]),
                                              regular_comparison = !.ipo_hw_at_bound(valid[[baseline]]) && !.ipo_hw_at_bound(valid[[enhanced]]),
                                              stringsAsFactors = FALSE)
    }
  }

  list(imv = if (length(rows)) do.call(rbind, rows) else .empty_df(),
       prediction = if (length(pred_rows)) do.call(rbind, pred_rows) else .empty_df())
}

.failed_rep_result <- function(arm, dgp, N, J, rep_id, seed, message) {
  health <- data.frame(arm = arm,
                       dgp = dgp,
                       N = as.integer(N),
                       J = as.integer(J),
                       rep = as.integer(rep_id),
                       seed = as.integer(seed),
                       model = "__replication__",
                       fit_ok = FALSE,
                       converged = FALSE,
                       logLik = NA_real_,
                       n_parameters = NA_integer_,
                       hw = NA_real_,
                       hw_at_bound = FALSE,
                       error = as.character(message),
                       stringsAsFactors = FALSE)
  list(health = health,
       parameter = .empty_df(),
       theta = .empty_df(),
       ic = .empty_df(),
       lr = .empty_df(),
       imv = .empty_df(),
       prediction = .empty_df())
}

### SECTION 8 - One-replication engines ========================================
.simulate_one_ipo_rep <- function(N, J, dgp, rep_id, condition_id, do_imv = TRUE) {
  seed <- BASE_SEED + condition_id * 100000L + rep_id
  meta <- list(arm = "IPO_recovery",
               dgp = paste0("IPO-", dgp),
               N = as.integer(N),
               J = as.integer(J),
               rep = as.integer(rep_id),
               seed = as.integer(seed))

  sim <- .simulate_ipo_pair(N, N_TEST, J, dgp, seed)
  fits <- .fit_models(sim$train, MODELS_MAIN)

  pred <- if (isTRUE(do_imv)) {
    .pairwise_imv_tables(fits, sim$test, meta)
  } else {
    list(imv = .empty_df(), prediction = .empty_df())
  }

  list(health = .fit_health_table(fits, meta),
       parameter = .parameter_recovery_table(fits, sim, meta),
       theta = .theta_recovery_table(fits, sim$test, sim$theta_test, meta),
       ic = .ic_table(fits, nobs = N, meta = meta),
       lr = .lr_u_k_table(fits, meta),
       imv = pred$imv,
       prediction = pred$prediction)
}

.simulate_one_logistic_rep <- function(dgp, rep_id, condition_id, N, J, do_imv = TRUE) {
  seed <- BASE_SEED + 500000000L + condition_id * 100000L + rep_id
  meta <- list(arm = "logistic_misspecification",
               dgp = dgp,
               N = as.integer(N),
               J = as.integer(J),
               rep = as.integer(rep_id),
               seed = as.integer(seed))

  sim <- .simulate_logistic_pair(N, N_TEST, J, dgp, seed)
  fits <- .fit_models(sim$train, MODELS_LOGISTIC_ARM)

  pred <- if (isTRUE(do_imv)) {
    .pairwise_imv_tables(fits, sim$test, meta)
  } else {
    list(imv = .empty_df(), prediction = .empty_df())
  }

  list(health = .fit_health_table(fits, meta),
       parameter = .empty_df(),
       theta = .theta_recovery_table(fits, sim$test, sim$theta_test, meta),
       ic = .ic_table(fits, nobs = N, meta = meta),
       lr = .lr_u_k_table(fits, meta),
       imv = pred$imv,
       prediction = pred$prediction)
}

### SECTION 9 - Parallel and result-combination helpers =========================
setup_future <- function(workers = NULL,
                         strategy = c("multisession", "multicore", "sequential")) {
  strategy <- match.arg(strategy)
  if (is.null(workers)) {
    workers <- max(1L, as.integer(future::availableCores()) - 1L)
  }

  options(future.globals.maxSize = 8 * 1024^3)

  if (strategy == "multisession") {
    future::plan(future::multisession, workers = workers)
  } else if (strategy == "multicore") {
    future::plan(future::multicore, workers = workers)
  } else {
    future::plan(future::sequential)
  }

  message("future plan: ", strategy, " | workers: ", workers)
  invisible(workers)
}

.bind_component <- function(rep_results, component) {
  xs <- lapply(rep_results, `[[`, component)
  xs <- xs[vapply(xs, function(x) is.data.frame(x) && nrow(x) > 0L, logical(1))]
  if (!length(xs)) return(.empty_df())
  do.call(rbind, xs)
}

.combine_rep_results <- function(rep_results) {
  list(health = .bind_component(rep_results, "health"),
       parameter = .bind_component(rep_results, "parameter"),
       theta = .bind_component(rep_results, "theta"),
       ic = .bind_component(rep_results, "ic"),
       lr = .bind_component(rep_results, "lr"),
       imv = .bind_component(rep_results, "imv"),
       prediction = .bind_component(rep_results, "prediction"))
}

.combine_condition_results <- function(condition_results) {
  components <- c("health", "parameter", "theta", "ic", "lr", "imv", "prediction")
  out <- setNames(vector("list", length(components)), components)
  for (component in components) {
    xs <- lapply(condition_results, `[[`, component)
    xs <- xs[vapply(xs, function(x) is.data.frame(x) && nrow(x) > 0L, logical(1))]
    out[[component]] <- if (length(xs)) do.call(rbind, xs) else .empty_df()
  }
  out
}

### SECTION 10 - Main IPO recovery / model-recovery study ======================
run_ipo_recovery_study <- function(
    n_reps = N_REPS,
    n_imv_reps = n_reps,
    overwrite = FALSE,
    results_dir = RESULTS_DIR) {

  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

  design <- expand.grid(N = N_VALUES, J = J_VALUES,
                        DGP = DGP_IPO,
                        stringsAsFactors = FALSE)
  design$condition_id <- seq_len(nrow(design))

  condition_results <- vector("list", nrow(design))

  for (cc in seq_len(nrow(design))) {
    d <- design[cc, ]
    file <- file.path(results_dir,
                      sprintf("IPO_recovery_DGP-%s_N-%d_J-%d.rds", d$DGP, d$N, d$J))

    if (file.exists(file) && !isTRUE(overwrite)) {
      message("Reading existing condition: ", basename(file))
      condition_results[[cc]] <- readRDS(file)
      next
    }

    message(sprintf("Running IPO recovery: DGP=%s, N=%d, J=%d, reps=%d",
                    d$DGP, d$N, d$J, n_reps))

    reps <- seq_len(n_reps)
    rr <- future.apply::future_lapply(
      reps,
      function(r) {
        tryCatch(.simulate_one_ipo_rep(N = d$N,
                                       J = d$J,
                                       dgp = d$DGP,
                                       rep_id = r,
                                       condition_id = d$condition_id,
                                       do_imv = r <= n_imv_reps),
                 error = function(e) {
                   seed <- BASE_SEED + d$condition_id * 100000L + r
                   .failed_rep_result(arm = "IPO_recovery", dgp = paste0("IPO-", d$DGP),
                                      N = d$N, J = d$J, rep_id = r, seed = seed,
                                      message = conditionMessage(e))})
      },
      future.seed = TRUE,
      future.packages = c("mirt", "imv", "statmod", "numDeriv")
    )

    condition_results[[cc]] <- .combine_rep_results(rr)
    saveRDS(condition_results[[cc]], file)
  }

  out <- .combine_condition_results(condition_results)
  attr(out, "design") <- design
  attr(out, "n_reps") <- n_reps
  attr(out, "n_imv_reps") <- n_imv_reps
  saveRDS(out, file.path(results_dir, "IPO_recovery_ALL.rds"))
  out
}

### SECTION 11 - Logistic-family misspecification arm ===========================
run_logistic_misspecification_study <- function(n_reps = N_REPS,
                                                n_imv_reps = n_reps,
                                                N = max(N_VALUES),
                                                J = max(J_VALUES),
                                                overwrite = FALSE,
                                                results_dir = RESULTS_DIR) {

  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  dgps <- c("Rasch", "2PL", "3PL")
  condition_results <- vector("list", length(dgps))

  for (cc in seq_along(dgps)) {
    dgp <- dgps[cc]
    file <- file.path(results_dir,
                      sprintf("logistic_truth_DGP-%s_N-%d_J-%d.rds", dgp, N, J))

    if (file.exists(file) && !isTRUE(overwrite)) {
      message("Reading existing condition: ", basename(file))
      condition_results[[cc]] <- readRDS(file)
      next
    }

    message(sprintf("Running logistic misspecification: DGP=%s, N=%d, J=%d, reps=%d",
                    dgp, N, J, n_reps))

    rr <- future.apply::future_lapply(seq_len(n_reps), function(r) {
        tryCatch(.simulate_one_logistic_rep(dgp = dgp,
                                            rep_id = r,
                                            condition_id = cc,
                                            N = N,
                                            J = J,
                                            do_imv = r <= n_imv_reps),
                 error = function(e) {
                   seed <- BASE_SEED + 500000000L + cc * 100000L + r
                   .failed_rep_result(arm = "logistic_misspecification", dgp = dgp,
                                      N = N, J = J, rep_id = r, seed = seed,
                                      message = conditionMessage(e))})
      },
      future.seed = TRUE,
      future.packages = c("mirt", "imv", "statmod", "numDeriv"))

    condition_results[[cc]] <- .combine_rep_results(rr)
    saveRDS(condition_results[[cc]], file)
  }

  out <- .combine_condition_results(condition_results)
  attr(out, "N") <- N
  attr(out, "J") <- J
  attr(out, "n_reps") <- n_reps
  attr(out, "n_imv_reps") <- n_imv_reps
  saveRDS(out, file.path(results_dir, "logistic_misspecification_ALL.rds"))
  out
}

### SECTION 12 - IMV prediction-protocol supplement ===================
run_imv_protocol_supplement <- function(n_reps = PROTOCOL_N_REPS,
                                        N_train = PROTOCOL_N_TRAIN,
                                        N_test = PROTOCOL_N_TEST,
                                        J_values = PROTOCOL_J_VALUES,
                                        overwrite = FALSE,
                                        results_dir = RESULTS_DIR) {

  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

  design <- expand.grid(truth = c("K", "2PL"),
                        J = J_values,
                        stringsAsFactors = FALSE)
  design$condition_id <- seq_len(nrow(design))
  all_rows <- vector("list", nrow(design))

  for (cc in seq_len(nrow(design))) {
    d <- design[cc, ]
    file <- file.path(results_dir,
                      sprintf("protocol_truth-%s_J-%d.rds", d$truth, d$J))

    if (file.exists(file) && !isTRUE(overwrite)) {
      message("Reading existing protocol condition: ", basename(file))
      all_rows[[cc]] <- readRDS(file)
      next
    }

    rr <- future.apply::future_lapply(seq_len(n_reps), function(r) {
      seed <- BASE_SEED + 800000000L + d$condition_id * 100000L + r
      meta <- list(arm = "IMV_protocol_supplement",
                   dgp = if (d$truth == "K") "IPO-K" else "2PL",
                   N = as.integer(N_train),
                   J = as.integer(d$J),
                   rep = as.integer(r),
                   seed = as.integer(seed))

        if (d$truth == "K") {
          sim <- .simulate_ipo_pair(N_train, N_test, d$J, "K", seed)
        } else {
          sim <- .simulate_logistic_pair(N_train, N_test, d$J, "2PL", seed)
        }

        fits <- .fit_models(sim$train, MODELS_MAIN)
        pr <- .pairwise_imv_tables(fits, sim$test, meta)
        list(health = .fit_health_table(fits, meta),
             parameter = .empty_df(),
             theta = .empty_df(),
             ic = .empty_df(),
             lr = .empty_df(),
             imv = pr$imv,
             prediction = pr$prediction)
      },
      future.seed = TRUE,
      future.packages = c("mirt", "imv", "statmod", "numDeriv")
    )

    all_rows[[cc]] <- .combine_rep_results(rr)
    saveRDS(all_rows[[cc]], file)
  }

  out <- .combine_condition_results(all_rows)
  attr(out, "design") <- design
  attr(out, "n_reps") <- n_reps
  saveRDS(out, file.path(results_dir, "IMV_protocol_supplement_ALL.rds"))
  out
}

### SECTION 13 - Simple summary functions ======================================
summarise_parameter_recovery <- function(x, correctly_specified_only = TRUE) {
  d <- x$parameter
  if (!nrow(d)) return(d)
  if (correctly_specified_only) d <- d[d$correctly_specified, , drop = FALSE]
  if (!nrow(d)) return(d)

  key <- interaction(d$dgp, d$N, d$J, d$model, d$parameter, drop = TRUE)
  out <- lapply(split(d, key), function(z) {
    data.frame(dgp = z$dgp[1],
               N = z$N[1],
               J = z$J[1],
               model = z$model[1],
               parameter = z$parameter[1],
               mean_bias = mean(z$bias, na.rm = TRUE),
               mean_rmse = mean(z$rmse, na.rm = TRUE),
               mean_correlation = mean(z$correlation, na.rm = TRUE),
               boundary_rate = mean(z$hw_at_bound, na.rm = TRUE),
               n_replications = nrow(z),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

summarise_theta_recovery <- function(x) {
  d <- x$theta
  if (!nrow(d)) return(d)
  key <- interaction(d$arm, d$dgp, d$N, d$J, d$model, drop = TRUE)
  out <- lapply(split(d, key), function(z) {
    data.frame(arm = z$arm[1],
               dgp = z$dgp[1],
               N = z$N[1],
               J = z$J[1],
               model = z$model[1],
               mean_theta_bias = mean(z$theta_bias, na.rm = TRUE),
               mean_theta_rmse = mean(z$theta_rmse, na.rm = TRUE),
               mean_theta_correlation = mean(z$theta_correlation, na.rm = TRUE),
               mean_posterior_sd = mean(z$mean_posterior_sd, na.rm = TRUE),
               boundary_rate = mean(z$hw_at_bound, na.rm = TRUE),
               n_replications = nrow(z),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

summarise_ic_selection <- function(x) {
  d <- x$ic
  if (!nrow(d)) return(d)
  key <- interaction(d$arm, d$dgp, d$N, d$J, d$criterion, d$model, drop = TRUE)
  out <- lapply(split(d, key), function(z) {
    data.frame(arm = z$arm[1],
               dgp = z$dgp[1],
               N = z$N[1],
               J = z$J[1],
               criterion = z$criterion[1],
               model = z$model[1],
               selection_rate = mean(z$selected, na.rm = TRUE),
               mean_deltaIC = mean(z$deltaIC, na.rm = TRUE),
               mean_weight = mean(z$model_weight, na.rm = TRUE),
               boundary_rate = mean(z$hw_at_bound, na.rm = TRUE),
               n_replications = nrow(z),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

summarise_imv <- function(x, regular_only = FALSE) {
  d <- x$imv
  if (!nrow(d)) return(d)
  if (regular_only) d <- d[d$regular_comparison, , drop = FALSE]
  if (!nrow(d)) return(d)

  key <- interaction(d$arm, d$dgp, d$N, d$J, d$baseline, d$enhanced,
                     d$protocol, drop = TRUE)
  out <- lapply(split(d, key), function(z) {
    data.frame(arm = z$arm[1],
               dgp = z$dgp[1],
               N = z$N[1],
               J = z$J[1],
               baseline = z$baseline[1],
               enhanced = z$enhanced[1],
               protocol = z$protocol[1],
               mean_IMV = mean(z$IMV, na.rm = TRUE),
               median_IMV = median(z$IMV, na.rm = TRUE),
               enhanced_win_rate = mean(z$IMV > 0, na.rm = TRUE),
               baseline_win_rate = mean(z$IMV < 0, na.rm = TRUE),
               irregular_rate = mean(!z$regular_comparison, na.rm = TRUE),
               n_replications = sum(is.finite(z$IMV)),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

### SECTION 14 - Convenience wrapper ===========================================
run_all_simulations <- function(workers = NULL, strategy = "multisession",
                                n_reps = N_REPS, n_imv_reps = n_reps,
                                protocol_n_reps = PROTOCOL_N_REPS,
                                run_protocol_supplement = TRUE,
                                overwrite = FALSE, results_dir = RESULTS_DIR) {

  setup_future(workers = workers, strategy = strategy)
  on.exit(future::plan(future::sequential), add = TRUE)

  main <- run_ipo_recovery_study(n_reps = n_reps, n_imv_reps = n_imv_reps,
                                 overwrite = overwrite, results_dir = results_dir)

  logistic <- run_logistic_misspecification_study(n_reps = n_reps,
                                                  n_imv_reps = n_imv_reps,
                                                  overwrite = overwrite,
                                                  results_dir = results_dir)

  protocol <- if (isTRUE(run_protocol_supplement)) {
    run_imv_protocol_supplement(n_reps = protocol_n_reps,
                                overwrite = overwrite,
                                results_dir = results_dir)
  } else {
    NULL
  }

  result <- list(main = main, logistic = logistic, protocol = protocol)
  saveRDS(result, file.path(results_dir, "IPO_IRT_SIMULATION_COMPLETE.rds"))
  result
}

### SECTION 15 - Publication figures ===========================================
.sim_plot_dir <- function(output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  normalizePath(output_dir, winslash = "/", mustWork = FALSE)
}

.sim_save_plot <- function(p, filename, output_dir, width, height, dpi = 600) {
  output_dir <- .sim_plot_dir(output_dir)
  ggplot2::ggsave(filename = file.path(output_dir, paste0(filename, ".png")),
                  plot = p, width = width, height = height, dpi = dpi, bg = "white")
  ggplot2::ggsave(filename = file.path(output_dir, paste0(filename, ".pdf")),
                  plot = p, width = width, height = height, device = "pdf")
  invisible(p)
}

.sim_theme <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "bottom",
                   strip.text = ggplot2::element_text(face = "bold"),
                   plot.title = ggplot2::element_text(face = "bold"),
                   plot.title.position = "plot")
}

.sim_long_metrics <- function(d, id_cols, metric_cols,
                              metric_name = "metric", value_name = "value") {
  out <- do.call(rbind, lapply(metric_cols, function(m) {
    z <- d[, id_cols, drop = FALSE]
    z[[metric_name]] <- m
    z[[value_name]] <- d[[m]]
    z
  }))
  rownames(out) <- NULL
  out
}

plot_parameter_recovery <- function(x,
                                    output_dir = file.path(RESULTS_DIR, "figures"),
                                    save = TRUE) {
  d <- summarise_parameter_recovery(x, correctly_specified_only = TRUE)
  if (!nrow(d)) stop("No parameter-recovery results were found.")
  
  # Figure 1: location and common-width recovery.
  d1 <- d[d$parameter %in% c("b", "hw"), , drop = FALSE]
  long1 <- .sim_long_metrics(d1,
                             id_cols = c("dgp", "N", "J", "model", "parameter"),
                             metric_cols = c("mean_bias", "mean_rmse"))
  long1$metric <- factor(long1$metric,
                         levels = c("mean_bias", "mean_rmse"),
                         labels = c("Bias", "RMSE"))
  long1$parameter <- factor(long1$parameter,
                            levels = c("b", "hw"),
                            labels = c("Item difficulty (b)", "Common half-width (h)"))
  long1$J <- factor(long1$J)
  
  p1 <- ggplot2::ggplot(long1,
                        ggplot2::aes(x = N, y = value, group = J,
                                     linetype = J, shape = J)) +
        ggplot2::geom_hline(data = subset(long1, metric == "Bias"),
                            yintercept = 0, linewidth = 0.35) +
        ggplot2::geom_line(linewidth = 0.65) +
        ggplot2::geom_point(size = 2.0) +
        ggplot2::facet_grid(rows = ggplot2::vars(parameter, metric),
                            cols = ggplot2::vars(dgp),
                            scales = "free_y") +
        ggplot2::scale_x_continuous(breaks = N_VALUES) +
        ggplot2::labs(title = "Parameter recovery under correctly specified IPO-IRT models",
                      x = "Training sample size (N)",
                      y = NULL,
                      linetype = "Items (J)",
                      shape = "Items (J)") +
        .sim_theme()
  
  # Figure 2: relative-sharpness recovery for IPO-K.
  d2 <- d[d$dgp == "IPO-K" &
            d$parameter %in% c("k", "u_log_relative_k"),,drop = FALSE]
  long2 <- .sim_long_metrics(d2, id_cols = c("dgp", "N", "J", "model", "parameter"),
                             metric_cols = c("mean_rmse", "mean_correlation"))
  long2$metric <- factor(long2$metric,
                         levels = c("mean_rmse", "mean_correlation"),
                         labels = c("RMSE", "Correlation with generating values"))
  long2$parameter <- factor(long2$parameter,
                            levels = c("k", "u_log_relative_k"),
                            labels = c("Sharpness (k)", "Relative log-sharpness (u)"))
  long2$J <- factor(long2$J)
  
  p2 <- ggplot2::ggplot(long2, ggplot2::aes(x = N, y = value, group = J,
                                            linetype = J, shape = J)) +
        ggplot2::geom_line(linewidth = 0.65) +
        ggplot2::geom_point(size = 2.0) +
        ggplot2::facet_grid(rows = ggplot2::vars(metric),
                            cols = ggplot2::vars(parameter),
                            scales = "free_y") +
        ggplot2::scale_x_continuous(breaks = N_VALUES) +
        ggplot2::labs(title = "Recovery of item-specific relative sharpness under IPO-K",
                      x = "Training sample size (N)",
                      y = NULL,
                      linetype = "Items (J)",
                      shape = "Items (J)") +
        .sim_theme()
  
  if (save) {
    .sim_save_plot(p1, "Figure_1_parameter_recovery_b_h",
                   output_dir, width = 9.0, height = 8.5)
    .sim_save_plot(p2, "Figure_2_relative_sharpness_recovery",
                   output_dir, width = 8.5, height = 6.5)
  }
  
  invisible(list(parameter = p1, sharpness = p2))
}

plot_theta_recovery <- function(x,
                                output_dir = file.path(RESULTS_DIR, "figures"),
                                save = TRUE) {
  d <- summarise_theta_recovery(x)
  if (!nrow(d)) stop("No theta-recovery results were found.")
  d$N <- factor(d$N, levels = N_VALUES)
  d$J <- as.numeric(d$J)
  
  long <- .sim_long_metrics(d,
                            id_cols = c("arm", "dgp", "N", "J", "model"),
                            metric_cols = c("mean_theta_rmse", "mean_theta_correlation"))
  long$metric <- factor(long$metric,
                        levels = c("mean_theta_rmse", "mean_theta_correlation"),
                        labels = c("Theta RMSE", "Theta correlation"))
  
  p <- ggplot2::ggplot(long, ggplot2::aes(x = J, y = value, group = model,
                                          linetype = model, shape = model)) +
       ggplot2::geom_line(linewidth = 0.55) +
       ggplot2::geom_point(size = 1.8) +
       ggplot2::facet_grid(rows = ggplot2::vars(metric, dgp),
                           cols = ggplot2::vars(N),
                           scales = "free_y") +
       ggplot2::scale_x_continuous(breaks = J_VALUES) +
       ggplot2::labs(title = "Recovery of latent ability on the independent test sample",
                     x = "Number of items (J)",
                     y = NULL,
                     linetype = "Fitted model",
                     shape = "Fitted model") +
       .sim_theme(base_size = 10)
  
  if (save) {
    .sim_save_plot(p, "Figure_3_theta_recovery",
                   output_dir, width = 12.5, height = 9.0)
  }
  invisible(p)
}

plot_ic_recovery <- function(x,
                             output_dir = file.path(RESULTS_DIR, "figures"),
                             save = TRUE) {
  d <- summarise_ic_selection(x)
  if (!nrow(d)) stop("No information-criterion results were found.")
  
  # Figure 4: probability that the generating IPO model is selected.
  d$generating_model <- ifelse(d$dgp == "IPO-U", "IPO-U", "IPO-K")
  gen <- d[d$model == d$generating_model, , drop = FALSE]
  gen$J <- factor(gen$J, levels = J_VALUES)
  gen$criterion <- factor(gen$criterion, levels = IC_CRITERIA)
  
  p_gen <- ggplot2::ggplot(gen, ggplot2::aes(x = N, y = selection_rate,
                                             linetype = criterion, shape = criterion,
                                             group = criterion)) +
           ggplot2::geom_hline(yintercept = 0.5, linewidth = 0.3, linetype = 3) +
           ggplot2::geom_line(linewidth = 0.65) +
           ggplot2::geom_point(size = 2.0) +
           ggplot2::facet_grid(rows = ggplot2::vars(dgp),
                               cols = ggplot2::vars(J)) +
           ggplot2::scale_x_continuous(breaks = N_VALUES) +
           ggplot2::scale_y_continuous(limits = c(0, 1)) +
           ggplot2::labs(title = "Information-criterion recovery of the generating IPO model",
                         x = "Training sample size (N)",
                         y = "Selection rate",
                         linetype = "Criterion",
                         shape = "Criterion") +
           .sim_theme()
  
  # Figure 5: under IPO-K truth, show which models absorb the evidence when K
  # is not selected. This makes the competition with IPO-U and the 2.5PLM visible.
  comp <- d[d$dgp == "IPO-K" &
            d$model %in% c("IPO-U", "IPO-K", "2.5PL"),,drop = FALSE]
  comp$J <- factor(comp$J, levels = J_VALUES)
  comp$criterion <- factor(comp$criterion, levels = IC_CRITERIA)
  
  p_comp <- ggplot2::ggplot(comp, ggplot2::aes(x = N, y = selection_rate,
                                               linetype = model, shape = model,
                                               group = model)) +
            ggplot2::geom_line(linewidth = 0.6) +
            ggplot2::geom_point(size = 1.8) +
            ggplot2::facet_grid(rows = ggplot2::vars(criterion),
                                cols = ggplot2::vars(J)) +
            ggplot2::scale_x_continuous(breaks = N_VALUES) +
            ggplot2::scale_y_continuous(limits = c(0, 1)) +
            ggplot2::labs(title = "Model-selection competition under IPO-K truth",
                          x = "Training sample size (N)",
                          y = "Selection rate",
                          linetype = "Selected model",
                          shape = "Selected model") +
            .sim_theme(base_size = 10)
  
  if (save) {
    .sim_save_plot(p_gen, "Figure_4_generating_model_IC_recovery",
                   output_dir, width = 10.5, height = 6.5)
    .sim_save_plot(p_comp, "Figure_5_IPOK_IC_competition",
                   output_dir, width = 11.5, height = 7.5)
  }
  
  invisible(list(generating = p_gen, competition = p_comp))
}

plot_imv_main <- function(x, protocol = "pp",
                          output_dir = file.path(RESULTS_DIR, "figures"),
                          save = TRUE) {
  d <- summarise_imv(x, regular_only = FALSE)
  if (!nrow(d)) stop("No IMV results were found.")
  d <- d[d$protocol == protocol, , drop = FALSE]
  
  focal <- rbind(d[d$baseline == "IPO-U" & d$enhanced == "IPO-K", , drop = FALSE],
                 d[d$baseline == "IPO-U" & d$enhanced == "2.5PL", , drop = FALSE],
                 d[d$baseline == "IPO-K" & d$enhanced == "2.5PL", , drop = FALSE],
                 d[d$baseline == "IPO-K" & d$enhanced == "2PL", , drop = FALSE],
                 d[d$baseline == "IPO-K" & d$enhanced == "Rasch", , drop = FALSE])
  if (!nrow(focal)) stop("The requested focal IMV pairs were not found.")
  
  focal$pair <- paste(focal$baseline, "\u2192", focal$enhanced)
  focal$N <- factor(focal$N, levels = N_VALUES)
  
  p <- ggplot2::ggplot(focal, ggplot2::aes(x = J, y = mean_IMV,
                                           linetype = N, shape = N,
                                           group = N)) +
       ggplot2::geom_hline(yintercept = 0, linewidth = 0.35) +
       ggplot2::geom_line(linewidth = 0.55) +
       ggplot2::geom_point(size = 1.8) +
       ggplot2::facet_grid(rows = ggplot2::vars(dgp),
                           cols = ggplot2::vars(pair),
                           scales = "free_y") +
       ggplot2::scale_x_continuous(breaks = J_VALUES) +
       ggplot2::labs(title = paste0("Pairwise IMV under the ",
                                    ifelse(protocol == "pp", "posterior-predictive", "EAP plug-in"),
                                    " protocol"),
                     subtitle = "Positive values favor the model to the right of the arrow; negative values favor the baseline.",
                     x = "Number of items (J)",
                     y = "Mean IMV",
                     linetype = "Training N",
                     shape = "Training N") +
       .sim_theme(base_size = 10)
  
  if (save) {
    .sim_save_plot(p, paste0("Figure_6_pairwise_IMV_", protocol),
                   output_dir, width = 14.0, height = 6.5)
  }
  invisible(p)
}

plot_imv_protocol_sensitivity <- function(x,
                                          output_dir = file.path(RESULTS_DIR, "figures"),
                                          save = TRUE) {
  d <- summarise_imv(x, regular_only = FALSE)
  if (!nrow(d)) stop("No IMV results were found.")
  
  focal <- d[(d$baseline == "IPO-U" & d$enhanced == "IPO-K") |
             (d$baseline == "IPO-K" & d$enhanced == "2.5PL"),,drop = FALSE]
  focal$pair <- paste(focal$baseline, "\u2192", focal$enhanced)
  focal$N <- factor(focal$N, levels = sort(unique(focal$N)))
  
  p <- ggplot2::ggplot(focal, ggplot2::aes(x = J, y = mean_IMV,
                                           linetype = protocol, shape = protocol,
                                           group = protocol)) +
       ggplot2::geom_hline(yintercept = 0, linewidth = 0.35) +
       ggplot2::geom_line(linewidth = 0.6) +
       ggplot2::geom_point(size = 1.8) +
       ggplot2::facet_grid(rows = ggplot2::vars(dgp, pair),
                           cols = ggplot2::vars(N),
                           scales = "free_y") +
       ggplot2::labs(title = "Sensitivity of pairwise IMV to the prediction protocol",
                     x = "Number of items (J)",
                     y = "Mean IMV",
                     linetype = "Protocol",
                     shape = "Protocol") +
       .sim_theme(base_size = 9)
  
  if (save) {
    .sim_save_plot(p, "Supplement_Figure_IMV_protocol_sensitivity",
                   output_dir, width = 12.5, height = 8.0)
  }
  invisible(p)
}

plot_protocol_supplement <- function(x,
                                     output_dir = file.path(RESULTS_DIR, "figures"),
                                     save = TRUE) {
  d <- summarise_imv(x, regular_only = FALSE)
  if (!nrow(d)) stop("No protocol-supplement IMV results were found.")
  
  focal <- d[(d$baseline == "IPO-U" & d$enhanced == "IPO-K") |
             (d$baseline == "IPO-K" & d$enhanced == "Rasch"),,drop = FALSE]
  focal$pair <- paste(focal$baseline, "\u2192", focal$enhanced)
  
  p <- ggplot2::ggplot(focal, ggplot2::aes(x = J, y = mean_IMV,
                                           linetype = protocol, shape = protocol,
                                           group = protocol)) +
       ggplot2::geom_hline(yintercept = 0, linewidth = 0.35) +
       ggplot2::geom_line(linewidth = 0.65) +
       ggplot2::geom_point(size = 2) +
       ggplot2::facet_grid(rows = ggplot2::vars(dgp),
                           cols = ggplot2::vars(pair),
                           scales = "free_y") +
       ggplot2::labs(title = "Prediction-protocol convergence as test length increases",
                     x = "Number of items (J)",
                     y = "Mean IMV",
                     linetype = "Protocol",
                     shape = "Protocol") +
       .sim_theme()
  
  if (save) {
    .sim_save_plot(p, "Supplement_Figure_protocol_length_sweep",
                   output_dir, width = 9.5, height = 6.0)
  }
  invisible(p)
}

plot_logistic_misspecification <- function(x,
                                           output_dir = file.path(RESULTS_DIR, "figures"),
                                           save = TRUE) {
  if (is.null(x)) stop("The logistic misspecification result object is NULL.")
  
  icd <- summarise_ic_selection(x)
  thd <- summarise_theta_recovery(x)
  imd <- summarise_imv(x, regular_only = FALSE)
  if (!nrow(icd) || !nrow(thd) || !nrow(imd)) {
    stop("Incomplete logistic-arm results.")
  }
  
  icd$criterion <- factor(icd$criterion, levels = IC_CRITERIA)
  p_ic <- ggplot2::ggplot(icd, ggplot2::aes(x = model, y = selection_rate,
                                            shape = criterion, group = criterion)) +
          ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.5),
                              size = 2.2) +
          ggplot2::facet_grid(rows = ggplot2::vars(dgp),
                              cols = ggplot2::vars(criterion)) +
          ggplot2::scale_y_continuous(limits = c(0, 1)) +
          ggplot2::labs(title = "Model selection under logistic-family data-generating processes",
                        x = "Fitted model", y = "Selection rate", shape = "Criterion") +
          .sim_theme() +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1))
  
  long_theta <- .sim_long_metrics(thd, id_cols = c("dgp", "model"),
                                  metric_cols = c("mean_theta_rmse", "mean_theta_correlation"))
  long_theta$metric <- factor(long_theta$metric,
                              levels = c("mean_theta_rmse", "mean_theta_correlation"),
                              labels = c("Theta RMSE", "Theta correlation"))
  
  p_theta <- ggplot2::ggplot(long_theta, ggplot2::aes(x = model, y = value, shape = model)) +
             ggplot2::geom_point(size = 2.4) +
             ggplot2::facet_grid(rows = ggplot2::vars(metric),
                                 cols = ggplot2::vars(dgp),
                                 scales = "free_y") +
             ggplot2::labs(title = "Ability recovery under logistic-family truth",
                           x = "Fitted model", y = NULL, shape = "Model") +
             .sim_theme() +
             ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1))
  
  imd <- imd[imd$protocol == "pp", , drop = FALSE]
  imd$pair <- paste(imd$baseline, "\u2192", imd$enhanced)
  p_imv <- ggplot2::ggplot(imd, ggplot2::aes(x = pair, y = mean_IMV, shape = pair)) +
           ggplot2::geom_hline(yintercept = 0, linewidth = 0.35) +
           ggplot2::geom_point(size = 2.2) +
           ggplot2::facet_wrap(ggplot2::vars(dgp), scales = "free_y") +
           ggplot2::labs(title = "Posterior-predictive pairwise IMV under logistic-family truth",
                         subtitle = "Positive values favor the model to the right of the arrow.",
                         x = NULL, y = "Mean IMV", shape = "Pair") +
           .sim_theme(base_size = 9) +
           ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 55, hjust = 1),
                          legend.position = "none")
  
  if (save) {
    .sim_save_plot(p_ic, "Supplement_Figure_logistic_IC",
                   output_dir, width = 11.0, height = 7.0)
    .sim_save_plot(p_theta, "Supplement_Figure_logistic_theta",
                   output_dir, width = 9.5, height = 6.0)
    .sim_save_plot(p_imv, "Supplement_Figure_logistic_IMV",
                   output_dir, width = 12.0, height = 6.5)
  }
  
  invisible(list(ic = p_ic, theta = p_theta, imv = p_imv))
}

generate_all_simulation_figures <- function(results,
                                            output_dir = file.path(RESULTS_DIR, "figures"),
                                            include_logistic = TRUE,
                                            include_protocol = TRUE) {
  
  figs <- list()
  figs$parameter <- plot_parameter_recovery(results$main, output_dir = output_dir,
                                            save = TRUE)
  figs$theta <- plot_theta_recovery(results$main, output_dir = output_dir,
                                    save = TRUE)
  figs$ic <- plot_ic_recovery(results$main, output_dir = output_dir, save = TRUE)
  figs$imv <- plot_imv_main(results$main, protocol = "pp",
                            output_dir = output_dir, save = TRUE)
  figs$imv_sensitivity <- plot_imv_protocol_sensitivity(results$main,
                                                        output_dir = output_dir,
                                                        save = TRUE)
  
  if (isTRUE(include_protocol) && !is.null(results$protocol)) {
    figs$protocol <- plot_protocol_supplement(results$protocol,
                                              output_dir = output_dir,
                                              save = TRUE)
  }
  
  if (isTRUE(include_logistic) && !is.null(results$logistic)) {
    figs$logistic <- plot_logistic_misspecification(results$logistic,
                                                    output_dir = output_dir,
                                                    save = TRUE)
  }
  
  invisible(figs)
}

################################################################################
### Full main simulation (posterior-predictive + EAP-plugin IMV)================
results <- run_all_simulations(workers = max(1, future::availableCores() - 1),
                               strategy = "multisession", n_reps = 1000,
                               n_imv_reps = 1000, run_protocol_supplement = TRUE,
                               overwrite = FALSE)

### Summaries===================================================================
## Tables (Main)
summarise_parameter_recovery(results$main)
summarise_theta_recovery(results$main)
summarise_ic_selection(results$main)
summarise_imv(results$main, regular_only = FALSE)

## Figures
figures <- generate_all_simulation_figures(results)

## Tables (Logistic)
summarise_theta_recovery(results$logistic)
summarise_ic_selection(results$logistic)
summarise_imv(results$logistic, regular_only = FALSE)

## IMV Sensitivity
summarise_imv(results$protocol, regular_only = FALSE)

################################################################################
