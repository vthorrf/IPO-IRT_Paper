### ICAR (SAPA) empirical example for the IPO-IRT manuscript
###
### Full-sample fits of Rasch, 2PL, 2.5PL (g = 1/6), IPO-U and IPO-K on the
### Letter-Number and Matrix Reasoning subtests (logLik, AIC/BIC/HQ, U vs K
### likelihood ratio, h with SE, boundary checks), a profile log-likelihood for
### h over [0.05, 4], and 5-fold cross-validated pairwise IMV under the
### posterior-predictive and EAP plug-in protocols.
###
### IPO fits are multi-start: the estimate is the best of free fits from six
### starting widths, a free fit refined from the best pinned-profile point, and
### that pinned point itself. A single default start can end in a local mode
### (Matrix Reasoning, IPO-K). CV folds add the full-sample solution as a start.
###
### CV folds, seed and prediction code come from IPO-IRT ICAR CV Protocol.R (its
### definitions are evaluated below). IC definitions match the simulation code:
### nobs = persons, AIC = -2LL + 2p, BIC = -2LL + log(N) p, HQ = -2LL + 2 log(log N) p.
###
### The reported profile also starts each pinned fit from the final estimate.
### The saved results were produced in stages with the ICAR_REPROFILE and
### ICAR_RECV modes below. A full run does the same steps in one pass and
### reproduces them exactly.
###
### Needs icar_sapa_cache.rds (built by IPO-IRT ICAR CV Protocol.R). About 3 hours on
### 10 cores; results are saved after each subtest.

suppressMessages({
  library(mirt)
  library(dplyr)
  library(tidyr)
  library(parallel)
  library(future)
  library(future.apply)
})

src <- readLines("IPO-IRT ICAR CV Protocol.R")
cut <- grep("^## ---- run ---", src)[1]
stopifnot(!is.na(cut))
eval(parse(text = paste(src[seq_len(cut - 1)], collapse = "\n")), envir = globalenv())

RESULTS_FILE <- "icar_empirical_results.txt"
FITS_FILE    <- "icar_empirical_fits.rds"
K_ALT     <- 6
HW_GRID   <- c(0.05, 0.15, 0.30, 0.45, 0.60, 0.80, 1.10, 1.50, 2.00, 2.50, 3.00, 3.50, 4.00)
HW_STARTS <- c(0.05, 0.30, 0.60, 1.00, 2.00, 3.50)   # full-sample fits
CV_STARTS <- c(0.05, 0.40, 0.90, 2.50)               # per fold, both IPO models
HW_BOUNDS <- c(0.01, 4)                              # ipo_mml defaults
K_BOUNDS  <- c(0.05, 20)                             # ipo_mml defaults
U_BOUND   <- 4                                       # relative log-k bound
N_CORES   <- max(1, parallel::detectCores() - 1)

## Parallel setup
TOP_PID <- Sys.getpid()
future::plan(future::multisession, workers = N_CORES)
in_worker <- function() Sys.getpid() != TOP_PID
#' Like mclapply, but Windows-safe and nesting-aware: sequential when there is
#' nothing to parallelize, or when we are already running inside a future.
par_lapply <- function(X, FUN, ..., seed = TRUE) {
  if (length(X) <= 1L || in_worker()) {
    lapply(X, FUN, ...)
  } else {
    future.apply::future_lapply(X, FUN, ..., future.seed = seed)
  }
}

## ICAR_REPORT_ONLY=1 rebuilds the results file and figures from FITS_FILE.
REPORT_ONLY <- identical(Sys.getenv("ICAR_REPORT_ONLY"), "1")
## ICAR_REPROFILE=1 recomputes only the h profiles from the saved fits.
REPROFILE <- identical(Sys.getenv("ICAR_REPROFILE"), "1")
## ICAR_RECV=1 reruns only the cross-validation from the saved fits and saves
## the fold predictions, so pair orientations can be changed without refitting.
RECV <- identical(Sys.getenv("ICAR_RECV"), "1")
PREDS_FILE <- "icar_empirical_cv_preds.rds"

## Pairs in the orientation the simulation uses (combn over this model order),
## so the empirical IMVs read the same way as Figures 6 and 9.
MS_ORDER <- c("IPO-U", "IPO-K", "Rasch", "2PL", "2.5PL")
ms_pairs <- combn(MS_ORDER, 2, simplify = FALSE)

emit <- function(...) {
  cat(..., sep = "")
  cat(..., sep = "", file = RESULTS_FILE, append = TRUE)
}

## ---- multi-start IPO fitting ------------------------------------------------
## starts: starting widths (default b, k); extra: full parameter starts.
fit_ipo_ms <- function(W, mode, K, starts, extra = list(), se = FALSE, cores = N_CORES) {
  st <- c(lapply(starts, function(h0) list(hw = h0)), extra)
  fits <- par_lapply(st, function(s0) {
    ipo_mml(W, mode = mode, K = K, se = FALSE, verbose = FALSE, start = s0)
  })
  starts <- c(starts, rep(NA_real_, length(extra)))
  ll <- sapply(fits, `[[`, "loglik")
  best <- which.max(ll)
  f <- fits[[best]]
  if (se) {
    ## SE via the numerical Hessian, restarted from the best solution
    g <- ipo_mml(W, mode = mode, K = K, se = TRUE, verbose = FALSE,
                 start = list(b = f$b, hw = f$hw, k = f$k))
    if (g$loglik >= f$loglik - 1e-6) f <- g else f$se_refit_worse <- TRUE
  }
  f$multistart <- data.frame(
    start_hw = starts,
    hw = sapply(fits, `[[`, "hw"),
    loglik = ll,
    conv = sapply(fits, `[[`, "convergence"),
    best = seq_along(starts) == best
  )
  f
}

fit_full <- function(W, K) {
  list(Rasch   = mirt(data.frame(W), 1, itemtype = "Rasch", SE = FALSE, verbose = FALSE),
       `2PL`   = mirt(data.frame(W), 1, itemtype = "2PL",   SE = FALSE, verbose = FALSE),
       `2.5PL` = mirt(data.frame(W), 1, itemtype = "2PL", guess = 1 / K, SE = FALSE, verbose = FALSE),
       `IPO-U` = fit_ipo_ms(W, "uniform",    K, HW_STARTS),
       `IPO-K` = fit_ipo_ms(W, "relative_k", K, HW_STARTS))
}

## A refit started at the optimum can return L-BFGS-B code 52 with the
## likelihood unchanged; treat that as converged.
ipo_converged <- function(fit) {
  if (isTRUE(fit$converged)) return(TRUE)
  identical(as.integer(fit$convergence), 52L) && !is.null(fit$route_logliks) &&
    fit$loglik >= max(fit$route_logliks) - 1e-6
}

model_info <- function(fit) {
  if (inherits(fit, "ipo_MML")) {
    list(loglik = fit$loglik, npar = fit$n_par, converged = ipo_converged(fit))
  } else {
    list(loglik = as.numeric(extract.mirt(fit, "logLik")),
         npar = as.integer(extract.mirt(fit, "nestpars")),
         converged = isTRUE(fit@OptimInfo$converged))
  }
}

ic_table <- function(fits, N) {
  info <- lapply(fits, model_info)
  ll <- sapply(info, `[[`, "loglik"); p <- sapply(info, `[[`, "npar")
  out <- data.frame(model = names(fits), logLik = ll, npar = p,
                    converged = sapply(info, `[[`, "converged"),
                    AIC = -2 * ll + 2 * p,
                    BIC = -2 * ll + log(N) * p,
                    HQ  = -2 * ll + 2 * log(log(N)) * p,
                    stringsAsFactors = FALSE)
  for (crit in c("AIC", "BIC", "HQ")) {
    d <- out[[crit]] - min(out[[crit]])
    out[[paste0("d", crit)]] <- d
    out[[paste0("w", crit)]] <- exp(-0.5 * d) / sum(exp(-0.5 * d))
  }
  out
}

ipo_diagnostics <- function(fit) {
  list(hw = fit$hw,
       hw_se = if (isTRUE(fit$hessian_ok)) fit$se$hw else NA_real_,
       hessian_ok = isTRUE(fit$hessian_ok),
       hessian_condition = fit$hessian_condition,
       hw_at_bound = any(abs(fit$hw - HW_BOUNDS) < 1e-6),
       u_at_bound = if (fit$mode == "relative_k") any(abs(abs(fit$u) - U_BOUND) < 1e-6) else FALSE,
       k_at_bound = if (fit$mode == "relative_k") any(fit$k <= K_BOUNDS[1] * 1.01 | fit$k >= K_BOUNDS[2] * 0.99) else FALSE,
       k_range = range(fit$k), k_geo_mean = fit$k_geometric_mean,
       b_range = range(fit$b))
}

## Share of the overlap state across items and theta in [-3, 3]; near 1 when
## h sits on the upper bound.
overlap_share <- function(fit, theta = seq(-3, 3, by = 0.1)) {
  J <- length(fit$b)
  po <- sapply(seq_len(J), function(j) {
    k_j <- if (length(fit$k) == 1L) fit$k else fit$k[j]
    ipo_state_prob(theta, fit$b[j], fit$hw, k_j)$pi_overlap
  })
  c(mean_overlap = mean(po), share_above_0.99 = mean(po > 0.99), share_below_0.01 = mean(po < 0.01))
}

## Profile log-likelihood for h: refit with h pinned to a hair-width interval.
## seed: optional (b, k) start for every pinned fit, tried alongside the
## default start; the better is kept.
profile_h <- function(W, mode, K, grid = HW_GRID, seed = NULL) {
  jobs <- lapply(grid, function(h) list(h = h, start = list(hw = h)))
  if (!is.null(seed)) {
    jobs <- c(jobs, lapply(grid, function(h) list(h = h, start = list(b = seed$b, hw = h, k = seed$k))))
  }
  all <- par_lapply(jobs, function(j) {
    ipo_mml(W, mode = mode, K = K, se = FALSE, verbose = FALSE,
            hw_bounds = c(j$h, j$h + 1e-8), start = j$start)
  })
  fits <- lapply(seq_along(grid), function(i) {
    cand <- all[c(i, if (!is.null(seed)) i + length(grid))]
    cand[[which.max(sapply(cand, `[[`, "loglik"))]]
  })
  P <- cbind(h = grid, loglik = sapply(fits, `[[`, "loglik"),
             conv = sapply(fits, `[[`, "convergence"))
  list(table = P, fits = fits)
}

## Best of: multi-start free fit, free fit refined from the best pinned point,
## the pinned point. SEs from a refit at the winner (dropped if it loses
## likelihood).
finalize_ipo <- function(W, mode, K, ms, prof) {
  pb  <- prof$fits[[which.max(prof$table[, "loglik"])]]
  ref <- ipo_mml(W, mode = mode, K = K, se = FALSE, verbose = FALSE,
                 start = list(b = pb$b, hw = pb$hw, k = pb$k))
  cand <- list(multistart = ms, profile_refined = ref, profile_pinned = pb)
  ll <- sapply(cand, `[[`, "loglik")
  f <- cand[[which.max(ll)]]
  g <- ipo_mml(W, mode = mode, K = K, se = TRUE, verbose = FALSE,
               start = list(b = f$b, hw = f$hw, k = f$k))
  if (g$loglik >= f$loglik - 1e-6) f <- g else f$se_refit_worse <- TRUE
  f$multistart <- ms$multistart
  f$route <- names(cand)[which.max(ll)]
  f$route_logliks <- ll
  f
}

## ---- cross-validation with multi-start IPO fits ------------------------------
## As run_subtest in icar_cv_protocol.R, with multi-start IPO fits.
run_subtest_ms <- function(W, K, full_fits, n_folds = 5, seed = 7) {
  extra_U <- list(list(b = full_fits[["IPO-U"]]$b, hw = full_fits[["IPO-U"]]$hw, k = full_fits[["IPO-U"]]$k))
  extra_K <- list(list(b = full_fits[["IPO-K"]]$b, hw = full_fits[["IPO-K"]]$hw, k = full_fits[["IPO-K"]]$k))
  set.seed(seed)
  fold <- sample(rep(seq_len(n_folds), length.out = nrow(W)))
  J <- ncol(W)
  one_fold <- function(f) {
    tr <- W[fold != f, ]
    te <- W[fold == f, ]
    fits <- list(
      Rasch   = mirt(data.frame(tr), 1, itemtype = "Rasch", SE = FALSE, verbose = FALSE),
      `2PL`   = mirt(data.frame(tr), 1, itemtype = "2PL",   SE = FALSE, verbose = FALSE),
      `2.5PL` = mirt(data.frame(tr), 1, itemtype = "2PL", guess = 1 / K, SE = FALSE, verbose = FALSE),
      `IPO-U` = fit_ipo_ms(tr, "uniform",    K, CV_STARTS, extra = extra_U, cores = 2),
      `IPO-K` = fit_ipo_ms(tr, "relative_k", K, CV_STARTS, extra = extra_K, cores = 2)
    )
    preds <- lapply(fits, function(ft) loo_predict(te, model_on_grid(ft, J)))
    y <- as.vector(te); keep <- !is.na(y)
    pv <- lapply(preds, function(pr) list(eap = as.vector(pr$eap)[keep], pp = as.vector(pr$pp)[keep]))
    list(
      row = c(cv_imv_row(y[keep], pv),
              hw_U = fits$`IPO-U`$hw, hw_K = fits$`IPO-K`$hw,
              hw_U_bound = any(abs(fits$`IPO-U`$hw - HW_BOUNDS) < 1e-6),
              hw_K_bound = any(abs(fits$`IPO-K`$hw - HW_BOUNDS) < 1e-6)),
      y = y[keep], preds = pv)
  }
  folds <- par_lapply(seq_len(n_folds), one_fold)
  M <- do.call(rbind, lapply(folds, `[[`, "row"))
  attr(M, "preds") <- lapply(folds, function(fl) list(y = fl$y, preds = fl$preds))
  M
}

## IMV for the July pairs and for every pair in the simulation's orientation.
cv_imv_row <- function(y, pv) {
  out <- c()
  for (p in unique(c(model_pairs, ms_pairs))) for (proto in c("eap", "pp")) {
    out[paste(p[1], p[2], proto, sep = "|")] <- imv_binary_compat(y, pv[[p[1]]][[proto]], pv[[p[2]]][[proto]])
  }
  out
}

## July numbers for the same pair, from icar_cv_results.txt.
july_values <- function(label) {
  if (!file.exists("icar_cv_results.txt")) return(NULL)
  july <- readLines("icar_cv_results.txt")
  start <- grep(paste0("^===== ", label), july)[1]
  block <- july[start + 1:6]
  out <- list()
  for (i in seq_along(model_pairs)) {
    ## the last four numbers on the line are mean/sd for eap then pp
    nums <- as.numeric(regmatches(block[i], gregexpr("[-+]?[0-9]*\\.?[0-9]+", block[i]))[[1]])
    out[[pair_name_key(model_pairs[[i]])]] <- tail(nums, 4)
  }
  out
}
pair_name_key <- function(p) paste(p[1], p[2], sep = "|")

## ---- reporting --------------------------------------------------------------
report_subtest <- function(W, K, label, saved = NULL) {
  N <- nrow(W); J <- ncol(W)
  emit(sprintf("\n\n===== %s: N = %d persons, J = %d items, K = %d, %d observed responses =====\n",
               label, N, J, K, sum(!is.na(W))))
  emit(sprintf("responses per person: mean %.2f, median %d; item p-values %.3f to %.3f\n",
               mean(rowSums(!is.na(W))), as.integer(median(rowSums(!is.na(W)))),
               min(colMeans(W, na.rm = TRUE)), max(colMeans(W, na.rm = TRUE))))

  t0 <- Sys.time()
  if (is.null(saved)) {
    fits <- fit_full(W, K)
    emit(sprintf("full-sample fits (IPO models from %d starting widths each): %.1f min\n",
                 length(HW_STARTS), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  } else {
    fits <- saved$fits
    emit(sprintf("full-sample fits (IPO models from %d starting widths each): from the saved run\n",
                 length(HW_STARTS)))
  }

  emit("\n--- Multi-start table (full sample, free fits) ---\n")
  for (nm in c("IPO-U", "IPO-K")) {
    ms <- fits[[nm]]$multistart
    emit(sprintf("%s: ", nm))
    emit(paste(sprintf("start %.2f -> h %.3f, logLik %.1f%s%s", ms$start_hw, ms$hw, ms$loglik,
                       ifelse(ms$conv != 0, " (nonconv)", ""), ifelse(ms$best, " [best]", "")),
               collapse = " | "), "\n")
  }

  emit("\n--- Profile log-likelihood for h (other parameters re-estimated at each h) ---\n")
  prof <- list()
  for (nm in c("uniform", "relative_k")) {
    t1 <- Sys.time()
    lab0 <- if (nm == "uniform") "IPO-U" else "IPO-K"
    pr <- if (is.null(saved)) profile_h(W, nm, K) else
          if (REPROFILE) profile_h(W, nm, K, seed = fits[[lab0]]) else saved$profile[[nm]]
    prof[[nm]] <- pr
    lab <- if (nm == "uniform") "IPO-U" else "IPO-K"
    P <- pr$table
    emit(sprintf("%s (%s), values relative to the profile maximum %.1f at h = %.2f:\n  ",
                 lab, if (is.null(saved) || REPROFILE) sprintf("%.1f min", as.numeric(difftime(Sys.time(), t1, units = "mins"))) else "saved",
                 max(P[, "loglik"]), P[which.max(P[, "loglik"]), "h"]))
    emit(paste(sprintf("h=%.2f %.1f%s", P[, "h"], P[, "loglik"] - max(P[, "loglik"]),
                       ifelse(P[, "conv"] != 0, "*", "")), collapse = " | "), "\n")
  }
  emit("(* marks a non-converged pinned refit)\n")
  for (nm in c("uniform", "relative_k")) {
    lab <- if (nm == "uniform") "IPO-U" else "IPO-K"
    if (!is.null(fits[[lab]]$loglik)) {
      gap <- max(prof[[nm]]$table[, "loglik"]) - fits[[lab]]$loglik
      emit(sprintf("%s: pinned-grid maximum minus final free fit = %+.1f%s\n", lab, gap,
                   if (gap < -2) "  (grid never reaches the free optimum: pinned fits are in a worse (b, k) basin)" else ""))
    }
  }

  emit("\n--- Final IPO estimates: best of multi-start, profile-refined, profile-pinned ---\n")
  for (nm in c("IPO-U", "IPO-K")) {
    mode <- if (nm == "IPO-U") "uniform" else "relative_k"
    if (is.null(saved)) fits[[nm]] <- finalize_ipo(W, mode, K, fits[[nm]], prof[[mode]])
    rl <- fits[[nm]]$route_logliks
    emit(sprintf("%s: winner = %s (multistart %.1f, profile_refined %.1f, profile_pinned %.1f); h = %.3f%s\n",
                 nm, fits[[nm]]$route, rl["multistart"], rl["profile_refined"], rl["profile_pinned"],
                 fits[[nm]]$hw,
                 if (isTRUE(fits[[nm]]$se_refit_worse)) "  (SE refit lost likelihood; no SE)" else ""))
  }

  ## Reported profile: each pinned fit also starts from the final estimate. The
  ## profile above was only used to choose the final estimate.
  if (is.null(saved)) {
    emit("\n--- Profile log-likelihood for h, restarted from the final estimates (reported) ---\n")
    for (nm in c("uniform", "relative_k")) {
      lab <- if (nm == "uniform") "IPO-U" else "IPO-K"
      prof[[nm]] <- profile_h(W, nm, K, seed = fits[[lab]])
      P <- prof[[nm]]$table
      emit(sprintf("%s, values relative to the profile maximum %.1f at h = %.2f:\n  ", lab,
                   max(P[, "loglik"]), P[which.max(P[, "loglik"]), "h"]))
      emit(paste(sprintf("h=%.2f %.1f%s", P[, "h"], P[, "loglik"] - max(P[, "loglik"]),
                         ifelse(P[, "conv"] != 0, "*", "")), collapse = " | "), "\n")
    }
  }

  ic <- ic_table(fits, N)
  emit("\n--- Model fit and information criteria (full sample) ---\n")
  emit(sprintf("%-6s %12s %5s %5s %12s %12s %12s %8s %8s %8s %6s %6s %6s\n",
               "model", "logLik", "npar", "conv", "AIC", "BIC", "HQ",
               "dAIC", "dBIC", "dHQ", "wAIC", "wBIC", "wHQ"))
  for (i in seq_len(nrow(ic))) {
    emit(sprintf("%-6s %12.1f %5d %5s %12.1f %12.1f %12.1f %8.1f %8.1f %8.1f %6.3f %6.3f %6.3f\n",
                 ic$model[i], ic$logLik[i], ic$npar[i], ic$converged[i],
                 ic$AIC[i], ic$BIC[i], ic$HQ[i], ic$dAIC[i], ic$dBIC[i], ic$dHQ[i],
                 ic$wAIC[i], ic$wBIC[i], ic$wHQ[i]))
  }
  for (crit in c("AIC", "BIC", "HQ")) {
    emit(sprintf("selected by %s: %s\n", crit, ic$model[which.min(ic[[crit]])]))
  }

  u <- model_info(fits[["IPO-U"]]); k <- model_info(fits[["IPO-K"]])
  lr <- 2 * (k$loglik - u$loglik); df <- k$npar - u$npar
  emit(sprintf("\n--- IPO-U vs IPO-K likelihood ratio ---\nLR = %.2f on %d df, p = %.3g\n",
               lr, df, pchisq(lr, df, lower.tail = FALSE)))

  emit("\n--- IPO parameter diagnostics (full sample, best of the starts) ---\n")
  for (nm in c("IPO-U", "IPO-K")) {
    d <- ipo_diagnostics(fits[[nm]])
    emit(sprintf("%s: h = %.3f (SE %s; Hessian %s, condition %.2g); h at bound: %s; u at bound: %s; k at bound: %s\n",
                 nm, d$hw, if (is.na(d$hw_se)) "n/a" else sprintf("%.3f", d$hw_se),
                 if (d$hessian_ok) "ok" else "NOT ok", d$hessian_condition,
                 d$hw_at_bound, d$u_at_bound, d$k_at_bound))
    emit(sprintf("       b range %.2f to %.2f; k range %.3f to %.3f (geometric mean %.3f)\n",
                 d$b_range[1], d$b_range[2], d$k_range[1], d$k_range[2], d$k_geo_mean))
    emit(sprintf("       optimiser: code %d, %s; counted as converged: %s\n",
                 fits[[nm]]$convergence, fits[[nm]]$message, ipo_converged(fits[[nm]])))
    if (d$k_at_bound) {
      ## Refit with the k bound raised to 200 to see whether the bound matters.
      if (is.null(fits[[nm]]$k_relax) && (is.null(saved) || REPROFILE)) {
        f0 <- fits[[nm]]
        g0 <- ipo_mml(W, mode = "relative_k", K = K, se = FALSE, verbose = FALSE,
                      start = list(b = f0$b, hw = f0$hw, k = f0$k), k_bounds = c(K_BOUNDS[1], 200))
        fits[[nm]]$k_relax <- c(items = sum(f0$k >= K_BOUNDS[2] * 0.99), k_max = max(g0$k),
                                hw = g0$hw, gain = g0$loglik - f0$loglik)
      }
      kr <- fits[[nm]]$k_relax
      if (!is.null(kr)) emit(sprintf("       k bound check: %d item(s) at k = %g; with the bound at 200, max k = %.1f, h = %.3f, logLik gain %+.1f\n",
                                     kr["items"], K_BOUNDS[2], kr["k_max"], kr["hw"], kr["gain"]))
    }
    ov <- overlap_share(fits[[nm]])
    emit(sprintf("       overlap state over theta in [-3, 3] and all items: mean %.3f; share > .99: %.3f; share < .01: %.3f\n",
                 ov["mean_overlap"], ov["share_above_0.99"], ov["share_below_0.01"]))
  }

  emit(sprintf("\n--- 5-fold CV pairwise IMV, both protocols; IPO fits from %d starting widths plus the full-sample solution, per fold ---\n",
               length(CV_STARTS)))
  t2 <- Sys.time()
  if (is.null(saved) || RECV) {
    M <- run_subtest_ms(W, K, fits)
    emit(sprintf("(%.1f min)\n", as.numeric(difftime(Sys.time(), t2, units = "mins"))))
  } else {
    M <- saved$cv
    emit("(from the saved run)\n")
  }
  jv <- july_values(label)
  emit(sprintf("%-6s -> %-6s   %-33s %-33s %s\n", "", "", "EAP-plugin mean (sd) [folds>0]",
               "post-pred mean (sd) [folds>0]", "change from July (eap, pp)"))
  for (p in model_pairs) {
    ce <- M[, paste(p[1], p[2], "eap", sep = "|")]
    cp <- M[, paste(p[1], p[2], "pp",  sep = "|")]
    jj <- jv[[pair_name_key(p)]]
    emit(sprintf("%-6s -> %-6s : %+.5f (%.5f) [%d/5]         %+.5f (%.5f) [%d/5]         %+.5f, %+.5f\n",
                 p[1], p[2], mean(ce), sd(ce), sum(ce > 0), mean(cp), sd(cp), sum(cp > 0),
                 mean(ce) - jj[1], mean(cp) - jj[3]))
  }
  if (all(paste(sapply(ms_pairs, `[`, 1), sapply(ms_pairs, `[`, 2), "pp", sep = "|") %in% colnames(M))) {
    emit("\nAll ten pairs in the simulation's orientation (as in Figures 6 and 9):\n")
    for (p in ms_pairs) {
      ce <- M[, paste(p[1], p[2], "eap", sep = "|")]
      cp <- M[, paste(p[1], p[2], "pp",  sep = "|")]
      emit(sprintf("%-6s -> %-6s : EAP-plugin %+.5f (%.5f) [%d/5] | post-pred %+.5f (%.5f) [%d/5]\n",
                   p[1], p[2], mean(ce), sd(ce), sum(ce > 0), mean(cp), sd(cp), sum(cp > 0)))
    }
  }
  emit(sprintf("fold h: IPO-U %s | IPO-K %s\n",
               paste(sprintf("%.3f%s", M[, "hw_U"], ifelse(M[, "hw_U_bound"] == 1, "b", "")), collapse = " "),
               paste(sprintf("%.3f%s", M[, "hw_K"], ifelse(M[, "hw_K_bound"] == 1, "b", "")), collapse = " ")))
  emit("(b marks a fold whose h sits on an optimiser bound)\n")

  list(fits = fits, ic = ic, lr = c(LR = lr, df = df), profile = prof, cv = M)
}

## ---- run --------------------------------------------------------------------

cat("", file = RESULTS_FILE)
emit("ICAR (SAPA) empirical example: model comparison, multi-start IPO fits\n")
emit("Doria Xiao, September 2026. Produced by icar_empirical.R.\n")
emit("Data: IRW table icar_sapa (Condon & Revelle SAPA), Letter-Number and Matrix\n")
emit("Reasoning subtests, six alternatives each, all persons with at least one response.\n")
emit(sprintf("%s; mirt %s; imv %s\n", R.version.string,
             as.character(packageVersion("mirt")), as.character(packageVersion("imv"))))
emit(sprintf("IPO starts: full sample %s; CV folds %s. Profile grid: %s.\n",
             paste(HW_STARTS, collapse = ","), paste(CV_STARTS, collapse = ","),
             paste(HW_GRID, collapse = ",")))

icar <- get_icar()
res <- list()
saved_all <- if (REPORT_ONLY || REPROFILE || RECV) readRDS(FITS_FILE) else NULL
if (RECV) emit("Cross-validation rerun from the saved fits (ICAR_RECV=1); all else saved.\n")
if (REPORT_ONLY) emit("Report rebuilt from the saved fits (ICAR_REPORT_ONLY=1); no refitting.\n")
if (REPROFILE)   emit("Profiles recomputed from the saved final fits (ICAR_REPROFILE=1); all else saved.\n")
for (spec in list(c("^LN", "Letter-Number"), c("^MR", "Matrix Reasoning"))) {
  W <- subtest_matrix(icar, spec[1])
  res[[spec[2]]] <- report_subtest(W, K_ALT, spec[2], saved = saved_all[[spec[2]]])
  if (REPROFILE) { saved_all[[spec[2]]]$profile <- res[[spec[2]]]$profile; saved_all[[spec[2]]]$fits <- res[[spec[2]]]$fits; saveRDS(saved_all, FITS_FILE) }
  if (RECV) {
    pr <- if (file.exists(PREDS_FILE)) readRDS(PREDS_FILE) else list()
    pr[[spec[2]]] <- attr(res[[spec[2]]]$cv, "preds"); saveRDS(pr, PREDS_FILE)
    attr(res[[spec[2]]]$cv, "preds") <- NULL
    saved_all[[spec[2]]]$cv <- res[[spec[2]]]$cv; saveRDS(saved_all, FITS_FILE)
  }
  if (!REPORT_ONLY && !REPROFILE && !RECV) {
    pr <- if (file.exists(PREDS_FILE)) readRDS(PREDS_FILE) else list()
    pr[[spec[2]]] <- attr(res[[spec[2]]]$cv, "preds"); saveRDS(pr, PREDS_FILE)
    attr(res[[spec[2]]]$cv, "preds") <- NULL
    saveRDS(res, FITS_FILE)
  }
}

## ---- manuscript outputs -----------------------------------------------------
## Table 1 (manuscript) and Table S5 (supplement), tab-separated.
TAB_FILE <- "icar_empirical_tables.tsv"
cat("", file = TAB_FILE)
tab <- function(...) cat(..., sep = "", file = TAB_FILE, append = TRUE)
tab("TABLE 1. Model fit on the ICAR subtests (full sample)\n")
tab("Subtest\tModel\tlogLik\tParameters\tdAIC\tdBIC\tdHQ\th\tSE(h)\tBoundary\n")
for (lab in names(res)) {
  ic <- res[[lab]]$ic
  for (i in seq_len(nrow(ic))) {
    m <- ic$model[i]; ft <- res[[lab]]$fits[[m]]
    isipo <- inherits(ft, "ipo_MML")
    d <- if (isipo) ipo_diagnostics(ft) else NULL
    flags <- if (isipo) paste(c(if (d$hw_at_bound) "h" else NULL, if (d$k_at_bound) "k" else NULL), collapse = ",") else ""
    tab(sprintf("%s\t%s\t%.1f\t%d\t%.1f\t%.1f\t%.1f\t%s\t%s\t%s\n", lab, m, ic$logLik[i], ic$npar[i],
                ic$dAIC[i], ic$dBIC[i], ic$dHQ[i],
                if (isipo) sprintf("%.6f", d$hw) else "", if (isipo && !is.na(d$hw_se) && !d$hw_at_bound) sprintf("%.6f", d$hw_se) else "",
                if (nchar(flags)) flags else if (isipo) "none" else ""))
  }
  tab(sprintf("%s\tIPO-U vs IPO-K\tLR = %.1f on %d df\n", lab, res[[lab]]$lr["LR"], res[[lab]]$lr["df"]))
}
tab("\nTABLE S5. Five-fold cross-validated pairwise IMV on the ICAR subtests, both prediction protocols\n")
tab("Subtest\tBaseline\tAlternative\tEAP-plugin mean\tEAP-plugin SD\tEAP folds favoring alternative\tPosterior-predictive mean\tPosterior-predictive SD\tPP folds favoring alternative\n")
for (lab in names(res)) {
  M <- res[[lab]]$cv
  for (p in ms_pairs) {
    ce <- M[, paste(p[1], p[2], "eap", sep = "|")]; cp <- M[, paste(p[1], p[2], "pp", sep = "|")]
    tab(sprintf("%s\t%s\t%s\t%+.4f\t%.4f\t%d/5\t%+.4f\t%.4f\t%d/5\n", lab, p[1], p[2],
                mean(ce), sd(ce), sum(ce > 0), mean(cp), sd(cp), sum(cp > 0)))
  }
}
emit(sprintf("\nManuscript tables written to %s.\n", TAB_FILE))

## Figures: ggplot2, theme_minimal, black and white. No titles inside the
## images; number, title and note are set in the manuscript text (APA).
suppressMessages({ library(ggplot2); library(ragg); library(scales) })
grDevices::pdf(NULL)   # ggplot measures text on the current device; without one, R writes Rplots.pdf
ARROW <- " \u2192 "   # right arrow, written as an escape to keep this file ASCII
fig_theme <- theme_minimal(base_size = 10) +
  theme(strip.text = element_text(face = "bold", size = 10),
        panel.grid.minor = element_blank(),
        legend.position = "bottom", legend.title = element_blank())
save_fig <- function(p, stem, h) {
  ggsave(paste0(stem, ".png"), p, width = 6.5, height = h, dpi = 300, device = ragg::agg_png, bg = "white")
  ## base pdf() cannot draw the arrow; use quartz on macOS and cairo elsewhere
  dev <- if (Sys.info()[["sysname"]] == "Darwin") {
    function(filename, width, height, ...) grDevices::quartz(type = "pdf", file = filename, width = width, height = height)
  } else grDevices::cairo_pdf
  ggsave(paste0(stem, ".pdf"), p, width = 6.5, height = h, device = dev)
}

## ---- pairwise IMV (horizontal dot plot) ------------------------------------------
imv_df <- do.call(rbind, lapply(names(res), function(lab) {
  M <- res[[lab]]$cv
  do.call(rbind, lapply(ms_pairs, function(p) do.call(rbind, lapply(c("pp", "eap"), function(proto) {
    v <- M[, paste(p[1], p[2], proto, sep = "|")]
    data.frame(subtest = lab, pair = paste0(p[1], ARROW, p[2]),
               protocol = if (proto == "pp") "Posterior-predictive" else "EAP plug-in",
               mean = mean(v), sd = sd(v))
  }))))
}))
imv_df$protocol <- factor(imv_df$protocol, levels = c("Posterior-predictive", "EAP plug-in"))
imv_df$pair <- factor(imv_df$pair, levels = rev(sort(unique(imv_df$pair))))   # alphabetical, top to bottom
dodge <- position_dodge(width = 0.6)
p_imv <- ggplot(imv_df, aes(x = mean, y = pair, shape = protocol)) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  geom_errorbarh(aes(xmin = mean - sd, xmax = mean + sd), height = 0.35, linewidth = 0.35, position = dodge) +
  geom_point(size = 2, stroke = 0.6, position = dodge) +
  scale_shape_manual(values = c(1, 2)) +
  facet_wrap(~ subtest) +
  labs(x = "Mean IMV (positive values favor the model to the right of the arrow)", y = NULL) +
  fig_theme + theme(panel.spacing.x = unit(1.2, "lines"))
save_fig(p_imv, "fig_icar_imv", 4.2)
emit("IMV figure written to fig_icar_imv.png/.pdf (means over 5 folds, bars are fold SDs).\n")

## ---- profile log-likelihood (one panel per model and subtest) ------------------------
prof_df <- do.call(rbind, lapply(names(res), function(lab) do.call(rbind, lapply(c("uniform", "relative_k"), function(nm) {
  P <- res[[lab]]$profile[[nm]]$table
  mod <- if (nm == "uniform") "IPO-U" else "IPO-K"
  data.frame(subtest = lab, model = mod, h = P[, "h"], y = P[, "loglik"] - res[[lab]]$fits[[mod]]$loglik)
}))))
est_df <- do.call(rbind, lapply(names(res), function(lab) data.frame(
  subtest = lab, model = c("IPO-U", "IPO-K"),
  h = c(res[[lab]]$fits[["IPO-U"]]$hw, res[[lab]]$fits[["IPO-K"]]$hw), y = 0)))
prof_df$model <- factor(prof_df$model, levels = c("IPO-U", "IPO-K"))
est_df$model <- factor(est_df$model, levels = c("IPO-U", "IPO-K"))
p_prof <- ggplot(prof_df, aes(h, y)) +
  geom_hline(yintercept = -1.92, linetype = "dotted", linewidth = 0.4) +
  geom_line(linewidth = 0.45) +
  geom_point(size = 1.5) +
  geom_point(data = est_df, shape = 8, size = 2.6, stroke = 0.7) +
  facet_grid(model ~ subtest) +
  scale_y_continuous(trans = pseudo_log_trans(sigma = 1, base = 10),
                     breaks = c(0, -3, -10, -30, -100, -300, -1000),
                     labels = c("0", "-3", "-10", "-30", "-100", "-300", "-1000")) +
  scale_x_continuous(breaks = 0:4) +
  labs(x = "Common half-width h", y = "Log-likelihood relative to the best fit") +
  fig_theme
save_fig(p_prof, "fig_icar_profile_h", 4.4)
emit("Profile figure written to fig_icar_profile_h.png/.pdf.\n")

## ---- cleanup ----------------------------------------------------------------
future::plan(future::sequential)
