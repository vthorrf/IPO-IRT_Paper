### ICAR (SAPA) cross-validated model comparison, with two prediction protocols
###
### Fits Rasch, 2PL and 2.5PL (mirt) and IPO-U and IPO-K (Utilities file) to
### Letter-Number and Matrix Reasoning, both six-option, with 5-fold CV by
### person, and computes pairwise out-of-sample IMV two ways:
###
###   "eap": EAP score from the person's other observed items, plugged into
###          the ICC (what imv_compare_irt does for both IPO and mirt fits).
###   "pp" : the ICC averaged over the person's posterior instead. Same
###          posterior and quadrature, just no plug-in.
###
### With 2 or 3 responses per person the posterior for theta is wide, and
### plugging in the EAP shifts the predicted probability, more so for models
### with a guessing floor. For IPO-K vs Rasch the two ways disagree; see the
### results file.
###
### Needs: IPO-IRT Utilities.R, mirt, imv, irw, dplyr, tidyr. The first irw 
### fetch asks you to authenticate with Redivis.
### Runtime: roughly 5 minutes per subtest with the folds run in parallel.

library(mirt)
library(dplyr)
library(tidyr)
library(parallel)  
library(future.apply)
source("IPO-IRT Utilities.R")

## ---- data -------------------------------------------------------------------

get_icar <- function(cache = "icar_sapa_cache.rds") {
  if (file.exists(cache)) return(readRDS(cache))
  d <- irw::irw_fetch("icar_sapa")
  saveRDS(d, cache)
  d
}

subtest_matrix <- function(d, prefix) {
  W <- d |>
    filter(grepl(prefix, item)) |>
    pivot_wider(id_cols = id, names_from = item, values_from = resp) |>
    select(-id) |>
    as.matrix()
  W <- W[rowSums(!is.na(W)) > 0, , drop = FALSE]
  storage.mode(W) <- "double"
  W
}

## ---- unified prediction machinery ------------------------------------------
## Each model is turned into an ICC on a common theta grid plus weights for
## its latent distribution, and predictions use the same code for all models,
## so a difference between protocols is not a difference in implementation.

theta_grid <- seq(-6, 6, length.out = 61)

model_on_grid <- function(fit, J) {
  if (inherits(fit, c("ipo_MML", "ipo_EM"))) {
    icc <- sapply(seq_len(J), function(j) {
      k_j <- if (length(fit$k) == 1L) fit$k else fit$k[j]
      ipo_p(theta_grid, fit$b[j], fit$hw, k_j, fit$alpha, fit$K)
    })
    w <- dnorm(theta_grid)                      # IPO fixes theta ~ N(0, 1)
  } else {
    icc <- sapply(seq_len(J), function(j) {
      tr <- mirt::probtrace(mirt::extract.item(fit, j), matrix(theta_grid))
      tr[, ncol(tr)]
    })
    gp <- mirt::coef(fit, simplify = TRUE)      # mirt Rasch estimates the
    m  <- gp$means[1]                           # latent variance, so use the
    s  <- sqrt(gp$cov[1, 1])                    # model's own distribution
    w <- dnorm(theta_grid, m, s)
  }
  list(icc = icc, w = w / sum(w))
}

## Leave-one-item-out predictions for every cell of Y, under both protocols.
## For item j, the person's posterior is computed from their other observed
## items only, by subtracting item j's contribution from the total loglik.
loo_predict <- function(Y, mg) {
  icc <- pmin(pmax(mg$icc, 1e-10), 1 - 1e-10)
  lw  <- log(mg$w)
  O  <- (!is.na(Y)) * 1
  Y0 <- Y; Y0[is.na(Y0)] <- 0
  liq <- Y0 %*% t(log(icc)) + (O - Y0) %*% t(log(1 - icc))
  P_pp <- P_eap <- matrix(NA_real_, nrow(Y), ncol(Y))
  for (j in seq_len(ncol(Y))) {
    contrib <- outer(Y0[, j], log(icc[, j])) +
               outer(O[, j] - Y0[, j], log(1 - icc[, j]))
    lp <- liq - contrib + matrix(lw, nrow(Y), length(theta_grid), byrow = TRUE)
    lp <- lp - apply(lp, 1, max)
    post <- exp(lp); post <- post / rowSums(post)
    P_pp[, j]  <- post %*% icc[, j]
    theta_eap  <- post %*% theta_grid
    P_eap[, j] <- approx(theta_grid, icc[, j], xout = theta_eap, rule = 2)$y
  }
  list(pp = P_pp, eap = P_eap)
}

## ---- cross-validation -------------------------------------------------------

## imv.binary takes p0/p1 in older versions and p1/p2 on current CRAN (this
## caused the Test.R error earlier); the wrapper handles both.
imv_binary_compat <- function(y, p_baseline, p_enhanced, sigma = 1e-4) {
  if ("p2" %in% names(formals(imv::imv.binary))) {
    imv::imv.binary(y = y, p1 = p_baseline, p2 = p_enhanced, sigma = sigma)
  } else {
    imv::imv.binary(y = y, p0 = p_baseline, p1 = p_enhanced, sigma = sigma)
  }
}

model_pairs <- list(c("Rasch", "IPO-K"), c("2PL", "IPO-K"), c("2.5PL", "IPO-K"),
                    c("IPO-U", "IPO-K"), c("Rasch", "2PL"), c("2.5PL", "IPO-U"))

run_subtest <- function(W, K, label, n_folds = max(1L, parallel::detectCores() - 1L), seed = 7) {
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
      `IPO-U` = ipo_mml(tr, mode = "uniform",    K = K, se = FALSE, verbose = FALSE),
      `IPO-K` = ipo_mml(tr, mode = "relative_k", K = K, se = FALSE, verbose = FALSE)
    )
    preds <- lapply(fits, function(ft) loo_predict(te, model_on_grid(ft, J)))
    y <- as.vector(te); keep <- !is.na(y)
    out <- list()
    for (p in model_pairs) for (proto in c("eap", "pp")) {
      out[[paste(p[1], p[2], proto, sep = "|")]] <- imv_binary_compat(
        y[keep],
        as.vector(preds[[p[1]]][[proto]])[keep],
        as.vector(preds[[p[2]]][[proto]])[keep]
      )
    }
    c(unlist(out), hw_U = fits$`IPO-U`$hw, hw_K = fits$`IPO-K`$hw)
  }

  # Set up parallel plan for Windows (multisession)
  old_plan <- future::plan()          # snapshot current plan
  on.exit(future::plan(old_plan), add = TRUE)   # restore on exit, always
  future::plan(future::multisession, workers = max(1L, n_folds))
  
  M <- do.call(rbind, future_lapply(seq_len(n_folds), one_fold,
                                    future.seed = TRUE))

  cat(sprintf("\n===== %s (K = %d), %d-fold CV, N = %d =====\n",
              label, K, n_folds, nrow(W)))
  for (p in model_pairs) {
    ce <- M[, paste(p[1], p[2], "eap", sep = "|")]
    cp <- M[, paste(p[1], p[2], "pp",  sep = "|")]
    cat(sprintf("%-6s -> %-6s : EAP-plugin %+.5f (sd %.5f) | post-pred %+.5f (sd %.5f)\n",
                p[1], p[2], mean(ce), sd(ce), mean(cp), sd(cp)))
  }
  cat(sprintf("hw estimates: IPO-U mean %.3f | IPO-K mean %.3f\n",
              mean(M[, "hw_U"]), mean(M[, "hw_K"])))
  invisible(M)
}

## ---- run --------------------------------------------------------------------

icar <- get_icar()
res_LN <- run_subtest(subtest_matrix(icar, "^LN"), K = 6, label = "Letter-Number")
future::plan()

res_MR <- run_subtest(subtest_matrix(icar, "^MR"), K = 6, label = "Matrix Reasoning")
future::plan()
