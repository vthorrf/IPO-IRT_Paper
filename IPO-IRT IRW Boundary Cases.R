### Why the art and Shipley datasets are not IPO-IRT empirical examples
###
### Recomputes, from the IRW tables, what the Empirical Examples section says
### about the Author Recognition Test (no chance floor, since responses are
### yes/no recognition rather than a choice among K alternatives, and what a
### 1/4 floor does to the IPO-U fit) and about Shipley Vocabulary (hardest items
### below the 1/4 chance rate for low-scoring respondents). The first fetch asks
### for a Redivis login; both tables are cached after that.

suppressMessages({library(dplyr); library(tidyr)})
OUT <- "irw_boundary_cases.txt"
cat("", file = OUT)
emit <- function(...) { cat(..., sep = ""); cat(..., sep = "", file = OUT, append = TRUE) }

fetch_cached <- function(name) {
  f <- sprintf("%s_cache.rds", name)
  if (file.exists(f)) return(readRDS(f))
  d <- irw::irw_fetch(name)
  if (is.null(d) || !nrow(d)) stop("irw_fetch returned nothing for '", name, "'")
  saveRDS(d, f); d
}
to_matrix <- function(d) {
  W <- d |> select(id, item, resp) |> pivot_wider(names_from = item, values_from = resp) |> select(-id) |> as.matrix()
  storage.mode(W) <- "double"; W[rowSums(!is.na(W)) > 0, , drop = FALSE]
}
describe <- function(W, label, floor) {
  p <- colMeans(W, na.rm = TRUE)
  tot <- rowMeans(W, na.rm = TRUE)
  dec <- dplyr::ntile(tot, 10)   # handles tied totals, which quantile breaks do not
  hardest5 <- names(sort(p))[1:5]
  emit(sprintf("\n%s: N = %d, J = %d, observed cells %d; response values %s\n", label, nrow(W), ncol(W), sum(!is.na(W)),
               paste(sort(unique(as.vector(W))), collapse = "/")))
  emit(sprintf("  item p-values: min %.3f, max %.3f; items below the %.2f floor: %d of %d\n",
               min(p), max(p), floor, sum(p < floor), length(p)))
  emit(sprintf("  items with zero correct responses: %s\n",
               if (any(p == 0)) paste(names(p)[p == 0], collapse = ", ") else "none"))
  for (k in 1:2) emit(sprintf("  decile %d respondents, mean accuracy on the five hardest items (%s): %.3f\n", k,
                              paste(hardest5, collapse = ", "), mean(W[dec == k, hardest5], na.rm = TRUE)))
  hardest5_legit <- names(sort(p[p > 0]))[1:5]
  if (any(p == 0)) for (k in 1:2) emit(sprintf("  decile %d, five hardest items with any correct response (%s): %.3f\n", k,
                                               paste(hardest5_legit, collapse = ", "), mean(W[dec == k, hardest5_legit], na.rm = TRUE)))
  invisible(p)
}

emit(sprintf("Boundary cases for the IPO-IRT empirical section. %s\n", format(Sys.time(), "%Y-%m-%d")))
setTimeLimit(elapsed = 900, transient = TRUE)
art <- to_matrix(fetch_cached("art"))
describe(art, "art (Author Recognition Test, IRW table 'art')", floor = 0.25)

## What the floor does to the fit: IPO-U with the 1/4 floor an earlier
## application assumed, IPO-U with a floor near zero, and the Rasch model.
suppressMessages(library(mirt)); source("IPO-IRT Utilities.R")
fit_floor <- function(K) ipo_mml(art, mode = "uniform", K = K, se = FALSE, verbose = FALSE)
f4 <- fit_floor(4); f0 <- fit_floor(1000)
rasch <- mirt(data.frame(art), 1, itemtype = "Rasch", SE = FALSE, verbose = FALSE)
emit(sprintf("  logLik: IPO-U with floor 1/4 = %.1f (h = %.3f); IPO-U with floor 1/1000 = %.1f (h = %.3f); Rasch = %.1f\n",
             f4$loglik, f4$hw, f0$loglik, f0$hw, as.numeric(extract.mirt(rasch, "logLik"))))
emit(sprintf("  gap to Rasch: %.1f with the 1/4 floor, %.1f with the floor removed\n",
             as.numeric(extract.mirt(rasch, "logLik")) - f4$loglik, as.numeric(extract.mirt(rasch, "logLik")) - f0$loglik))
shp <- to_matrix(fetch_cached("himmelstein-shipley_vocabulary-2025"))
describe(shp, "Shipley Vocabulary (IRW table 'himmelstein-shipley_vocabulary-2025', four alternatives)", floor = 0.25)
