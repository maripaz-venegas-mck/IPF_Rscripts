# ipf_reliability_helpers.R
# Shared helpers for IPF exploration and benchmark scripts.
#
# Roles (read once):
# - Census-style table = IPF seed (starting counts).
# - Survey-style targets = margins IPF must match (often from a newer or parallel source).
# - Benchmark table = a known later census (synthetic here) used only to score the fit.
#
# Ipfp expects margins that agree with a single joint table (same implied grand total).
# If you merge sources with different totals, rescale or harmonize before Ipfp.
#
# Weaker targets (e.g. only one-way margins) fit many joint tables — IPF still returns one
# solution guided by the seed, but cell uncertainty vs the real world can be large.
#
# Later: add a Region dimension to the array and extend target.list the same way;
#        metrics stay cellwise on as.vector(array).

# ---- METRIC_GUIDE (plain language) ----
# Print with: cat(METRIC_GUIDE, sep = "\n")
METRIC_GUIDE <- c(
  "MAE ............... Typical error in people per cell.",
  "RMSE .............. Like MAE, but a few very wrong cells count more.",
  "MaxAbs ............ Worst single cell — where IPF is most stretched.",
  "Cor_cells ......... Do high cells stay high and low stay low? (1 = same ranking pattern)",
  "TVD ............... Total variation distance on shares (0–1). How different the overall split is.",
  "RelativeMAE ....... MAE divided by mean cell size — handy when tables differ in scale or sparsity.",
  "JS_bits ........... Jensen–Shannon (log base 2). Symmetric; best for comparing scenarios, not a pass/fail cutoff.",
  "KL (optional) ..... Asymmetric; can explode with zeros — interpret with care."
)

# Plain language: how IPF targets relate to real published tables
CONSTRAINT_GUIDE <- c(
  "age_edu_plus_emp .. Age×Education cross-tab + Employment margins (typical when a publication gives an interaction table plus employment totals).",
  "main_effects ...... Only one-way totals (each variable alone). Common when agencies release separate marginal tables without full interactions. Many different joint tables match the same one-way margins — IPF still returns one solution, but it is driven more by the census seed; error vs a true joint table can be large.",
  "age_edu_emp_urban . Age×Education + Urban (dim 3) + Employment margins (4D: 2222, 3322, …)."
)

#' All metrics we report: compare a reference table (e.g. simulated truth) to the IPF fit.
compute_metrics <- function(reference, fitted, include_kl = FALSE) {
  reference <- as.array(reference)
  fitted <- as.array(fitted)
  stopifnot(length(reference) == length(fitted))

  diff <- as.vector(reference - fitted)
  ref_v <- as.vector(reference)
  fit_v <- as.vector(fitted)

  mae <- mean(abs(diff))
  rmse <- sqrt(mean(diff^2))
  max_abs <- max(abs(diff))

  cor_cells <- stats::cor(ref_v, fit_v, method = "pearson")
  if (is.na(cor_cells)) cor_cells <- NA_real_

  s_ref <- sum(ref_v)
  s_fit <- sum(fit_v)
  p_ref <- ref_v / s_ref
  p_fit <- fit_v / s_fit
  tvd <- 0.5 * sum(abs(p_ref - p_fit))

  mean_ref <- mean(ref_v)
  relative_mae <- if (mean_ref > 0) mae / mean_ref else NA_real_

  js_bits <- jensen_shannon_bits(p_ref, p_fit)

  out <- c(
    MAE = mae,
    RMSE = rmse,
    MaxAbs = max_abs,
    Cor_cells = cor_cells,
    TVD = tvd,
    RelativeMAE = relative_mae,
    JS_bits = js_bits
  )

  if (isTRUE(include_kl)) {
    eps <- 1e-12
    p_f <- pmax(p_fit, eps)
    p_r <- pmax(p_ref, eps)
    kl <- sum(p_r * (log(p_r) - log(p_f)))
    out <- c(out, KL = kl)
  }
  out
}

#' Jensen–Shannon divergence, base-2 log, in [0, 1] for probability vectors.
#' Use to rank scenarios; do not treat a single cutoff as "good" or "bad."
jensen_shannon_bits <- function(p, q, eps = 1e-12) {
  p <- as.numeric(p)
  q <- as.numeric(q)
  p <- p / sum(p)
  q <- q / sum(q)
  p <- (p + eps) / sum(p + eps)
  q <- (q + eps) / sum(q + eps)
  m <- 0.5 * (p + q)
  kl_pm <- sum(p * (log2(p) - log2(m)))
  kl_qm <- sum(q * (log2(q) - log2(m)))
  0.5 * (kl_pm + kl_qm)
}

# ---- Census seeds (synthetic, national) ----

#' Small national table: Age x Education x Employment (employed vs not).
generate_base_222 <- function() {
  array(
    c(
      300, 200,
      250, 150,
      400, 100,
      350, 150
    ),
    dim = c(2, 2, 2),
    dimnames = list(
      Age = c("Young", "Adult"),
      Education = c("Low", "High"),
      Employment = c("Employed", "NotEmployed")
    )
  )
}

#' Larger table 3 x 3 x 2 — same variables, more categories.
generate_base_332 <- function() {
  d1 <- 3L
  d2 <- 3L
  d3 <- 2L
  x <- array(0, dim = c(d1, d2, d3))
  for (i in seq_len(d1)) {
    for (j in seq_len(d2)) {
      for (k in seq_len(d3)) {
        x[i, j, k] <- 80 + 12 * i + 7 * j + 25 * k + ((i + j + k) %% 3) * 5
      }
    }
  }
  dimnames(x) <- list(
    Age = c("Young", "Mid", "Older"),
    Education = c("Low", "Mid", "High"),
    Employment = c("Employed", "NotEmployed")
  )
  x
}

#' Age x Education x Urban x Employment (2 levels each). Urban before Employment in dims 3–4.
generate_base_2222 <- function() {
  x <- array(50, dim = c(2, 2, 2, 2))
  for (i in 1:2) for (j in 1:2) for (u in 1:2) for (e in 1:2) {
    x[i, j, u, e] <- 40 + 15 * i + 10 * j + 8 * u + 20 * e
  }
  dimnames(x) <- list(
    Age = c("Young", "Adult"),
    Education = c("Low", "High"),
    Urban = c("Urban", "Rural"),
    Employment = c("Employed", "NotEmployed")
  )
  x
}

#' Age x Education x Urban x Employment — 3x3x2x2 (finer Age/Edu than 2222).
generate_base_3322 <- function() {
  d1 <- 3L
  d2 <- 3L
  d3 <- 2L
  d4 <- 2L
  x <- array(0, dim = c(d1, d2, d3, d4))
  for (i in seq_len(d1)) {
    for (j in seq_len(d2)) {
      for (u in seq_len(d3)) {
        for (e in seq_len(d4)) {
          x[i, j, u, e] <- 35 + 11 * i + 8 * j + 7 * u + 18 * e + ((i + j + u + e) %% 4L) * 4
        }
      }
    }
  }
  dimnames(x) <- list(
    Age = c("Young", "Mid", "Older"),
    Education = c("Low", "Mid", "High"),
    Urban = c("Urban", "Rural"),
    Employment = c("Employed", "NotEmployed")
  )
  x
}

# ---- Drift: shift employment share inside demographic slices (informality-style story) ----

#' Last dimension = Employment. Redistribute counts within each slice over all earlier dims
#' (story: employment/informality-style shifts inside demographic cells).
apply_drift_employment_last <- function(seed_table, drift_intensity) {
  x <- array(as.numeric(seed_table), dim = dim(seed_table), dimnames = dimnames(seed_table))
  d <- dim(x)
  nd <- length(d)
  if (nd < 2L) stop("Table needs at least Age/Education + Employment dimensions.")
  L <- d[nd]
  idx_grid <- do.call(expand.grid, lapply(d[seq_len(nd - 1L)], seq_len))
  for (r in seq_len(nrow(idx_grid))) {
    ei <- idx_grid[r, , drop = FALSE]
    args_get <- c(list(x), as.list(ei), list(seq_len(L)))
    v <- as.numeric(do.call(`[`, args_get))
    tot <- sum(v)
    if (tot <= 0) next
    p <- v / tot
    if (L == 2L) {
      p0 <- p[1]
      shift <- drift_intensity * stats::runif(1, -1, 1)
      p1 <- min(max(p0 + shift, 0.05), 0.95)
      nv <- tot * c(p1, 1 - p1)
    } else {
      delta <- drift_intensity * stats::runif(L, -1 / L, 1 / L)
      p2 <- pmax(p + delta, 0.01 / L)
      p2 <- p2 / sum(p2)
      nv <- tot * p2
    }
    for (ee in seq_len(L)) {
      args_set <- c(list(x), as.list(ei), list(ee), list(value = nv[ee]))
      x <- do.call(`[<-`, args_set)
    }
  }
  x
}

# ---- Targets from a full table (survey world) ----

#' Build Ipfp inputs from a reference table and a constraint pattern.
build_targets <- function(ref_table, constraint_set = c("age_edu_plus_emp", "main_effects", "age_edu_emp_urban")) {
  constraint_set <- match.arg(constraint_set)
  x <- ref_table
  d <- dim(x)
  nd <- length(d)

  if (constraint_set == "age_edu_plus_emp") {
    stopifnot(nd == 3L)
    return(list(
      target.list = list(c(1L, 2L), 3L),
      target.data = list(apply(x, c(1, 2), sum), apply(x, 3, sum))
    ))
  }

  # One-way margins only; see CONSTRAINT_GUIDE in this file for real-world implications.
  if (constraint_set == "main_effects") {
    tl <- lapply(seq_len(nd), identity)
    td <- lapply(seq_len(nd), function(i) apply(x, i, sum))
    return(list(target.list = tl, target.data = td))
  }

  if (constraint_set == "age_edu_emp_urban") {
    stopifnot(nd == 4L)
    return(list(
      target.list = list(c(1L, 2L), 3L, 4L),
      target.data = list(
        apply(x, c(1, 2), sum),
        apply(x, 3, sum),
        apply(x, 4, sum)
      )
    ))
  }

  stop("Unknown constraint_set")
}

#' Run IPF; returns fitted array (same dim as seed).
run_ipf <- function(seed, target.list, target.data, iter = 2000L, tol = 1e-10) {
  res <- mipfp::Ipfp(
    seed = seed,
    target.list = target.list,
    target.data = target.data,
    iter = iter,
    tol = tol
  )
  array(res$x.hat, dim = dim(seed), dimnames = dimnames(seed))
}

#' Check that fitted margins match targets (max absolute error on each constraint).
constraint_check <- function(fitted, target.list, target.data) {
  errs <- numeric(length(target.list))
  for (i in seq_along(target.list)) {
    mi <- apply(fitted, target.list[[i]], sum)
    errs[i] <- max(abs(mi - target.data[[i]]))
  }
  errs
}

# ---- Scenario grid (crossed factors, profile x constraint compatibility) ----

#' Number of dimensions for a synthetic census profile string.
profile_ndim <- function(profile) {
  switch(as.character(profile),
    "222" = 3L,
    "332" = 3L,
    "2222" = 4L,
    "3322" = 4L,
    stop("Unknown profile: ", profile)
  )
}

#' Constraint sets valid for Ipfp given table dimensionality.
#' age_edu_plus_emp: 3D only. age_edu_emp_urban: 4D only. main_effects: any.
constraints_for_profile <- function(profile) {
  nd <- profile_ndim(profile)
  if (nd == 3L) {
    c("age_edu_plus_emp", "main_effects")
  } else {
    c("main_effects", "age_edu_emp_urban")
  }
}

constraint_abbr <- function(constraint_set) {
  m <- c(
    age_edu_plus_emp = "aepe",
    main_effects = "main",
    age_edu_emp_urban = "aeu"
  )
  unname(m[as.character(constraint_set)])
}

#' Build a filtered full factorial over drift, profile, and constraint.
#'
#' Rows with incompatible (profile, constraint) are dropped. Targets are always
#' exact margins of the simulated true table (no sampling noise).
#'
#' @param drifts Numeric vector of drift intensities (passed to apply_drift_employment_last).
#' @param profiles Character vector: 222, 332, 2222, 3322, etc.
#' @param constraints If NULL, all constraints that appear for any profile are candidates;
#'   rows are still filtered by compatibility.
expand_ipf_scenarios <- function(
    drifts,
    profiles = c("222", "332", "2222", "3322"),
    constraints = NULL) {
  drifts <- sort(unique(as.numeric(drifts)))
  profiles <- unique(as.character(profiles))
  if (is.null(constraints)) {
    constraints <- unique(unlist(lapply(profiles, constraints_for_profile), use.names = FALSE))
  } else {
    constraints <- unique(as.character(constraints))
  }
  out <- expand.grid(
    drift = drifts,
    profile = profiles,
    constraint = constraints,
    stringsAsFactors = FALSE
  )
  ok <- mapply(function(p, cst) cst %in% constraints_for_profile(p), out$profile, out$constraint)
  out <- out[ok, , drop = FALSE]
  rownames(out) <- NULL
  if (nrow(out) == 0L) {
    stop("expand_ipf_scenarios: no valid (profile, constraint) pairs after filtering.")
  }
  drift_key <- sprintf("%.4f", out$drift)
  drift_key <- gsub("\\.", "p", drift_key, fixed = FALSE)
  abbr <- constraint_abbr(out$constraint)
  out$scenario_id <- paste(out$profile, paste0("d", drift_key), abbr, sep = "_")
  out$scenario_label <- paste0(out$profile, " drift=", out$drift, " ", out$constraint)
  out
}

#' Median and IQR for each metric, grouped by scenario_id (or custom columns).
summarize_scenarios <- function(results_df, metric_names, group_col = "scenario_id") {
  if (nrow(results_df) == 0L) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  stopifnot(group_col %in% names(results_df))
  # Character grouping avoids factor-level / == mismatches and empty subsets.
  gvec <- as.character(results_df[[group_col]])
  rows <- list()
  for (sid in unique(gvec)) {
    if (is.na(sid)) {
      next
    }
    sub <- results_df[gvec == sid, , drop = FALSE]
    if (nrow(sub) == 0L) {
      next
    }
    # One-row frame: assigning into data.frame() (0 rows) breaks `[[<-` ("replacement has 1 row, data has 0").
    one <- stats::setNames(
      data.frame(sid, stringsAsFactors = FALSE),
      group_col
    )
    for (nm in metric_names) {
      v <- sub[[nm]]
      qs <- stats::quantile(v, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
      one[[paste0(nm, "_Q1")]] <- unname(qs[1L])
      one[[paste0(nm, "_median")]] <- unname(qs[2L])
      one[[paste0(nm, "_Q3")]] <- unname(qs[3L])
    }
    rows[[length(rows) + 1L]] <- one
  }
  if (length(rows) == 0L) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}

#' Median of a metric by factor combinations (for heatmaps / line plots).
summarize_by_factors <- function(results_df, metric, group_cols) {
  stopifnot(metric %in% names(results_df))
  g <- interaction(results_df[group_cols], drop = TRUE, sep = " | ")
  spl <- split(seq_len(nrow(results_df)), g)
  out <- data.frame(stringsAsFactors = FALSE)
  med_col <- paste0(metric, "_median")
  for (nm in names(spl)) {
    idx <- spl[[nm]]
    v <- results_df[[metric]][idx]
    first <- idx[1L]
    row <- results_df[first, group_cols, drop = FALSE]
    row[[med_col]] <- stats::median(v, na.rm = TRUE)
    row$n_draws <- length(idx)
    out <- rbind(out, row)
  }
  rownames(out) <- NULL
  out
}
