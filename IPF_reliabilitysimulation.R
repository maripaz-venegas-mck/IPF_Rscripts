# IPF_reliabilitysimulation.R
#
# Exploration-first script: census-style seed + survey-style targets + a simulated "real"
# population. Stress-tests IPF under crossed factors (drift x profile x constraint), not
# one-at-a-time scenarios. Targets are exact margins of the simulated truth.
#
# Workflow: run this script (Rscript or source), then explore saved results with
# shiny::runApp("ipf_reliability_app") after results are written to
# ipf_reliability_output/results_all.rds. Static PDFs are off by default (SAVE_PDF_PLOTS).
#
# Requires: mipfp, ggplot2
#
# Optional: reduce n_sim or drifts while developing (grid size grows quickly).

library(mipfp)
library(ggplot2)

# Load helpers from the same folder as this script (Rscript sets --file=...)
args_all <- commandArgs(trailingOnly = FALSE)
file_idx <- grep("^--file=", args_all)
if (length(file_idx) > 0L) {
  script_path <- sub("^--file=", "", args_all[file_idx[1L]])
  script_dir <- dirname(normalizePath(script_path))
  helper_path <- file.path(script_dir, "ipf_reliability_helpers.R")
} else {
  script_dir <- getwd()
  helper_path <- file.path(script_dir, "ipf_reliability_helpers.R")
}
if (!file.exists(helper_path)) {
  stop("Cannot find ipf_reliability_helpers.R. Tried: ", helper_path)
}
source(helper_path)

# ---- Run controls ----
n_sim <- 200L
set.seed(12231233)

SHOW_TOY_RUN <- TRUE
SHOW_SCATTER_ONE_REP <- TRUE

# Crossed grid knobs (trim for faster runs)
drift_levels <- c(0.05, 0.15, 0.30)
profile_levels <- c("222", "332", "2222", "3322")
# Legacy six-row table (reference only; set USE_LEGACY_SCENARIOS <- TRUE to run it instead)
USE_LEGACY_SCENARIOS <- FALSE
legacy_scenarios_tbl <- data.frame(
  scenario_name = c(
    "Gentle demographic shift",
    "Moderate demographic shift",
    "Strong demographic shift",
    "Finer age and education categories",
    "Weaker survey: only one-way totals",
    "Urban/rural dimension added"
  ),
  profile = c("222", "222", "222", "332", "222", "2222"),
  drift = c(0.05, 0.15, 0.30, 0.15, 0.15, 0.12),
  constraint = c(
    "age_edu_plus_emp", "age_edu_plus_emp", "age_edu_plus_emp",
    "age_edu_plus_emp", "main_effects", "age_edu_emp_urban"
  ),
  stringsAsFactors = FALSE
)

# ---- What each printed number tries to tell you ----
cat("\n--- METRIC_GUIDE (read once) ---\n")
cat(paste(METRIC_GUIDE, collapse = "\n"), "\n\n")
cat("\n--- CONSTRAINT_GUIDE (read once) ---\n")
cat(paste(CONSTRAINT_GUIDE, collapse = "\n"), "\n\n")

if (isTRUE(USE_LEGACY_SCENARIOS)) {
  scenarios_tbl <- legacy_scenarios_tbl
  scenarios_tbl$scenario_id <- scenarios_tbl$scenario_name
  scenarios_tbl$scenario_label <- scenarios_tbl$scenario_name
} else {
  scenarios_tbl <- expand_ipf_scenarios(
    drifts = drift_levels,
    profiles = profile_levels,
    constraints = NULL
  )
}

pick_seed_fn <- function(profile) {
  switch(as.character(profile),
    "222" = generate_base_222,
    "332" = generate_base_332,
    "2222" = generate_base_2222,
    "3322" = generate_base_3322,
    stop("Unknown profile: ", profile)
  )
}

metric_names <- c(
  "MAE", "RMSE", "MaxAbs", "Cor_cells", "TVD", "RelativeMAE", "JS_bits", "max_margin_error"
)

#' One Monte Carlo draw for a scenario row (single-row data.frame).
run_one_draw <- function(row) {
  seed_fn <- pick_seed_fn(row$profile)
  seed_table <- seed_fn()
  true_table <- apply_drift_employment_last(seed_table, row$drift)
  tg <- build_targets(true_table, constraint_set = row$constraint)

  fitted <- run_ipf(seed_table, tg$target.list, tg$target.data)
  chk <- constraint_check(fitted, tg$target.list, tg$target.data)

  m <- compute_metrics(true_table, fitted, include_kl = FALSE)
  out <- c(m, max_margin_error = max(chk))
  out[metric_names]
}

# ---- Toy run: one replication (moderate drift, 2x2x2) ----
if (isTRUE(SHOW_TOY_RUN)) {
  cat("\n========== TOY RUN (one draw, moderate shift, 2x2x2) ==========\n")
  if (isTRUE(USE_LEGACY_SCENARIOS)) {
    row0 <- scenarios_tbl[scenarios_tbl$scenario_name == "Moderate demographic shift", ][1L, , drop = FALSE]
  } else {
    row0 <- scenarios_tbl[
      scenarios_tbl$profile == "222" &
        abs(scenarios_tbl$drift - 0.15) < 1e-9 &
        scenarios_tbl$constraint == "age_edu_plus_emp",
    ][1L, , drop = FALSE]
  }
  seed_fn <- pick_seed_fn(row0$profile)
  seed_table <- seed_fn()
  true_table <- apply_drift_employment_last(seed_table, row0$drift)
  tg <- build_targets(true_table, constraint_set = row0$constraint)
  fitted <- run_ipf(seed_table, tg$target.list, tg$target.data)

  cat("Total people (census seed): ", sum(seed_table), "\n", sep = "")
  cat("Total people (simulated 'true' world): ", sum(true_table), "\n", sep = "")
  cat("Total people (after IPF): ", sum(fitted), "\n", sep = "")
  cat(
    "Largest mismatch on any fitted margin vs survey target: ",
    max(constraint_check(fitted, tg$target.list, tg$target.data)),
    " (should be near numerical tolerance)\n",
    sep = ""
  )
  cat("\nMetrics vs simulated truth for this single draw:\n")
  print(round(compute_metrics(true_table, fitted), 6))

  if (isTRUE(SHOW_SCATTER_ONE_REP)) {
    cat("\nPlot: each point is one cell — truth vs IPF (toy run).\n")
    sc <- data.frame(
      truth = as.vector(true_table),
      ipf = as.vector(fitted)
    )
  #   print(
  #     ggplot(sc, aes(truth, ipf)) +
  #       geom_point(alpha = 0.7, size = 2) +
  #       geom_abline(slope = 1, intercept = 0, linetype = 2) +
  #       labs(
  #         title = "One replication: cell counts, simulated truth vs IPF",
  #         subtitle = "Points on the dashed line would match exactly in every cell.",
  #         x = "True cell count (synthetic)",
  #         y = "IPF fitted cell count"
  #       ) +
  #       theme_minimal()
  #   )
   }
  cat("========== END TOY RUN ==========\n\n")
}

# ---- Monte Carlo by scenario ----
output_dir <- file.path(script_dir, "ipf_reliability_output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

results_list <- vector("list", nrow(scenarios_tbl))
for (s in seq_len(nrow(scenarios_tbl))) {
  row_s <- scenarios_tbl[s, , drop = FALSE]
  lab <- if ("scenario_label" %in% names(row_s)) row_s$scenario_label else row_s$scenario_name
  cat("Running scenario: ", lab, " ...\n", sep = "")
  mat <- t(vapply(seq_len(n_sim), function(...) run_one_draw(row_s), numeric(length(metric_names))))
  colnames(mat) <- metric_names
  df_s <- as.data.frame(mat)
  df_s$scenario_id <- row_s$scenario_id
  if ("scenario_label" %in% names(row_s)) {
    df_s$scenario_label <- row_s$scenario_label
  } else {
    df_s$scenario_label <- row_s$scenario_name
  }
  df_s$drift <- row_s$drift
  df_s$profile <- row_s$profile
  df_s$constraint <- row_s$constraint
  results_list[[s]] <- df_s
}

results_all <- do.call(rbind, results_list)
rownames(results_all) <- NULL

results_all$profile <- factor(
  results_all$profile,
  levels = c("222", "332", "2222", "3322")
)
results_all$constraint <- factor(
  results_all$constraint,
  levels = c("age_edu_plus_emp", "main_effects", "age_edu_emp_urban")
)
results_all$drift_label <- factor(
  sprintf("%g", results_all$drift),
  levels = sprintf("%g", sort(unique(as.numeric(as.character(results_all$drift)))))
)

saveRDS(results_all, file.path(output_dir, "results_all.rds"))
utils::write.csv(
  scenarios_tbl,
  file.path(output_dir, "scenario_factor_levels.csv"),
  row.names = FALSE
)

summary_by_scenario <- summarize_scenarios(results_all, metric_names, group_col = "scenario_id")
meta_cols <- intersect(
  c("scenario_id", "scenario_label", "drift", "profile", "constraint"),
  names(results_all)
)
scenario_meta <- unique(results_all[, meta_cols, drop = FALSE])
summary_by_scenario <- merge(scenario_meta, summary_by_scenario, by = "scenario_id", all.y = TRUE)
utils::write.csv(summary_by_scenario, file.path(output_dir, "summary_by_scenario.csv"), row.names = FALSE)

cat("\nWrote: ", file.path(output_dir, "results_all.rds"), "\n", sep = "")
cat("Wrote: ", file.path(output_dir, "scenario_factor_levels.csv"), "\n", sep = "")
cat("Wrote: ", file.path(output_dir, "summary_by_scenario.csv"), "\n", sep = "")

# ---- Printed summary: median and spread (IQR) ----
# cat("\n--- Summary across simulations (median [Q1, Q3]) ---\n")
# for (sid in unique(results_all$scenario_id)) {
#   sub <- results_all[results_all$scenario_id == sid, , drop = FALSE]
#   cat("\n** ", as.character(sub$scenario_label[1L]), " **\n", sep = "")
#   for (nm in c("MAE", "RMSE", "RelativeMAE", "JS_bits", "TVD", "Cor_cells")) {
#     v <- sub[[nm]]
#     qs <- stats::quantile(v, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
#     cat(sprintf(
#       "  %-14s %8.4f  [%8.4f , %8.4f]\n",
#       nm, qs[2], qs[1], qs[3]
#     ))
#   }
# }

# ---- Static plots: crossed-factor exploration ----
prof_levels <- c("222", "332", "2222", "3322")
heat_df <- aggregate(
  JS_bits ~ drift + profile + constraint,
  data = results_all,
  FUN = stats::median
)
heat_df$profile <- factor(heat_df$profile, levels = prof_levels)
heat_df$drift_f <- factor(heat_df$drift, levels = sort(unique(heat_df$drift)))

p_heat <- ggplot(heat_df, aes(drift_f, profile, fill = JS_bits)) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_wrap(~constraint, nrow = 1L) +
  scale_fill_viridis_c(option = "C", end = 0.9) +
  labs(
    title = "Median JS (bits) by drift and table profile",
    subtitle = "Crossed factors; warmer = farther from simulated truth on average.",
    x = "Drift intensity",
    y = "Profile (dimension ladder)",
    fill = "median\nJS bits"
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank())
#print(p_heat)

line_df <- aggregate(
  JS_bits ~ drift + profile + constraint,
  data = results_all,
  FUN = stats::median
)
line_df$profile <- factor(line_df$profile, levels = prof_levels)

p_line <- ggplot(line_df, aes(drift, JS_bits, color = profile, group = profile)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~constraint, nrow = 1L) +
  labs(
    title = "Drift vs median Jensen–Shannon (bits)",
    subtitle = "Lines show whether coarser or finer tables diverge more as drift increases.",
    x = "Drift intensity",
    y = "Median JS (bits)",
    color = "Profile"
  ) +
  theme_minimal()
#print(p_line)

box_df <- results_all
p_box <- ggplot(box_df, aes(profile, JS_bits, fill = profile)) +
  geom_boxplot(outlier.alpha = 0.2) +
  facet_grid(constraint ~ drift_label, scales = "free_y", labeller = label_both) +
  labs(
    title = "Distribution of JS (bits) across replications",
    subtitle = "Facet columns = drift; rows = constraint pattern.",
    x = "Profile",
    y = "JS (bits)",
    fill = "Profile"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
#print(p_box)

plot_cols <- c("MAE", "RMSE", "RelativeMAE", "JS_bits", "TVD")
long_df <- data.frame()
for (pc in plot_cols) {
  long_df <- rbind(long_df, data.frame(
    scenario_label = results_all$scenario_label,
    metric = pc,
    value = results_all[[pc]],
    stringsAsFactors = FALSE
  ))
}
long_df$scenario_label <- factor(long_df$scenario_label, levels = unique(results_all$scenario_label))

# print(
#   ggplot(long_df, aes(scenario_label, value)) +
#     geom_boxplot(outlier.alpha = 0.35) +
#     facet_wrap(~metric, scales = "free_y", ncol = 2) +
#     labs(
#       title = "Metrics by scenario (full grid)",
#       subtitle = "Census seed + survey-style margins; truth is simulated.",
#       x = NULL,
#       y = NULL
#     ) +
#     theme_minimal() +
#     theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 6))
# )

# print(
#   ggplot(results_all, aes(scenario_label, JS_bits)) +
#     geom_boxplot(outlier.alpha = 0.35) +
#     labs(
#       title = "Jensen–Shannon distance (bits) by scenario",
#       subtitle = "Symmetric measure of how different fitted shares are from simulated truth.",
#       x = NULL,
#       y = "JS (bits)"
#     ) +
#     theme_minimal() +
#     theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 6))
# )

# PDF exports optional (Shiny app is primary for exploration). Set TRUE to write heatmap/lines/box PDFs.
SAVE_PDF_PLOTS <- FALSE
if (isTRUE(SAVE_PDF_PLOTS)) {
  ggsave(file.path(output_dir, "ipf_reliability_heatmap.pdf"), p_heat, width = 10, height = 7)
  ggsave(file.path(output_dir, "ipf_reliability_lines.pdf"), p_line, width = 10, height = 7)
  ggsave(file.path(output_dir, "ipf_reliability_boxplots.pdf"), p_box, width = 11, height = 8)
  cat("\nSaved PDFs under: ", output_dir, "\n", sep = "")
}
