# ─────────────────────────────────────────────────────────────
# Compare 2-param vs 4-param Stan fits (median curves, 12 months)
# Uses here::here() with your new layout:
#   2-param-fit/   and   4-param-fit/
# Outputs (under function-plots/model-compare/):
#   - model_compare_monthly_long.csv
#   - model_compare_metrics_by_region.csv
#   - model_compare_overlay.png   (optional, uncomment to save)
# ─────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(rstan)
  library(stringr)
  library(here)
})

# ----- paths --------------------------------------------------
root_dir  <- here::here()
fit_dir_2 <- here("2-param-fit")
fit_dir_4 <- here("4-param-fit")

# Candidate locations for the regional CSVs. Adjust/add if your data live elsewhere.
data_path_opts <- c(
  here("RegionalData"),
  here("R", "fits", "RegionalData")
)
data_path <- data_path_opts[dir.exists(data_path_opts)][1]
if (is.na(data_path)) stop("Could not find RegionalData folder. Update `data_path_opts`.")

out_dir <- here("function-plots", "model-compare")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ----- config -------------------------------------------------
regions  <- c("Centre","Adamaoua","Nord","Littoral","Sud","Ouest","NordOuest","SudOuest","Est")
n_months <- 12

# ----- helpers ------------------------------------------------
find_fit_rds <- function(dir_path, region) {
  files <- list.files(dir_path, pattern = "\\.RDS$", full.names = TRUE, ignore.case = TRUE)
  if (!length(files)) return(NA_character_)
  hits <- files[str_detect(tolower(basename(files)), tolower(region))]
  if (length(hits) >= 1) return(hits[[1]])
  if (length(files) == 1) return(files[[1]])
  NA_character_
}

get_obs_monthly_cases <- function(region, n_take = n_months) {
  fn <- file.path(data_path, paste0("Region ", region, ".csv"))
  if (!file.exists(fn)) stop("Missing regional data CSV: ", fn)
  
  raw <- read.csv(fn) %>% mutate(Date = as.Date(Date))
  raw_2019 <- raw %>% filter(format(Date, "%Y") == "2019")
  
  l5 <- raw_2019 %>% filter(Age_group == "Age_L5") %>% pull(Cases)
  g5 <- raw_2019 %>% filter(Age_group == "Age_G5") %>% pull(Cases)
  
  m   <- min(length(l5), length(g5), n_take)
  y   <- as.numeric(l5[seq_len(m)] + g5[seq_len(m)])
  mon <- raw_2019 %>%
    transmute(Month = floor_date(Date, "month")) %>%
    distinct() %>% slice_head(n = m) %>% pull(Month)
  
  list(y = y, months = mon)
}

get_pred_median <- function(fit_obj, n = n_months) {
  par_candidates <- c("pred_cases","y_hat","yhat","pred_y","cases_pred")
  for (p in par_candidates) {
    if (p %in% names(rstan::extract(fit_obj))) {
      arr <- rstan::extract(fit_obj, pars = p)[[1]]  # draws x time
      med <- apply(arr, 2, median)
      return(as.numeric(med[seq_len(min(n, length(med)))]))
    }
  }
  for (p in par_candidates) {
    sm <- try(rstan::summary(fit_obj, pars = p, probs = 0.5)$summary, silent = TRUE)
    if (!inherits(sm, "try-error")) {
      med <- sm[, "50%"] %||% sm[, "X50."]
      return(as.numeric(med[seq_len(min(n, length(med)))]))
    }
  }
  stop("Could not find a predicted-cases parameter in stanfit (tried: ",
       paste(par_candidates, collapse = ", "), ").")
}

metrics_vec <- function(y, yhat) {
  stopifnot(length(y) == length(yhat))
  sse  <- sum((y - yhat)^2)
  sst  <- sum((y - mean(y))^2)
  tibble(
    n        = length(y),
    RMSE     = sqrt(mean((y - yhat)^2)),
    MAE      = mean(abs(y - yhat)),
    MAPE_pct = mean(abs((y - yhat) / ifelse(y == 0, NA, y)), na.rm = TRUE) * 100,
    R2_SSE   = if (sst == 0) NA_real_ else 1 - sse / sst,
    R2_corr  = suppressWarnings(cor(y, yhat)^2)
  )
}

# ----- main loop ----------------------------------------------
monthly_long   <- list()
metrics_region <- list()

for (region in regions) {
  obs <- get_obs_monthly_cases(region, n_take = n_months)
  y   <- obs$y
  mo  <- obs$months
  
  if (length(y) < 3) {
    message(sprintf("[SKIP] %s has <3 months of data.", region))
    next
  }
  
  f2_path <- find_fit_rds(fit_dir_2, region)
  f4_path <- find_fit_rds(fit_dir_4, region)
  if (is.na(f2_path) || is.na(f4_path)) {
    message(sprintf("[SKIP] Missing fit(s) for %s | 2-param: %s | 4-param: %s",
                    region, !is.na(f2_path), !is.na(f4_path)))
    next
  }
  
  fit2 <- readRDS(f2_path)
  fit4 <- readRDS(f4_path)
  
  yhat2 <- get_pred_median(fit2, n = length(y))
  yhat4 <- get_pred_median(fit4, n = length(y))
  
  monthly_long[[region]] <- tibble(
    region = region,
    month  = mo,
    obs    = y,
    `2param` = yhat2,
    `4param` = yhat4
  ) %>%
    pivot_longer(cols = c(`2param`,`4param`),
                 names_to = "model", values_to = "pred")
  
  m2 <- metrics_vec(y, yhat2) %>% mutate(model = "2param", region = region, .before = 1)
  m4 <- metrics_vec(y, yhat4) %>% mutate(model = "4param", region = region, .before = 1)
  
  metrics_region[[region]] <- bind_rows(m2, m4) %>%
    pivot_wider(names_from = model, values_from = c(RMSE, MAE, MAPE_pct, R2_SSE, R2_corr, n)) %>%
    mutate(
      d_RMSE     = RMSE_4param     - RMSE_2param,
      d_MAE      = MAE_4param      - MAE_2param,
      d_MAPE_pct = MAPE_pct_4param - MAPE_pct_2param,
      d_R2_SSE   = R2_SSE_4param   - R2_SSE_2param,
      d_R2_corr  = R2_corr_4param  - R2_corr_2param
    )
}

monthly_long_tbl   <- bind_rows(monthly_long)
metrics_region_tbl <- bind_rows(metrics_region)

# Pooled metrics (simple pooling across all regions/months)
pooled_by_model <- monthly_long_tbl %>%
  group_by(model) %>%
  summarise(
    RMSE     = sqrt(mean((obs - pred)^2)),
    MAE      = mean(abs(obs - pred)),
    MAPE_pct = mean(abs((obs - pred) / ifelse(obs == 0, NA, obs)), na.rm = TRUE) * 100,
    R2_SSE   = { y <- obs; yhat <- pred; sst <- sum((y - mean(y))^2); ifelse(sst == 0, NA_real_, 1 - sum((y - yhat)^2)/sst) },
    R2_corr  = suppressWarnings(cor(obs, pred)^2),
    .groups = "drop"
  )

# ----- save outputs -------------------------------------------
readr::write_csv(monthly_long_tbl,   file.path(out_dir, "model_compare_monthly_long.csv"))
readr::write_csv(metrics_region_tbl, file.path(out_dir, "model_compare_metrics_by_region.csv"))

print(pooled_by_model)
message("Wrote:\n - ", file.path(out_dir, "model_compare_monthly_long.csv"),
        "\n - ", file.path(out_dir, "model_compare_metrics_by_region.csv"))

# ----- optional quick plot ------------------------------------
# g <- monthly_long_tbl %>%
#   mutate(model = factor(model, levels = c("2param","4param"))) %>%
#   ggplot(aes(month, pred, group = model)) +
#   geom_line(aes(linetype = model)) +
#   geom_point(aes(y = obs), inherit.aes = FALSE) +
#   facet_wrap(~region, scales = "free_y") +
#   labs(x = NULL, y = "Cases", linetype = "Model",
#        title = "Observed vs Predicted (median)\n2-param vs 4-param") +
#   theme_bw()
# ggsave(file.path(out_dir, "model_compare_overlay.png"), g, width = 12, height = 8, dpi = 200)
