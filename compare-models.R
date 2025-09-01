# ─────────────────────────────────────────────────────────────
# Compare 2-param vs 4-param Stan fits (median curves, 12 months)
# Uses here::here() with your layout:
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

# (Optional) support for cmdstanr fits if you happen to have them saved
.has_cmdstanr <- requireNamespace("cmdstanr", quietly = TRUE)
.has_posterior <- requireNamespace("posterior", quietly = TRUE)

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

# ----- helpers: data -----------------------------------------
get_obs_monthly_cases <- function(region, n_take = n_months) {
  fn <- file.path(data_path, paste0("Region ", region, ".csv"))
  if (!file.exists(fn)) stop("Missing regional data CSV: ", fn)
  
  raw <- readr::read_csv(fn, show_col_types = FALSE) %>%
    mutate(Date = as.Date(Date))
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

# ----- helpers: fit loading ----------------------------------
# Safely coerce a loaded object into a usable fit (stanfit or CmdStanMCMC), else NULL
coerce_to_fit <- function(obj) {
  # Direct stanfit
  if (inherits(obj, "stanfit")) return(obj)
  
  # CmdStanR object
  if (.has_cmdstanr && inherits(obj, "CmdStanMCMC")) return(obj)
  
  # Common patterns where the fit is nested inside a list
  if (is.list(obj) && length(obj)) {
    # try a few name guesses
    candidates <- c("fit", "stanfit", "result", "object", "mod", "mcmc")
    for (nm in candidates) {
      if (!is.null(obj[[nm]])) {
        inner <- obj[[nm]]
        if (inherits(inner, "stanfit")) return(inner)
        if (.has_cmdstanr && inherits(inner, "CmdStanMCMC")) return(inner)
      }
    }
    # otherwise search all elements
    for (el in obj) {
      if (inherits(el, "stanfit")) return(el)
      if (.has_cmdstanr && inherits(el, "CmdStanMCMC")) return(el)
    }
  }
  
  # If it's a stanmodel (compiled only), not usable for extraction — return NULL
  if (inherits(obj, "stanmodel")) return(NULL)
  
  NULL
}

# Load a fit for a given region from a directory. Returns stanfit/CmdStanMCMC or NULL.
load_fit_for_region <- function(dir_path, region) {
  if (!dir.exists(dir_path)) return(NULL)
  files <- list.files(dir_path, pattern = "\\.(RDS|rds|RData)$", full.names = TRUE, ignore.case = TRUE)
  if (!length(files)) return(NULL)
  
  # Prefer files whose basename contains the region; fall back to all files
  score <- str_detect(tolower(basename(files)), tolower(region))
  files <- c(files[score], files[!score])
  
  for (f in files) {
    obj <- try({
      if (grepl("\\.RData$", f, ignore.case = TRUE)) {
        e <- new.env(parent = emptyenv())
        load(f, envir = e)
        # look for a usable fit in the environment
        objs <- mget(ls(e), envir = e)
        # Try each object
        found <- NULL
        for (nm in names(objs)) {
          found <- coerce_to_fit(objs[[nm]])
          if (!is.null(found)) break
        }
        found
      } else {
        coerce_to_fit(readRDS(f))
      }
    }, silent = TRUE)
    
    if (!inherits(obj, "try-error") && !is.null(obj)) return(obj)
  }
  NULL
}

# ----- helpers: extract median predictions --------------------
get_pred_median <- function(fit_obj, par_candidates = c("pred_cases","y_hat","yhat","pred_y","cases_pred"),
                            n = n_months) {
  # Case A: rstan::stanfit
  if (inherits(fit_obj, "stanfit")) {
    # first try extracting draws
    for (p in par_candidates) {
      ex <- try(rstan::extract(fit_obj, pars = p), silent = TRUE)
      if (!inherits(ex, "try-error") && !is.null(ex) && !is.null(ex[[1]])) {
        arr <- ex[[1]]
        # arr can be vector, matrix (iter x T), or 3D; reduce to time dim
        med <- if (is.null(dim(arr))) {
          median(arr, na.rm = TRUE)
        } else if (length(dim(arr)) == 2) {
          apply(arr, 2, median, na.rm = TRUE)
        } else {
          # assume dims: iterations x chains x time OR iterations x time x ...
          time_dim <- which.max(dim(arr)) # heuristic: take longest as time
          med <- apply(arr, time_dim, median, na.rm = TRUE)
          med
        }
        return(as.numeric(med[seq_len(min(n, length(med)))]))
      }
    }
    # fallback to summary
    for (p in par_candidates) {
      sm <- try(rstan::summary(fit_obj, pars = p, probs = 0.5)$summary, silent = TRUE)
      if (!inherits(sm, "try-error") && !is.null(sm) && nrow(sm) > 0) {
        cn <- colnames(sm)
        if ("50%" %in% cn) {
          med <- sm[, "50%"]
        } else if ("X50." %in% cn) {
          med <- sm[, "X50."]
        } else {
          next
        }
        return(as.numeric(med[seq_len(min(n, length(med)))]))
      }
    }
    stop("Could not find predicted-cases parameter in stanfit. Tried: ",
         paste(par_candidates, collapse = ", "))
  }
  
  # Case B: CmdStanR
  if (.has_cmdstanr && inherits(fit_obj, "CmdStanMCMC")) {
    if (!.has_posterior) stop("Found CmdStanMCMC fit but 'posterior' package is not installed.")
    for (p in par_candidates) {
      has_var <- try(p %in% fit_obj$metadata()$variables(), silent = TRUE)
      if (!inherits(has_var, "try-error") && isTRUE(has_var)) {
        draws <- fit_obj$draws(variables = p)  # draws_array
        # Convert to matrix: iterations*chains x time
        dm <- posterior::as_draws_matrix(draws)
        # If p is vector/time: columns are p[1], p[2], ...
        med <- apply(dm, 2, median, na.rm = TRUE)
        return(as.numeric(med[seq_len(min(n, length(med)))]))
      }
    }
    stop("Could not find predicted-cases parameter in CmdStanMCMC fit. Tried: ",
         paste(par_candidates, collapse = ", "))
  }
  
  stop("Unsupported fit object (need stanfit or CmdStanMCMC). You may have saved a compiled 'stanmodel' instead of a fitted object.")
}

# ----- metrics ------------------------------------------------
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
  
  fit2 <- load_fit_for_region(fit_dir_2, region)
  fit4 <- load_fit_for_region(fit_dir_4, region)
  
  if (is.null(fit2) || is.null(fit4)) {
    message(sprintf("[SKIP] Missing usable fit(s) for %s | 2-param usable: %s | 4-param usable: %s",
                    region, !is.null(fit2), !is.null(fit4)))
    next
  }
  
  yhat2 <- try(get_pred_median(fit2, n = length(y)), silent = TRUE)
  yhat4 <- try(get_pred_median(fit4, n = length(y)), silent = TRUE)
  
  if (inherits(yhat2, "try-error") || inherits(yhat4, "try-error")) {
    message(sprintf("[SKIP] Could not extract predictions for %s. 2-param ok: %s | 4-param ok: %s",
                    region, !inherits(yhat2, "try-error"), !inherits(yhat4, "try-error")))
    next
  }
  
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

# ----- bind & pool --------------------------------------------
if (length(monthly_long) == 0L) {
  warning("No regions produced results. Check that your *.RDS/*.RData files contain fitted objects (stanfit or CmdStanMCMC) and that filenames include region names.")
  monthly_long_tbl   <- tibble(region=character(), month=as.Date(character()),
                               obs=numeric(), model=character(), pred=numeric())
  metrics_region_tbl <- tibble()
  pooled_by_model    <- tibble()
} else {
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
}

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
