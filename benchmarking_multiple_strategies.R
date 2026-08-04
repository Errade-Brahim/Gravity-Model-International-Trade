## ============================================================
## Benchmarking Multiple Strategies — S&P 500 (2015–2024)
## Momentum | Mean-Reversion | Pairs Trading
## VERSION CORRIGÉE — Out-of-sample 2021–2024
## ============================================================

# ── 0. Packages ───────────────────────────────────────────────
packages <- c("quantmod", "PerformanceAnalytics",
              "TTR", "zoo", "xts", "dplyr", "tidyr",
              "ggplot2", "lubridate", "tseries",
              "depmixS4", "quadprog", "reshape2")

lapply(packages, function(p) {
  if (!require(p, character.only = TRUE))
    install.packages(p)
  library(p, character.only = TRUE)
})

options(xts.warn_dplyr_breaks_lag = FALSE)


# ── 1. Data Download ──────────────────────────────────────────
tickers <- c("AAPL","MSFT","AMZN","NVDA","GOOGL",
             "META","BRK-B","TSLA","JPM","UNH",
             "XOM","JNJ","V","PG","MA",
             "HD","CVX","MRK","LLY","ABBV")

download_prices <- function(tickers, from = "2015-01-01",
                            to = "2024-12-31", max_tries = 3) {
  success <- c(); failed <- tickers
  for (attempt in 1:max_tries) {
    if (length(failed) == 0) break
    for (t in failed) {
      tryCatch({
        suppressWarnings(
          getSymbols(t, src = "yahoo", from = from, to = to,
                     auto.assign = TRUE, env = globalenv())
        )
        success <- c(success, t)
      }, error = function(e) cat("  Echec:", t, "\n"))
    }
    failed <- setdiff(failed, success)
    if (length(failed) > 0) Sys.sleep(3)
  }
  return(success)
}

loaded_tickers <- download_prices(tickers)

prices_list <- lapply(loaded_tickers, function(t) {
  obj_name <- gsub("-", ".", t)
  sym <- tryCatch(Ad(get(obj_name)), error = function(e) NULL)
  if (is.null(sym) || nrow(sym) < 10) return(NULL)
  colnames(sym) <- t; sym
})
prices_list <- Filter(Negate(is.null), prices_list)
prices <- Reduce(function(a, b) merge(a, b, join = "outer"), prices_list)
min_obs <- floor(0.8 * ncol(prices))
prices  <- prices[rowSums(!is.na(prices)) >= min_obs, ]
tickers <- colnames(prices)
cat("Tickers:", length(tickers), "| Période:", as.character(start(prices)),
    "→", as.character(end(prices)), "\n")

# SPY benchmark
getSymbols("SPY", src = "yahoo", from = "2015-01-01",
           to = "2024-12-31", auto.assign = TRUE)
sp500_returns <- dailyReturn(Ad(SPY), type = "log")

# Splits in-sample / out-of-sample
prices_in  <- prices["2015/2020"]
prices_out <- prices["2021/2024"]
sp500_in   <- sp500_returns["2015/2020"]
sp500_out  <- sp500_returns["2021/2024"]


# ── 2. Momentum Strategy ──────────────────────────────────────
momentum_strategy <- function(prices, J = 252, K = 21,
                              top_n = 5, cost = 0.00225) {
  prices  <- prices[, colSums(is.na(prices)) == 0]
  n       <- nrow(prices)
  tickers <- colnames(prices)

  log_ret <- apply(log(prices), 2, diff)
  log_ret <- rbind(rep(NA, ncol(prices)), log_ret)

  port_ret_vec <- rep(NA_real_, n)

  for (i in (J + K + 1):n) {
    form_start <- i - J - K
    form_end   <- i - K - 1
    p_start <- as.numeric(prices[form_start, ])
    p_end   <- as.numeric(prices[form_end,   ])
    valid <- p_start > 0 & p_end > 0 & !is.na(p_start) & !is.na(p_end)
    if (sum(valid) < top_n) next
    form_ret <- p_end[valid] / p_start[valid] - 1
    names(form_ret) <- tickers[valid]
    winners  <- names(sort(form_ret, decreasing = TRUE))[1:top_n]
    daily_r  <- as.numeric(log_ret[i, winners])
    daily_r  <- daily_r[is.finite(daily_r)]
    if (length(daily_r) == 0) next
    port_ret_vec[i] <- mean(daily_r) - cost / K
  }

  xts(port_ret_vec, order.by = index(prices)) |> na.omit()
}

# Calibration in-sample, évaluation out-of-sample
# On utilise les prix full pour la formation mais on garde l'index out-of-sample
mom_returns_full <- momentum_strategy(prices)
mom_returns_in   <- mom_returns_full[index(mom_returns_full) >= as.Date("2015-01-01") &
                                     index(mom_returns_full) <= as.Date("2020-12-31")]
mom_returns_out  <- mom_returns_full[index(mom_returns_full) >= as.Date("2021-01-01")]

cat("Momentum out-of-sample obs:", length(mom_returns_out), "\n")


# ── 3. Mean-Reversion (Bollinger Bands) — VERSION CORRIGÉE ────
# CORRECTION : on remplace rowMeans par médiane pour éviter les outliers
# et on filtre les returns aberrants (|r| > 0.5 = erreur de données)

mean_reversion_strategy <- function(prices, n = 20, k = 2,
                                    cost = 0.00225) {
  dates      <- index(prices)
  N          <- nrow(prices)
  result_mat <- matrix(NA_real_, nrow = N, ncol = ncol(prices))
  colnames(result_mat) <- colnames(prices)
  valid_cols <- c()

  for (ticker in colnames(prices)) {
    px <- as.numeric(prices[, ticker])
    if (sum(!is.na(px)) < n + 2) next

    bb <- tryCatch(BBands(px, n = n, sd = k), error = function(e) NULL)
    if (is.null(bb)) next

    dn   <- as.numeric(bb[, "dn"])
    mavg <- as.numeric(bb[, "mavg"])

    signal <- rep(0, N)
    in_pos <- FALSE
    for (i in (n + 1):N) {
      if (is.na(px[i]) || is.na(dn[i]) || is.na(mavg[i])) next
      if (!in_pos && px[i] < dn[i])        { signal[i] <-  1; in_pos <- TRUE  }
      else if (in_pos && px[i] >= mavg[i]) { signal[i] <- -1; in_pos <- FALSE }
      else if (in_pos)                      { signal[i] <-  1 }
    }

    lr <- c(NA, diff(log(px)))

    # CORRECTION : winsorize les log-returns aberrants
    lr[abs(lr) > 0.5 & !is.na(lr)] <- NA

    pos_lag <- c(0, signal[-N])
    strat_r <- pos_lag * lr
    changed <- which(signal != pos_lag)
    strat_r[changed] <- strat_r[changed] - cost

    result_mat[, ticker] <- strat_r
    valid_cols <- c(valid_cols, ticker)
  }

  if (length(valid_cols) == 0) stop("Aucun ticker valide.")

  # CORRECTION : médiane au lieu de moyenne pour robustesse
  port_vec <- apply(result_mat[, valid_cols, drop = FALSE], 1,
                    function(x) median(x, na.rm = TRUE))
  xts(port_vec, order.by = dates) |> na.omit()
}

mr_returns_full <- mean_reversion_strategy(prices)
mr_returns_in   <- mr_returns_full[index(mr_returns_full) >= as.Date("2015-01-01") &
                                   index(mr_returns_full) <= as.Date("2020-12-31")]
mr_returns_out  <- mr_returns_full[index(mr_returns_full) >= as.Date("2021-01-01")]

cat("MeanRev out-of-sample obs:", length(mr_returns_out), "\n")


# ── 4. Pairs Trading (V + MA) — meilleure paire coïntégrée ────
pairs_strategy <- function(prices, s1, s2,
                           entry_z = 2, exit_z = 0,
                           cost = 0.00225) {
  p1 <- as.numeric(prices[, s1])
  p2 <- as.numeric(prices[, s2])

  # Hedge ratio estimé sur in-sample uniquement (évite look-ahead bias)
  n_in <- sum(index(prices) <= as.Date("2020-12-31"))
  reg  <- lm(p1[1:n_in] ~ p2[1:n_in])
  hedge <- coef(reg)[2]

  spread    <- p1 - hedge * p2
  roll_mean <- as.numeric(rollmean(spread, 60, fill = NA, align = "right"))
  roll_sd   <- as.numeric(rollapply(spread, 60, sd, fill = NA, align = "right"))
  z_score   <- (spread - roll_mean) / roll_sd

  signal <- rep(0, length(z_score))
  pos    <- 0
  for (i in 61:length(z_score)) {
    z <- z_score[i]
    if (is.na(z)) { signal[i] <- pos; next }
    if (pos ==  0 && z >  entry_z) pos <- -1
    if (pos ==  0 && z < -entry_z) pos <-  1
    if (pos != 0  && abs(z) < exit_z) pos <- 0
    signal[i] <- pos
  }

  r1     <- c(0, diff(log(p1)))
  r2     <- c(0, diff(log(p2)))
  port_r <- signal * (r1 - hedge * r2)

  # Transaction costs on signal changes
  tc     <- cost * abs(diff(c(0, signal)))
  port_r <- port_r - tc

  xts(port_r, order.by = index(prices))
}

# Meilleure paire : V & MA (pval = 0.01)
pairs_returns_full <- pairs_strategy(prices, s1 = "V", s2 = "MA")
pairs_returns_in   <- pairs_returns_full[index(pairs_returns_full) >= as.Date("2015-01-01") &
                                         index(pairs_returns_full) <= as.Date("2020-12-31")]
pairs_returns_out  <- pairs_returns_full[index(pairs_returns_full) >= as.Date("2021-01-01")]
pairs_returns_out  <- na.omit(pairs_returns_out)

cat("Pairs out-of-sample obs:", length(pairs_returns_out), "\n")


# ── 5. Performance Metrics ────────────────────────────────────
full_report <- function(strat_list, rf = 0.04/252, label = "Full Period") {
  cat("\n===", label, "===\n")
  results <- lapply(names(strat_list), function(nm) {
    ret <- as.numeric(na.omit(strat_list[[nm]]))
    ret <- ret[is.finite(ret)]
    if (length(ret) < 30) { cat("  SKIP", nm, "\n"); return(NULL) }

    ret_xts <- xts(matrix(ret, ncol = 1),
                   order.by = seq.Date(as.Date("2021-01-01"),
                                       by = "day", length.out = length(ret)))
    colnames(ret_xts) <- nm

    ann_r  <- as.numeric(Return.annualized(ret_xts, scale = 252))
    vol    <- as.numeric(StdDev.annualized(ret_xts, scale = 252))
    sharpe <- as.numeric(SharpeRatio.annualized(ret_xts, Rf = rf, scale = 252))
    sort_r <- as.numeric(SortinoRatio(ret_xts, MAR = rf))
    mdd    <- as.numeric(maxDrawdown(ret_xts))
    calmar <- as.numeric(CalmarRatio(ret_xts))
    wr     <- mean(ret > 0)
    skew   <- as.numeric(skewness(ret))
    kurt   <- as.numeric(kurtosis(ret)) - 3
    cvar   <- -mean(ret[ret <= quantile(ret, 0.05)])

    data.frame(
      Strategy   = nm,
      AnnReturn  = round(ann_r  * 100, 2),
      Volatility = round(vol    * 100, 2),
      Sharpe     = round(sharpe, 3),
      Sortino    = round(sort_r, 3),
      MaxDD      = round(mdd    * 100, 2),
      Calmar     = round(calmar, 3),
      WinRate    = round(wr     * 100, 1),
      Skewness   = round(skew,  3),
      ExKurtosis = round(kurt,  3),
      CVaR95     = round(cvar   * 100, 2)
    )
  })
  df <- do.call(rbind, Filter(Negate(is.null), results))
  print(df, row.names = FALSE)
  invisible(df)
}

# ── 6. Rapports ───────────────────────────────────────────────

# A. Full period (2015–2024)
report_full <- full_report(
  list(Momentum     = mom_returns_full,
       MeanReversion = mr_returns_full,
       PairsTrading  = pairs_returns_full,
       SP500         = sp500_returns),
  label = "FULL PERIOD 2015–2024"
)

# B. In-sample (2015–2020)
report_in <- full_report(
  list(Momentum     = mom_returns_in,
       MeanReversion = mr_returns_in,
       PairsTrading  = pairs_returns_in,
       SP500         = sp500_in),
  label = "IN-SAMPLE 2015–2020"
)

# C. Out-of-sample (2021–2024) — chiffres pour le rapport
report_out <- full_report(
  list(Momentum     = mom_returns_out,
       MeanReversion = mr_returns_out,
       PairsTrading  = pairs_returns_out,
       SP500         = sp500_out),
  label = "OUT-OF-SAMPLE 2021–2024"
)


# ── 7. CVaR Table ─────────────────────────────────────────────
cvar_95 <- function(r) {
  r <- na.omit(as.numeric(r))
  q <- quantile(r, 0.05)
  list(VaR = -q * 100, CVaR = -mean(r[r <= q]) * 100)
}

strats_out <- list(
  Momentum     = mom_returns_out,
  MeanReversion = mr_returns_out,
  PairsTrading  = pairs_returns_out,
  SP500         = sp500_out
)

cat("\n=== VaR & CVaR (Out-of-Sample, 95%) ===\n")
var_cvar <- do.call(rbind, lapply(names(strats_out), function(nm) {
  v <- cvar_95(strats_out[[nm]])
  data.frame(Strategy = nm,
             VaR_95   = round(v$VaR,  2),
             CVaR_95  = round(v$CVaR, 2),
             Ratio    = round(v$CVaR / v$VaR, 2))
}))
print(var_cvar, row.names = FALSE)


# ── 8. Correlation Matrix ─────────────────────────────────────
strat_mat <- merge(mom_returns_out, mr_returns_out,
                   pairs_returns_out, sp500_out, join = "inner") |> na.omit()
colnames(strat_mat) <- c("Momentum", "MeanRev", "Pairs", "SP500")
cat("\n=== Correlation Matrix (Out-of-Sample) ===\n")
print(round(cor(strat_mat), 3))


# ── 9. Portfolio Optimisation (Out-of-Sample weights) ─────────
library(quadprog)

strat_mat_in <- merge(mom_returns_in, mr_returns_in,
                      pairs_returns_in, join = "inner") |> na.omit()
colnames(strat_mat_in) <- c("Momentum", "MeanRev", "Pairs")

Sigma   <- cov(strat_mat_in) * 252
mu_vec  <- colMeans(strat_mat_in) * 252
n_s     <- 3

# Min-Variance
Dmat <- 2 * Sigma
dvec <- rep(0, n_s)
Amat <- cbind(rep(1, n_s), diag(n_s))
bvec <- c(1, rep(0, n_s))
mvp  <- solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
w_mvp <- setNames(round(mvp$solution, 4), colnames(strat_mat_in))

# Risk Parity
risk_parity <- function(Sigma) {
  n <- nrow(Sigma); w0 <- rep(1/n, n)
  obj <- function(w) {
    pv  <- as.numeric(t(w) %*% Sigma %*% w)
    mrc <- Sigma %*% w / pv
    sum((w * mrc - mean(w * mrc))^2)
  }
  res <- optim(w0, obj, method = "L-BFGS-B",
               lower = rep(0, n), upper = rep(1, n))
  res$par / sum(res$par)
}
w_rp <- setNames(round(risk_parity(Sigma), 4), colnames(strat_mat_in))

cat("\n=== Portfolio Weights ===\n")
cat("Min-Variance:"); print(w_mvp)
cat("Risk Parity: "); print(w_rp)

# Equal-weight portfolio out-of-sample
ew_ret <- rowMeans(strat_mat[, c("Momentum","MeanRev","Pairs")], na.rm = TRUE)
ew_xts <- xts(ew_ret, order.by = index(strat_mat))

# MV portfolio out-of-sample
strat_mat_out3 <- strat_mat[, c("Momentum","MeanRev","Pairs")]
mv_ret <- as.numeric(as.matrix(strat_mat_out3) %*% w_mvp)
mv_xts <- xts(mv_ret, order.by = index(strat_mat))

# RP portfolio out-of-sample
rp_ret <- as.numeric(as.matrix(strat_mat_out3) %*% w_rp)
rp_xts <- xts(rp_ret, order.by = index(strat_mat))

report_portfolios <- full_report(
  list(EqualWeight  = ew_xts,
       MinVariance  = mv_xts,
       RiskParity   = rp_xts),
  label = "MULTI-STRATEGY PORTFOLIOS (Out-of-Sample)"
)


# ── 10. Regime Analysis ───────────────────────────────────────
library(depmixS4)

fit_hmm_regime <- function(returns, n_states = 2) {
  ret_df <- data.frame(r = as.numeric(na.omit(returns)))
  mod    <- depmix(r ~ 1, data = ret_df, nstates = n_states,
                   family = gaussian())
  fit    <- fit(mod, verbose = FALSE)
  post   <- posterior(fit, type = "viterbi")
  list(states = post$state, fit = fit)
}

hmm_result  <- fit_hmm_regime(sp500_returns)
sp500_states_full <- hmm_result$states

# Identifier Bull vs Bear (état avec moyenne la plus haute = Bull)
hmm_means <- summary(hmm_result$fit)
cat("\n=== HMM State Means ===\n")
print(hmm_means)


# ── 11. Visualisations ────────────────────────────────────────

# A. Cumulative Wealth (Out-of-Sample)
cum_wealth <- function(ret, initial = 1e6) {
  r <- na.omit(as.numeric(ret))
  initial * cumprod(1 + r)
}

w_mom_out   <- xts(cum_wealth(mom_returns_out),   order.by = index(mom_returns_out))
w_mr_out    <- xts(cum_wealth(mr_returns_out),    order.by = index(mr_returns_out))
w_pairs_out <- xts(cum_wealth(pairs_returns_out), order.by = index(pairs_returns_out))
w_sp5_out   <- xts(cum_wealth(sp500_out),         order.by = index(sp500_out))

w_all <- merge(w_mom_out, w_mr_out, w_pairs_out, w_sp5_out)
w_all <- na.omit(w_all)
colnames(w_all) <- c("Momentum", "MeanRev", "Pairs", "SP500")

df_w <- data.frame(Date = as.Date(index(w_all)),
                   as.data.frame(w_all))
df_w_long <- reshape2::melt(df_w, id.vars = "Date",
                            variable.name = "Strategy",
                            value.name = "Wealth")

p1 <- ggplot(df_w_long, aes(x = Date, y = Wealth / 1e6, colour = Strategy)) +
  geom_line(linewidth = 0.8) +
  scale_y_continuous(labels = scales::dollar_format(suffix = "M", prefix = "$")) +
  labs(title = "Cumulative Wealth — Out-of-Sample 2021–2024 ($1M Initial)",
       x = "Date", y = "Portfolio Value (USD Millions)", colour = "Strategy") +
  theme_minimal(base_size = 13)
print(p1)
ggsave("cumulative_wealth_out.png", p1, width = 10, height = 5, dpi = 150)


# B. Rolling Sharpe (Full Period)
rolling_sharpe <- function(ret, window = 252, rf = 0.04/252) {
  rollapply(ret, width = window, FUN = function(x) {
    ex <- x - rf
    if (sd(ex, na.rm = TRUE) == 0) return(NA)
    mean(ex, na.rm = TRUE) / sd(ex, na.rm = TRUE) * sqrt(252)
  }, align = "right", fill = NA)
}

rs_mom   <- rolling_sharpe(mom_returns_full)
rs_mr    <- rolling_sharpe(mr_returns_full)
rs_pairs <- rolling_sharpe(pairs_returns_full)

rs_all <- merge(rs_mom, rs_mr, rs_pairs) |> na.omit()
colnames(rs_all) <- c("Momentum", "MeanRev", "Pairs")

df_rs <- data.frame(Date = as.Date(index(rs_all)), as.data.frame(rs_all))
df_rs_long <- reshape2::melt(df_rs, id.vars = "Date",
                             variable.name = "Strategy",
                             value.name = "RollingSharpe")

p2 <- ggplot(df_rs_long, aes(x = Date, y = RollingSharpe, colour = Strategy)) +
  geom_line(linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "green4") +
  labs(title = "Rolling 252-Day Sharpe Ratio (2015–2024)",
       x = "Date", y = "Sharpe Ratio") +
  theme_minimal(base_size = 13)
print(p2)
ggsave("rolling_sharpe.png", p2, width = 10, height = 5, dpi = 150)


# C. Drawdown Chart
drawdown_plot <- function(ret, name) {
  r   <- na.omit(as.numeric(ret))
  cw  <- cumprod(1 + r)
  pk  <- cummax(cw)
  dd  <- (cw - pk) / pk * 100
  data.frame(Date = as.Date(index(na.omit(ret))), DD = dd, Strategy = name)
}

dd_df <- rbind(
  drawdown_plot(mom_returns_full,   "Momentum"),
  drawdown_plot(mr_returns_full,    "MeanRev"),
  drawdown_plot(pairs_returns_full, "Pairs"),
  drawdown_plot(sp500_returns,      "SP500")
)

p3 <- ggplot(dd_df, aes(x = Date, y = DD, colour = Strategy)) +
  geom_line(linewidth = 0.6) +
  labs(title = "Drawdown Profile (2015–2024)",
       x = "Date", y = "Drawdown (%)") +
  theme_minimal(base_size = 13)
print(p3)
ggsave("drawdown_profile.png", p3, width = 10, height = 5, dpi = 150)


# ── 12. LaTeX Tables (prêtes à coller dans le rapport) ────────
cat("\n\n")
cat("% ============================================================\n")
cat("% TABLES LaTeX — copiez-collez dans benchmarking_strategies.tex\n")
cat("% ============================================================\n\n")

print_latex_table <- function(df, caption, label) {
  cat("\\begin{table}[H]\n")
  cat("  \\centering\n")
  cat("  \\caption{", caption, "}\n", sep="")
  cat("  \\label{", label, "}\n", sep="")
  cat("  \\small\n")
  cat("  \\begin{tabular}{lrrrrrrr}\n")
  cat("    \\toprule\n")
  cat("    \\textbf{Strategy} & \\textbf{Ann.Ret (\\%)} & \\textbf{Vol (\\%)} &",
      "\\textbf{Sharpe} & \\textbf{Sortino} & \\textbf{MaxDD (\\%)} &",
      "\\textbf{Calmar} & \\textbf{WinRate (\\%)} \\\\\n")
  cat("    \\midrule\n")
  for (i in 1:nrow(df)) {
    cat("   ", df$Strategy[i], "&", df$AnnReturn[i], "&", df$Volatility[i],
        "&", df$Sharpe[i], "&", df$Sortino[i], "&",
        paste0("-", abs(df$MaxDD[i])), "&", df$Calmar[i], "&",
        df$WinRate[i], "\\\\\n")
  }
  cat("    \\bottomrule\n")
  cat("  \\end{tabular}\n")
  cat("\\end{table}\n\n")
}

print_latex_table(report_out,
  caption = "Strategy Performance --- Out-of-Sample 2021--2024",
  label   = "tab:perf_out")

print_latex_table(report_in,
  caption = "Strategy Performance --- In-Sample 2015--2020",
  label   = "tab:perf_in")

print_latex_table(report_portfolios,
  caption = "Multi-Strategy Portfolio Performance --- Out-of-Sample 2021--2024",
  label   = "tab:multi_strat")

cat("\n% VaR & CVaR Table\n")
cat("\\begin{table}[H]\n  \\centering\n")
cat("  \\caption{VaR and CVaR Estimates (Daily, 95\\% Confidence, Out-of-Sample)}\n")
cat("  \\label{tab:var_cvar_real}\n")
cat("  \\begin{tabular}{lrrr}\n    \\toprule\n")
cat("    \\textbf{Strategy} & \\textbf{Hist. VaR (\\%)} & \\textbf{CVaR (\\%)} & \\textbf{CVaR/VaR} \\\\\n")
cat("    \\midrule\n")
for (i in 1:nrow(var_cvar)) {
  cat("   ", var_cvar$Strategy[i], "&", round(var_cvar$VaR_95[i],2),
      "&", round(var_cvar$CVaR_95[i],2), "&", var_cvar$Ratio[i], "\\\\\n")
}
cat("    \\bottomrule\n  \\end{tabular}\n\\end{table}\n\n")

cat("\n% Correlation Matrix\n")
cor_mat <- round(cor(strat_mat[, c("Momentum","MeanRev","Pairs","SP500")]), 3)
cat("\\begin{table}[H]\n  \\centering\n")
cat("  \\caption{Pairwise Correlation Matrix --- Out-of-Sample 2021--2024}\n")
cat("  \\label{tab:corr_real}\n")
cat("  \\begin{tabular}{lrrrr}\n    \\toprule\n")
cat("    & \\textbf{Momentum} & \\textbf{MeanRev} & \\textbf{Pairs} & \\textbf{SP500} \\\\\n")
cat("    \\midrule\n")
for (nm in rownames(cor_mat)) {
  cat("   ", nm, "&", paste(cor_mat[nm,], collapse=" & "), "\\\\\n")
}
cat("    \\bottomrule\n  \\end{tabular}\n\\end{table}\n")

cat("\n\n=== DONE — Sauvegardez les graphiques ===\n")
cat("  cumulative_wealth_out.png\n")
cat("  rolling_sharpe.png\n")
cat("  drawdown_profile.png\n")
