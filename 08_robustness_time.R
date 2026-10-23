# Online Appendix OA.1: forecast horizon, decay factors, Minnesota prior and rolling-window VAR (time domain)
source("R/setup.R")

base <- readRDS(file.path(out_dir, "estimates", "baseline.rds"))
x <- base$x; p <- base$p; cn <- colnames(x)
dates <- as.Date(rownames(base$time$TCI))
summ <- function(L) {
  avg <- sapply(L, function(d) c(colMeans(d$TO), colMeans(d$FROM), colMeans(d$NET), TSI = mean(d$TCI)))
  k <- length(cn)
  f <- function(i) sprintf("%.2f [%.2f, %.2f]", median(avg[i, ]), min(avg[i, ]), max(avg[i, ]))
  out <- data.frame(TO = sapply(1:k, f), FROM = sapply(k + 1:k, f), NET = sapply(2 * k + 1:k, f), row.names = cn)
  rbind(out, TSI = c("", f(3 * k + 1), ""))
}
range_df <- function(M) data.frame(date = dates, lo = apply(M, 1, min), hi = apply(M, 1, max), med = apply(M, 1, median))
fig_range <- function(d, name, w = 8.15, h = 4.55, extra = NULL) {
  g <- ggplot(d, aes(date)) + geom_ribbon(aes(ymin = lo, ymax = hi), fill = "lightblue") +
    geom_line(aes(y = med), colour = "blue", linewidth = 0.6) + extra +
    scale_x_date(date_breaks = "4 years", date_labels = "%Y", expand = c(0.01, 0)) + theme_paper
  if ("country" %in% names(d)) g <- g + geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) + facet_wrap(~country, ncol = 3, scales = "free_y")
  save_fig(g, name, w, h)
}

# forecast horizon H = 1, ..., 24 (Table OA.1, Figures OA.1-OA.2)
hz <- lapply(1:24, function(h) conn_time(base$fit, h))
tab <- summ(hz); print(tab); write_table(tab, "tableOA1_horizon")
d <- cbind(range_df(sapply(hz[-1], function(d) d$TCI[, 1])), h1 = hz[[1]]$TCI[, 1])
fig_range(d, "figOA1_tsi_horizon", extra = geom_line(aes(y = h1), colour = "red", linewidth = 0.6))
d <- do.call(rbind, lapply(cn, function(k) cbind(range_df(sapply(hz[-1], function(d) d$NET[, k])), h1 = hz[[1]]$NET[, k], country = k)))
d$country <- factor(d$country, cn)
fig_range(d, "figOA2_net_horizon", extra = geom_line(aes(y = h1), colour = "red", linewidth = 0.5))

# decay factors (Table OA.2, Figures OA.3-OA.4); frequency results kept for 09
grid <- expand.grid(k1 = c(0.96, 0.97, 0.98, 0.99), k2 = c(0.96, 0.97, 0.98, 0.99))
dk <- lapply(seq_len(nrow(grid)), function(i) {
  fit <- tvpvar(x, p, l = c(grid$k1[i], grid$k2[i]))
  list(time = conn_time(fit), freq = conn_freq(fit))
})
saveRDS(list(grid = grid, res = dk), file.path(out_dir, "estimates", "decay.rds"))
tab <- summ(lapply(dk, `[[`, "time")); print(tab); write_table(tab, "tableOA2_decay")
fig_range(range_df(sapply(dk, function(d) d$time$TCI[, 1])), "figOA3_tsi_decay")
d <- do.call(rbind, lapply(cn, function(k) cbind(range_df(sapply(dk, function(d) d$time$NET[, k])), country = k)))
d$country <- factor(d$country, cn)
fig_range(d, "figOA4_net_decay")

# Minnesota prior (Table OA.3, Figures OA.5-OA.6)
fit_m <- tvpvar(x, p, prior = "minnesota")
mn <- list(time = conn_time(fit_m), freq = conn_freq(fit_m))
saveRDS(mn, file.path(out_dir, "estimates", "minnesota.rds"))
tab <- avg_table(mn$time$CT); print(tab); write_table(tab, "tableOA3_minnesota")
g <- ggplot(data.frame(date = dates, v = mn$time$TCI[, 1]), aes(date, v)) + geom_line(linewidth = 0.7) + theme_paper
save_fig(g, "figOA5_tsi_minnesota")
g <- ggplot(long_panel(mn$time$NET, dates), aes(date, v)) + geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
  geom_line(linewidth = 0.5) + facet_wrap(~country, ncol = 3) + theme_paper
save_fig(g, "figOA6_net_minnesota")

# rolling-window VAR, 200-month windows (Figures OA.7-OA.8)
rw <- ConnectednessApproach(x, nlag = p, nfore = H, window.size = prior_size + p, model = "VAR", connectedness = "Time",
                            corrected = FALSE, Connectedness_config = list(TimeConnectedness = list(generalized = TRUE)))
dw <- as.Date(rownames(rw$TCI))
g <- ggplot(data.frame(date = dw, v = rw$TCI[, 1]), aes(date, v)) + geom_line(linewidth = 0.7) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = c(0.01, 0)) + theme_paper
save_fig(g, "figOA7_tsi_rolling")
g <- ggplot(long_panel(rw$NET, dw), aes(date, v)) + geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
  geom_line(linewidth = 0.5) + facet_wrap(~country, ncol = 3) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y", expand = c(0.01, 0)) + theme_paper
save_fig(g, "figOA8_net_rolling")
cat("rolling window: first estimate", format(dw[1], "%YM%m"), "| TSI pre/post 2020M03:",
    round(mean(rw$TCI[dw < covid, 1]), 2), round(mean(rw$TCI[dw >= covid, 1]), 2), "\n")
