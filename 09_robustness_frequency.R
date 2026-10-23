# Online Appendix OA.2: frequency-domain robustness
source("R/setup.R")

base <- readRDS(file.path(out_dir, "estimates", "baseline.rds"))
x <- base$x; p <- base$p; cn <- colnames(x)
dates <- as.Date(rownames(base$time$TCI))
yrs <- function(k) scale_x_date(date_breaks = paste(k, "years"), date_labels = "%Y", expand = c(0.01, 0))
two_bands <- function(S, L, f = dates) {
  d <- rbind(data.frame(date = f, band = "Short run", S), data.frame(date = f, band = "Long run", L))
  d$band <- factor(d$band, c("Short run", "Long run")); d
}
fig_fsi <- function(S, L, name, f = dates) {
  g <- ggplot(two_bands(data.frame(v = S), data.frame(v = L), f), aes(date, v, colour = band)) + geom_line(linewidth = 0.7) +
    scale_colour_manual(values = col_bands) + yrs(3) + theme_paper
  save_fig(g, name, 7.55, 4.45)
}
fig_net <- function(NS, NL, name, f = dates) {
  d <- rbind(long_panel(NS, f, "Short run"), long_panel(NL, f, "Long run")); d$band <- factor(d$band, c("Short run", "Long run"))
  g <- ggplot(d, aes(date, v, colour = band)) + geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
    geom_line(linewidth = 0.5) + scale_colour_manual(values = col_bands) + facet_wrap(~country, ncol = 3) + yrs(4) + theme_paper
  save_fig(g, name, 7.85, 4.45)
}
rng <- function(M) data.frame(lo = apply(M, 1, min), hi = apply(M, 1, max), med = apply(M, 1, median))

# decay factors (Figures OA.9-OA.10)
dk <- readRDS(file.path(out_dir, "estimates", "decay.rds"))$res
d <- two_bands(rng(sapply(dk, function(r) r$freq$TCI[, 2])), rng(sapply(dk, function(r) r$freq$TCI[, 3])))
g <- ggplot(d, aes(date, colour = band, fill = band)) + geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, colour = NA) +
  geom_line(aes(y = med), linewidth = 0.7) + scale_colour_manual(values = col_bands) + scale_fill_manual(values = col_bands) + yrs(3) + theme_paper
save_fig(g, "figOA9_fsi_decay", 7.55, 4.45)
d <- do.call(rbind, lapply(cn, function(k) cbind(two_bands(rng(sapply(dk, function(r) r$freq$NET[, k, 2])),
                                                           rng(sapply(dk, function(r) r$freq$NET[, k, 3]))), country = k)))
d$country <- factor(d$country, cn)
g <- ggplot(d, aes(date, colour = band, fill = band)) + geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, colour = NA) + geom_line(aes(y = med), linewidth = 0.5) +
  scale_colour_manual(values = col_bands) + scale_fill_manual(values = col_bands) + facet_wrap(~country, ncol = 3) + yrs(4) + theme_paper
save_fig(g, "figOA10_net_freq_decay", 7.85, 4.45)

# Minnesota prior (Figures OA.11-OA.12)
mn <- readRDS(file.path(out_dir, "estimates", "minnesota.rds"))$freq
fig_fsi(mn$TCI[, 2], mn$TCI[, 3], "figOA11_fsi_minnesota")
fig_net(mn$NET[, , 2], mn$NET[, , 3], "figOA12_net_freq_minnesota")

# three bands: 2-8, 8-20 and > 20 months (Table OA.4, Figures OA.13-OA.14)
f3 <- conn_freq(base$fit, bands3)
nb <- c("Short run", "Medium run", "Long run")
samples <- list("(a) Full sample" = rep(TRUE, length(dates)), "(b) Pre-COVID" = dates < covid, "(c) Post-COVID" = dates >= covid)
tab <- do.call(rbind, lapply(names(samples), function(s) do.call(rbind, lapply(1:3, function(b) {
  M <- apply(f3$CT[, , samples[[s]], b + 1], 1:2, mean) * 100; off <- M; diag(off) <- 0
  rbind(data.frame(sample = s, band = nb[b], measure = c("TO", "FROM", "NET"),
                   round(rbind(colSums(off), rowSums(off), colSums(off) - rowSums(off)), 2), check.names = FALSE),
        data.frame(sample = s, band = nb[b], measure = "FSI", t(setNames(c(round(sum(off) / ncol(M), 2), rep(NA, ncol(M) - 1)), cn)), check.names = FALSE))
}))))
tab <- tab[order(match(tab$sample, names(samples)), match(tab$measure, c("TO", "FROM", "NET", "FSI")), match(tab$band, nb)), ]
print(tab, row.names = FALSE); write_table(tab, "tableOA4_three_bands", row.names = FALSE)
d <- do.call(rbind, lapply(1:3, function(b) data.frame(date = dates, band = nb[b], v = f3$TCI[, b + 1])))
d$band <- factor(d$band, nb)
save_fig(ggplot(d, aes(date, v, colour = band)) + geom_line(linewidth = 0.7) + scale_colour_manual(values = col_bands) + yrs(3) + theme_paper,
         "figOA13_fsi_three_bands", 7.55, 4.45)
d <- do.call(rbind, lapply(1:3, function(b) long_panel(f3$NET[, , b + 1], dates, nb[b]))); d$band <- factor(d$band, nb)
g <- ggplot(d, aes(date, v, colour = band)) + geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
  geom_line(linewidth = 0.5) + scale_colour_manual(values = col_bands) + facet_wrap(~country, ncol = 3) + yrs(4) + theme_paper
save_fig(g, "figOA14_net_three_bands", 7.85, 4.45)

# excluding the United States (Figures OA.15-OA.16)
xs <- x[, setdiff(cn, "USA")]
fs <- conn_freq(tvpvar(xs, p))
for (b in 2:3) {
  d <- rbind(data.frame(date = dates, network = "With USA", v = base$freq$TCI[, b]), data.frame(date = dates, network = "Without USA", v = fs$TCI[, b]))
  g <- ggplot(d, aes(date, v, colour = network)) + geom_line(linewidth = 0.7) + scale_colour_manual(values = c("With USA" = "blue", "Without USA" = "red")) + yrs(3) + theme_paper
  save_fig(g, c("", "figOA15_fsi_short_without_usa", "figOA16_fsi_long_without_usa")[b], 7.55, 4.45)
}

# China in place of the United States (Table OA.5, Figure OA.17)
xc <- read_cpi("cpi_china")
fc <- conn_freq(tvpvar(xc, bic_lag(xc)))
for (b in 2:3) { tab <- avg_table(fc$CT[, , , b]); print(tab); write_table(tab, c("", "tableOA5_china_short", "tableOA5_china_long")[b]) }
fig_net(fc$NET[, , 2], fc$NET[, , 3], "figOA17_net_freq_china")

# rolling-window VAR (Figures OA.18-OA.19)
rw <- ConnectednessApproach(x, nlag = p, nfore = H_freq, window.size = prior_size + p, model = "VAR", connectedness = "Frequency",
                            corrected = FALSE, Connectedness_config = list(FrequencyConnectedness = list(partition = bands2, generalized = TRUE, scenario = "ABS")))
dw <- as.Date(rownames(rw$TCI))
fig_fsi(rw$TCI[, 2], rw$TCI[, 3], "figOA18_fsi_rolling", dw)
fig_net(rw$NET[, , 2], rw$NET[, , 3], "figOA19_net_freq_rolling", dw)
