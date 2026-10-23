# Figures 2-9, 11, 12, A6, A7 with +/- 1 s.d. bootstrap bands; NET intervals in Tables 3 and 4
source("R/setup.R")
tb <- new.env(); sys.source("R/functions_tvpvarx.R", envir = tb)

base <- readRDS(file.path(out_dir, "estimates", "baseline.rds"))
ctrl <- readRDS(file.path(out_dir, "estimates", "controls.rds"))
files <- list.files(file.path(out_dir, "bootstrap", "replications"), full.names = TRUE)
R <- tb$juntar_replicas(unlist(lapply(files, readRDS), recursive = FALSE))
dates <- as.Date(rownames(base$time$TCI)); cn <- colnames(base$x); N <- length(cn)
cat("replications:", R$B, "\n")

band <- function(x, D) { s <- apply(D, 1, sd); x <- as.numeric(x); data.frame(est = x, lo = x - s, hi = x + s) }
yrs <- function(k) scale_x_date(date_breaks = paste(k, "years"), date_labels = "%Y", expand = c(0.01, 0))

fig_index <- function(x, D, name) {
  d <- cbind(date = dates, band(x, D))
  g <- ggplot(d, aes(date)) + geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey65", alpha = 0.6) +
    geom_line(aes(y = est), linewidth = 0.7) + yrs(2) + theme_paper
  save_fig(g, name)
}
fig_country <- function(X, D, name) {
  d <- do.call(rbind, lapply(1:N, function(i) cbind(date = dates, country = cn[i], band(X[, i], D[, i, ]))))
  d$country <- factor(d$country, cn)
  g <- ggplot(d, aes(date)) + geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey65", alpha = 0.6) +
    geom_line(aes(y = est), linewidth = 0.5) + facet_wrap(~country, ncol = 3) + yrs(4) + theme_paper
  save_fig(g, name)
}
fig_index_freq <- function(X, D, name) {
  d <- rbind(cbind(date = dates, band = "Short run", band(X[, 2], D[, 2, ])),
             cbind(date = dates, band = "Long run", band(X[, 3], D[, 3, ])))
  d$band <- factor(d$band, c("Long run", "Short run"))
  g <- ggplot(d, aes(date, colour = band, fill = band)) + geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, colour = NA) +
    geom_line(aes(y = est), linewidth = 0.8) + scale_colour_manual(values = col_bands) + scale_fill_manual(values = col_bands) +
    yrs(3) + theme_paper
  save_fig(g, name, 7.55, 4.45)
}
fig_country_freq <- function(X, D, name) {
  d <- do.call(rbind, lapply(1:N, function(i) rbind(
    cbind(date = dates, country = cn[i], band = "Short run", band(X[, i, 2], D[, i, 2, ])),
    cbind(date = dates, country = cn[i], band = "Long run", band(X[, i, 3], D[, i, 3, ])))))
  d$country <- factor(d$country, cn); d$band <- factor(d$band, c("Short run", "Long run"))
  g <- ggplot(d, aes(date, colour = band, fill = band)) + geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, colour = NA) + geom_line(aes(y = est), linewidth = 0.6) +
    scale_colour_manual(values = col_bands) + scale_fill_manual(values = col_bands) + facet_wrap(~country, ncol = 3) + yrs(4) + theme_paper
  save_fig(g, name, 7.85, 4.45)
}

# baseline
r0 <- R$M0; t0 <- base$time; f0 <- base$freq
fig_index(t0$TCI, r0$TCI, "fig02_tsi")
fig_country(t0$TO, r0$TO, "fig03_to")
fig_country(t0$FROM, r0$FROM, "fig04_from")
fig_country(t0$NET, r0$NET, "fig05_net")
fig_index_freq(f0$TCI, r0$TCI_f, "fig06_fsi")
fig_country_freq(f0$TO, r0$TO_f, "fig07_to_freq")
fig_country_freq(f0$FROM, r0$FROM_f, "fig08_from_freq")
fig_country_freq(f0$NET, r0$NET_f, "fig09_net_freq")

# controls
r1 <- R$M1; t1 <- ctrl$time; f1 <- ctrl$freq
fig_index_freq(f1$TCI, r1$TCI_f, "fig11_fsi_controls")
fig_country_freq(f1$NET, r1$NET_f, "fig12_net_freq_controls")
fig_index(t1$TCI, r1$TCI, "figA6_tsi_controls")
fig_country(t1$NET, r1$NET, "figA7_net_controls")

# NET +/- 1 s.d. of the time-averaged NET (Tables 3 and 4)
ci <- function(X, D) { m <- colMeans(X); s <- apply(apply(D, c(2, 3), mean), 1, sd); rbind(NET = m, lower = m - s, upper = m + s) }
net_ci <- rbind(cbind(measure = "time", round(ci(t0$NET, r0$NET), 2)),
                cbind(measure = "short run", round(ci(f0$NET[, , 2], r0$NET_f[, , 2, ]), 2)),
                cbind(measure = "long run", round(ci(f0$NET[, , 3], r0$NET_f[, , 3, ]), 2)))
colnames(net_ci) <- c("measure", cn)
print(net_ci, quote = FALSE)
write_table(net_ci, "tables3_4_net_intervals")
