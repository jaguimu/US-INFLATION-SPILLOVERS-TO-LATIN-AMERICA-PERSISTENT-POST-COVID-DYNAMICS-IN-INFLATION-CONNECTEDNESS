# Section 5.3: Bai-Perron breaks (Figure 10, Figure A5) and bootstrap probabilities (Table 5)
source("R/setup.R")
source("R/breaks.R")
tb <- new.env(); sys.source("R/functions_tvpvarx.R", envir = tb)

base <- readRDS(file.path(out_dir, "estimates", "baseline.rds"))
t0 <- base$time; f0 <- base$freq
dates <- as.Date(rownames(t0$TCI)); fm <- function(i) format(dates[i], "%YM%m")
series <- list("TSI (time domain)" = t0$TCI[, 1], "FSI (short run)" = f0$TCI[, 2], "FSI (long run)" = f0$TCI[, 3],
               "NET USA (time domain)" = t0$NET[, "USA"], "NET USA (short run)" = f0$NET[, "USA", 2], "NET USA (long run)" = f0$NET[, "USA", 3])

bp_mean <- lapply(series, bai_perron)
bp_trend <- lapply(series, bai_perron, trend = TRUE)
tab <- do.call(rbind, lapply(names(series), function(s) data.frame(series = s, UDmax = round(bp_mean[[s]]$ud, 2),
  breaks_mean = paste(fm(bp_mean[[s]]$b), collapse = ", "), regime_means = paste(sprintf("%.2f", bp_mean[[s]]$means), collapse = " / "),
  breaks_mean_trend = paste(fm(bp_trend[[s]]$b), collapse = ", "))))
print(tab, right = FALSE, row.names = FALSE)
write_table(tab, "bai_perron", row.names = FALSE)
save_fig(bp_plot(series, bp_mean, dates), "fig10_breaks_mean", 10, 5.2)
save_fig(bp_plot(series, bp_trend, dates), "figA5_breaks_trend", 10, 5.2)

# Table 5
files <- list.files(file.path(out_dir, "bootstrap", "replications"), full.names = TRUE)
r0 <- tb$juntar_replicas(unlist(lapply(files, readRDS), recursive = FALSE))$M0
reps <- list(r0$TCI, r0$TCI_f[, 2, ], r0$TCI_f[, 3, ], r0$NET[, "USA", ], r0$NET_f[, "USA", 2, ], r0$NET_f[, "USA", 3, ])
e <- match(covid, dates); j <- c(0, 1, 3, 6, 12, 24, 36)
for (ref in c(1, 3)) {
  tab5 <- t(sapply(reps, function(D) round(100 * gn_prob(D, e, j, ref), 1)))
  dimnames(tab5) <- list(names(series), paste0("e+", j))
  cat("\nreference: average of", ref, "month(s) before", fm(e), "\n"); print(tab5)
  write_table(tab5, if (ref == 1) "table5_gn" else "table5_gn_ref3")
}
