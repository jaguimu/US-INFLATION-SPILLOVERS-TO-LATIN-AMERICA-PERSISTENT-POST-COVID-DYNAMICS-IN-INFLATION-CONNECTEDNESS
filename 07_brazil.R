# Section 5.4: Brazil's NET spillover. Figure A8, Table 6, Figure A9
source("R/setup.R")
source("R/breaks.R")
tb <- new.env(); sys.source("R/functions_tvpvarx.R", envir = tb)

base <- readRDS(file.path(out_dir, "estimates", "baseline.rds"))
t0 <- base$time; f0 <- base$freq; cn <- colnames(base$x)
dates <- as.Date(rownames(t0$TCI)); n <- length(dates); fm <- function(i) format(dates[i], "%YM%m")
files <- list.files(file.path(out_dir, "bootstrap", "replications"), full.names = TRUE)
r0 <- tb$juntar_replicas(unlist(lapply(files, readRDS), recursive = FALSE))$M0

# breaks and Table 6 (time domain and long run)
for (dom in c("time", "long run")) {
  y <- if (dom == "time") t0$NET[, "Brazil"] else f0$NET[, "Brazil", 3]
  D <- if (dom == "time") r0$NET[, "Brazil", ] else r0$NET_f[, "Brazil", 3, ]
  bp <- bai_perron(y)
  e <- bp$b[which(head(bp$means, -1) < 0 & bp$means[-1] > 0)[1]]
  j <- sort(unique(c(seq(0, min(120, n - e), 12), n - e)))
  pr <- round(100 * gn_prob(D, e, j), 1)
  cat(sprintf("\nBrazil NET (%s): breaks %s | means %s | event %s\n", dom, paste(fm(bp$b), collapse = ", "),
              paste(sprintf("%.2f", bp$means), collapse = " / "), fm(e)))
  print(setNames(pr, paste0("e+", j)))
  if (dom == "time") {
    write_table(data.frame(j = j, date = fm(e + j), probability = pr), "table6_gn_brazil", row.names = FALSE)
    g <- bp_plot(list(Brazil = y), list(Brazil = bp), dates, ncol = 1, step = "2 years") +
      geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) + theme(strip.text = element_blank())
    save_fig(g, "figA8_breaks_brazil", 8.15, 4)
    regimes <- cut(1:n, c(0, bp$b, n), labels = FALSE)
  }
}

# pairwise NET of Brazil: share of j's FEV due to Brazil minus share of Brazil's FEV due to j
CT <- t0$CT * 100
pw <- sapply(setdiff(cn, "Brazil"), function(j) CT[j, "Brazil", ] - CT["Brazil", j, ])
stopifnot(max(abs(rowSums(pw) - t0$NET[, "Brazil"])) < 1e-8)
avg <- round(rbind(apply(pw, 2, function(v) tapply(v, regimes, mean))), 2)
rownames(avg) <- tapply(fm(1:n), regimes, function(v) paste(v[1], v[length(v)], sep = "-"))
print(avg)
write_table(avg, "brazil_pairwise_regimes")
d <- long_panel(pw, dates)
g <- ggplot() + geom_col(data = d, aes(date, v, fill = country), width = 28, position = "stack") +
  geom_line(data = data.frame(date = dates, v = t0$NET[, "Brazil"]), aes(date, v), linewidth = 0.6) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y", expand = c(0.01, 0)) + theme_paper
save_fig(g, "figA9_pairwise_brazil")
