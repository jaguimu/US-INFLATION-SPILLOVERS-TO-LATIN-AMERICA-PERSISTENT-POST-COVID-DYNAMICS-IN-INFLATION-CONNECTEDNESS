# Common settings, data and helpers. Run scripts from the codes_pack folder.

suppressPackageStartupMessages({
  library(ConnectednessApproach)
  library(zoo)
  library(openxlsx)
  library(seasonal)
  library(vars)
  library(ggplot2)
})

data_file  <- "data/data.xlsx"
out_dir    <- "output"
p_max      <- 6
H          <- 12
H_freq     <- 100
kappa      <- c(0.99, 0.96)
prior_size <- 200
bands2     <- c(pi + 0.00001, pi / 3, 0)
bands3     <- c(pi + 0.00001, pi / 4, pi / 10, 0)
covid      <- as.Date("2020-03-01")
for (d in c("estimates", "tables", "figures")) dir.create(file.path(out_dir, d), recursive = TRUE, showWarnings = FALSE)

# CPI -> X-13 -> 100 * dlog
read_cpi <- function(sheet = "cpi") {
  d <- read.xlsx(data_file, sheet = sheet, detectDates = TRUE)
  start <- c(as.integer(format(d$date[1], "%Y")), as.integer(format(d$date[1], "%m")))
  sa <- sapply(d[, -1], function(s) as.numeric(final(seas(ts(s, frequency = 12, start = start)))))
  x <- zoo(100 * diff(log(sa)), d$date[-1])
  colnames(x) <- names(d)[-1]
  x
}

# controls -> 100 * dlog (not seasonally adjusted)
read_controls <- function() {
  d <- read.xlsx(data_file, sheet = "controls", detectDates = TRUE)
  z <- zoo(100 * diff(log(as.matrix(d[, -1]))), d$date[-1])
  colnames(z) <- names(d)[-1]
  z
}

bic_lag <- function(x) as.integer(VARselect(x, lag.max = p_max, type = "const")$selection["SC(n)"])

tvpvar <- function(x, p, l = kappa, prior = c("bayes", "minnesota")) {
  prior <- match.arg(prior)
  pr <- if (prior == "bayes") BayesPrior(x = x, size = prior_size - p, nlag = p) else MinnesotaPrior(gamma = 0.1, k = ncol(x), nlag = p)
  TVPVAR(x, configuration = list(l = l, nlag = p, prior = pr))
}
conn_time <- function(fit, h = H) TimeConnectedness(Phi = fit$B_t, Sigma = fit$Q_t, nfore = h, generalized = TRUE, corrected = FALSE)
conn_freq <- function(fit, partition = bands2) FrequencyConnectedness(Phi = fit$B_t, Sigma = fit$Q_t, nfore = H_freq, partition = partition,
                                                                     generalized = TRUE, scenario = "ABS", corrected = FALSE)

# average connectedness table (pp)
avg_table <- function(CT) {
  M <- apply(CT, 1:2, mean) * 100
  off <- M; diag(off) <- 0
  tab <- cbind(round(M, 2), FROM = round(rowSums(off), 2))
  rbind(tab, TO = c(round(colSums(off), 2), NA), NET = c(round(colSums(off) - rowSums(off), 2), round(sum(off) / ncol(M), 2)))
}
write_table <- function(x, name, row.names = TRUE) write.csv(x, file.path(out_dir, "tables", paste0(name, ".csv")), row.names = row.names, na = "")

# plots
theme_paper <- theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"), legend.position = "bottom",
        legend.title = element_blank(), axis.title = element_blank())
col_bands <- c("Short run" = "red", "Medium run" = "black", "Long run" = "blue")
save_fig <- function(g, name, w = 8.15, h = 4.55) ggsave(file.path(out_dir, "figures", paste0(name, ".png")), g, width = w, height = h, dpi = 150, bg = "white")
long_panel <- function(M, dates, band = NULL) {
  d <- do.call(rbind, lapply(colnames(M), function(k) data.frame(date = dates, country = k, v = as.numeric(M[, k]))))
  d$country <- factor(d$country, colnames(M)); if (!is.null(band)) d$band <- band
  d
}
