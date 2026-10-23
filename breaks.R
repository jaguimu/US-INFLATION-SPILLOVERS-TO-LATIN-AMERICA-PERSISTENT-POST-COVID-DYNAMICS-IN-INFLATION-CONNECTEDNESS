# Bai-Perron (1998, 2003): number of breaks by UDmax and sequential sup F(l+1|l) at 5%
library(mbreaks)
bai_perron <- function(y, trend = FALSE, eps1 = 0.15, m = 5, signif = 2) {
  n <- length(y); df <- data.frame(y = as.numeric(y), tend = (1:n) / 12)
  zn <- if (trend) "tend" else NULL
  invisible(capture.output(r <- suppressWarnings(mdl("y", z_name = zn, data = df, eps1 = eps1, m = m,
                                                     prewhit = 1, robust = 1, hetvar = 1, signif = signif))))
  ud <- as.numeric(r$sbtests$UDmax[1, 1]); cv <- as.numeric(r$sbtests$cv_Dmax[, 1])
  sq <- as.numeric(r$seqtests$supfl); cvs <- r$seqtests$cv
  nb <- 0
  if (ud > cv[signif]) { nb <- length(sq); for (l in seq_len(nb - 1)) if (!(sq[l + 1] > cvs[signif, l + 1])) { nb <- l; break } }
  if (nb == 0) return(list(ud = ud, cv1 = cv[4], nb = 0, b = integer(0), fit = rep(mean(df$y), n), means = mean(df$y)))
  invisible(capture.output(f <- suppressWarnings(dofix("y", z_name = zn, data = df, fixn = nb, eps1 = eps1,
                                                       prewhit = 1, robust = 1, hetvar = 1))))
  b <- as.integer(f$date); reg <- findInterval(1:n, b + 1)
  list(ud = ud, cv1 = cv[4], nb = nb, b = b, reg = reg, fit = as.numeric(f$fitted.values),
       means = as.numeric(tapply(df$y, reg, mean)))
}

# Greenwood-Nimmo et al. (2024), eq. (5): Pr(S_{e+j} > mean_b S_{e-1})
gn_prob <- function(D, e, j, ref = 1) {
  Sbar <- mean(D[(e - ref):(e - 1), ])
  sapply(j, function(h) mean(D[e + h, ] - Sbar > 0))
}

bp_plot <- function(series, res, dates, ncol = 3, step = "4 years") {
  d <- do.call(rbind, lapply(names(series), function(s) data.frame(date = dates, series = s, y = as.numeric(series[[s]]),
                                                                    fit = res[[s]]$fit, reg = findInterval(seq_along(dates), res[[s]]$b + 1))))
  v <- do.call(rbind, lapply(names(series), function(s) if (res[[s]]$nb > 0) data.frame(series = s, x = dates[res[[s]]$b])))
  d$series <- factor(d$series, names(series)); v$series <- factor(v$series, names(series))
  ggplot(d, aes(date)) + geom_vline(data = v, aes(xintercept = x), linetype = "dashed", colour = "grey30") +
    geom_line(aes(y = y), linewidth = 0.5) + geom_line(aes(y = fit, group = reg), colour = "red", linewidth = 0.9) +
    facet_wrap(~series, ncol = ncol, scales = "free_y") +
    scale_x_date(date_breaks = step, date_labels = "%Y", expand = c(0.01, 0)) +
    theme_minimal(base_size = 11) + theme(strip.text = element_text(face = "bold"), axis.title = element_blank())
}
