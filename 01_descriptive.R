# Table 1, Table 2 and Figure 1
source("R/setup.R")

x <- read_cpi()
x <- x[, c("USA", "Brazil", "Chile", "Colombia", "Mexico", "Peru")]
saveRDS(x, file.path(out_dir, "estimates", "inflation.rds"))

# Table 1
ss <- SummaryStatistics(x, portmanteau = "Ljung-Box", nlag = 1)
ss <- ss[seq_len(which(rownames(ss) == "Q2(1)")), ]
tab1 <- ss[!grepl("ERS", rownames(ss)) & rownames(ss) != "", ]
print(tab1, quote = FALSE)
write_table(tab1, "table1_descriptive")

# Table 2
cr <- cor(coredata(x)); n <- nrow(x)
pv <- 2 * pt(-abs(cr * sqrt((n - 2) / (1 - cr^2))), n - 2)
tab2 <- matrix(paste0(sprintf("%.3f", cr), ifelse(pv < 0.01, "***", ifelse(pv < 0.05, "**", ifelse(pv < 0.1, "*", "")))), ncol(x), dimnames = dimnames(cr))
diag(tab2) <- "1.000"; tab2[upper.tri(tab2)] <- ""
print(tab2, quote = FALSE)
write_table(tab2, "table2_correlations")

# Figure 1
g <- ggplot(long_panel(x, index(x)), aes(date, v)) + geom_line(linewidth = 0.4) +
  facet_wrap(~country, ncol = 3) + theme_paper
save_fig(g, "fig01_inflation")
