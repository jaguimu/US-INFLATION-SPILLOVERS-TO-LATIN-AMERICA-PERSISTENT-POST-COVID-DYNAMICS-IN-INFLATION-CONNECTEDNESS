# Baseline TVP-VAR: time (H = 12) and frequency (H = 100) connectedness. Tables 3 and 4
source("R/setup.R")

x <- read_cpi()
p <- bic_lag(x)
fit <- tvpvar(x, p)
tc <- conn_time(fit)
fc <- conn_freq(fit)
saveRDS(list(x = x, p = p, fit = fit, time = tc, freq = fc), file.path(out_dir, "estimates", "baseline.rds"))

tab3 <- avg_table(tc$CT)
tab4 <- lapply(2:3, function(b) avg_table(fc$CT[, , , b]))
names(tab4) <- c("short", "long")
print(tab3); print(tab4)
write_table(tab3, "table3_time")
write_table(tab4$short, "table4_short_run")
write_table(tab4$long, "table4_long_run")
cat(sprintf("p = %d | TSI = %.2f | FSI short = %.2f | FSI long = %.2f\n", p, mean(tc$TCI), mean(fc$TCI[, 2]), mean(fc$TCI[, 3])))
