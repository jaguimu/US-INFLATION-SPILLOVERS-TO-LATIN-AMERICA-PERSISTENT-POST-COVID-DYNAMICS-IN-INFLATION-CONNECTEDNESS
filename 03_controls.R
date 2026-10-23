# TVP-VARX with energy prices and the broad USD as exogenous controls (Section 5.5)
source("R/setup.R")
tb <- new.env(); sys.source("R/functions_tvpvarx.R", envir = tb)

x <- read_cpi()
z <- read_controls()
lags <- tb$seleccion_rezagos_X(x, z, p_max = p_max, s_max = 3, verbose = FALSE)$elegido$BIC
vc <- list(TVPVAR = list(kappa1 = kappa[1], kappa2 = kappa[2], prior = "BayesPrior", gamma = 0.01))
tc <- tb$ConnectednessApproachX(x, z, nlag = lags[["nlag"]], zlag = lags[["zlag"]], nfore = H, window.size = prior_size,
                                corrected = FALSE, connectedness = "Time", VAR_config = vc, verbose = FALSE)
fc <- tb$FrequencyConnectedness(Phi = tc$B_t, Sigma = tc$Q_t, nfore = H_freq, partition = bands2,
                                generalized = TRUE, scenario = "ABS", corrected = FALSE, verbose = FALSE)
saveRDS(list(lags = lags, time = tc, freq = fc), file.path(out_dir, "estimates", "controls.rds"))

base <- readRDS(file.path(out_dir, "estimates", "baseline.rds"))
res <- rbind(baseline = c(TSI = mean(base$time$TCI), FSI_short = mean(base$freq$TCI[, 2]), FSI_long = mean(base$freq$TCI[, 3]),
                          NET_USA = mean(base$time$NET[, "USA"]), NET_USA_short = mean(base$freq$NET[, "USA", 2]), NET_USA_long = mean(base$freq$NET[, "USA", 3])),
             controls = c(mean(tc$TCI), mean(fc$TCI[, 2]), mean(fc$TCI[, 3]),
                          mean(tc$NET[, "USA"]), mean(fc$NET[, "USA", 2]), mean(fc$NET[, "USA", 3])))
print(round(res, 2))
write_table(round(res, 2), "controls_averages")
