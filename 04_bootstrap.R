# Recursive-design wild bootstrap (Goncalves and Kilian, 2004) for the baseline (M0) and controls (M1) models.
# Replications are saved in batches; rerunning continues from the last saved batch.
source("R/setup.R")
library(parallel)
tb <- new.env(); sys.source("R/functions_tvpvarx.R", envir = tb)

B <- 1000
batch_size <- 50
n_cores <- max(1, detectCores() - 1)
seed <- 20251001
boot_dir <- file.path(out_dir, "bootstrap")
dir.create(file.path(boot_dir, "replications"), recursive = TRUE, showWarnings = FALSE)

x <- read_cpi()
z <- read_controls()
p0 <- bic_lag(x)
l1 <- tb$seleccion_rezagos_X(x, z, p_max = p_max, s_max = 3, verbose = FALSE)$elegido$BIC
cfg <- tb$config_bootstrap(nlag0 = p0, nlag1 = l1[["nlag"]], zlag1 = l1[["zlag"]], nfore = H, nfore_frec = H_freq,
                           window_size = prior_size, kappa1 = kappa[1], kappa2 = kappa[2], particion = bands2,
                           modelos = c("M0", "M1"))

est <- tb$estimar_modelos(x, z, cfg, guardar_kalman = TRUE)
saveRDS(list(est = est, cfg = cfg, seed = seed), file.path(boot_dir, "estimation.rds"))

one <- function(b) tryCatch(tb$replica_bootstrap(b, est, x, z, cfg, semilla = seed), error = function(e) NULL)
batches <- split(1:B, ceiling((1:B) / batch_size))
batch_file <- function(i) file.path(boot_dir, "replications", sprintf("batch_%03d.rds", i))
todo <- which(!file.exists(sapply(seq_along(batches), batch_file)))

cl <- NULL
if (length(todo) > 0 && n_cores > 1) {
  cl <- makeCluster(n_cores)
  wd <- getwd()
  clusterExport(cl, c("wd", "est", "x", "z", "cfg", "seed"))
  invisible(clusterEvalQ(cl, { setwd(wd); library(zoo); tb <- new.env(); sys.source("R/functions_tvpvarx.R", envir = tb); NULL }))
}
for (i in todo) {
  reps <- if (is.null(cl)) lapply(batches[[i]], one) else parLapply(cl, batches[[i]], one)
  saveRDS(Filter(Negate(is.null), reps), batch_file(i))
  cat(sprintf("batch %d/%d done (%s)\n", i, length(batches), format(Sys.time(), "%H:%M")))
}
if (!is.null(cl)) stopCluster(cl)
