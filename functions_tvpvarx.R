# TVP-VAR / TVP-VARX with exogenous controls and recursive wild bootstrap.
# TVPVAR, BayesPrior and connectedness functions follow ConnectednessApproach 1.0.4.

suppressPackageStartupMessages({
  library(zoo)
  library(MASS)
})

UninformativePrior <- function(k, nlag) {
  if (k <= 1) stop("k must be an integer larger than 1")
  if (nlag <= 0) stop("nlag must be a positive integer")
  m <- nlag * (k^2)
  bprior <- c(cbind(0 * diag(k), matrix(0, ncol = (nlag - 1) * k, nrow = k)))
  V_i <- matrix(4, nrow = (m / k), ncol = k)
  Vprior <- diag(c(t(V_i)))
  list(bprior = bprior, Vprior = Vprior)
}

MinnesotaPrior <- function(gamma = 0.1, k, nlag) {
  if (k <= 1) stop("k must be an integer larger than 1")
  if (nlag <= 0) stop("nlag must be a positive integer")
  m <- nlag * (k^2)
  bprior <- c(cbind(0 * diag(k), matrix(0, ncol = (nlag - 1) * k, nrow = k)))
  V_i <- matrix(0, nrow = (m / k), ncol = k)
  for (i in 1:k) {
    for (j in 1:(m / k)) {
      V_i[j, i] <- gamma / (ceiling(j / k)^2)
    }
  }
  Vprior <- diag(c(t(V_i)))
  list(bprior = bprior, Vprior = Vprior)
}

BayesPrior <- function(x, size = NULL, nlag) {
  if (!inherits(x, "zoo")) stop("Data must be a 'zoo' object")
  if (nlag <= 0) stop("nlag must be a positive integer")
  if (is.null(size)) size <- nrow(x)
  if ((size - nlag) <= 0) stop("size must be larger than nlag")

  p <- ncol(x)
  size <- size - nlag
  yt <- t(x[(nlag + 1):(size + nlag), ])
  m <- p + nlag * (p^2)

  Zt <- NULL
  for (i in (nlag + 1):(size + nlag)) {
    ztemp <- diag(p)
    for (j in 1:nlag) {
      xlag <- x[(i - j), 1:p]
      xtemp <- matrix(0, p, p * p)
      for (jj in 1:p) {
        xtemp[jj, ((jj - 1) * p + 1):(jj * p)] <- xlag
      }
      ztemp <- cbind(ztemp, xtemp)
    }
    Zt <- rbind(Zt, ztemp)
  }

  vbar <- matrix(0, m, m)
  xhy <- matrix(0, m, 1)
  for (i in 1:size) {
    zhat1 <- Zt[((i - 1) * p + 1):(i * p), ]
    vbar <- vbar + t(zhat1) %*% zhat1
    xhy <- xhy + t(zhat1) %*% as.matrix(yt[, i])
  }
  vbar <- ginv(vbar)
  aols <- vbar %*% xhy

  sse2 <- matrix(0, p, p)
  for (i in 1:size) {
    zhat1 <- Zt[((i - 1) * p + 1):(i * p), ]
    sse2 <- sse2 + (yt[, i] - zhat1 %*% aols) %*% t(yt[, i] - zhat1 %*% aols)
  }
  hbar <- sse2 / size

  list(bprior = aols[-c(1:p)], Vprior = vbar[-c(1:p), -c(1:p)], Q = hbar)
}

TVPVAR <- function(x, configuration = list(l = c(0.99, 0.99), nlag = 1, prior = NULL),
                   verbose = TRUE, guardar_kalman = FALSE) {
  l <- as.numeric(configuration$l)
  nlag <- configuration$nlag
  prior <- configuration$prior
  if (nlag <= 0) stop("nlag must be a positive integer")
  if (!inherits(x, "zoo")) stop("Data must be a 'zoo' object")
  if (l[1] <= 0 || l[1] >= 1) stop("kappa1 must be in (0, 1)")
  if (l[2] <= 0 || l[2] >= 1) stop("kappa2 must be in (0, 1)")

  k <- ncol(x)
  if (is.null(prior)) prior <- UninformativePrior(k, nlag)
  NAMES <- colnames(x)
  if (is.null(NAMES)) NAMES <- 1:k

  beta_0.mean <- prior$bprior
  beta_0.var <- prior$Vprior
  Q <- prior$Q
  if (is.null(Q)) Q <- cov(x)

  create_RHS_NI <- function(templag, r, nlag, t) {
    K <- nlag * (r^2)
    x_t <- matrix(0, (t - nlag) * r, K)
    for (i in 1:(t - nlag)) {
      ztemp <- NULL
      for (j in 1:nlag) {
        xtemp <- templag[i, ((j - 1) * r + 1):(j * r)]
        xtemp <- t(kronecker(diag(r), xtemp))
        ztemp <- cbind(ztemp, xtemp)
      }
      x_t[((i - 1) * r + 1):(i * r), ] <- ztemp
    }
    list(x_t = x_t, K = K)
  }

  media_x <- colMeans(x)
  x <- scale(x, TRUE, FALSE)
  r <- n <- ncol(x)
  m <- nlag * (r^2)
  k <- nlag * r
  t <- nrow(x)

  beta_pred <- matrix(0, m, t)
  beta_update <- matrix(0, m, t)
  Rb_t <- array(0, c(m, m, t))
  Sb_t <- array(0, c(m, m, t))
  beta_t <- array(0, c(k, k, t))
  Q_t <- array(0, c(r, r, t), dimnames = list(NAMES, NAMES, as.character(zoo::index(x))))

  l_2 <- l[2]
  l_4 <- l[1]

  yy <- x[(nlag + 1):t, ]
  xx <- embed(x, nlag + 1)[, -c(1:r)]
  templag <- embed(x, nlag + 1)[, -c(1:r)]
  Flagtemp <- create_RHS_NI(templag, r, nlag, t)$x_t
  Flag <- rbind(matrix(0, k, m), Flagtemp)

  if (guardar_kalman) {
    K_t <- array(0, c(m, r, t)); e_t <- matrix(0, r, t)
  }
  if (verbose) pb <- txtProgressBar(min = 0, max = t, style = 3)
  for (irep in 1:t) {

    if (irep == 1) {
      Q_t[, , irep] <- Q
    } else {
      if (irep <= (nlag + 1)) {
        Gf_t <- 0.1 * (t(matrix(x[irep, ], nrow = 1)) %*% (x[irep, ]))
      } else {
        res <- yy[(irep - nlag), ] - xx[(irep - nlag), ] %*% t(B[1:r, 1:k])
        Gf_t <- t(res) %*% res
      }
      Q_t[, , irep] <- l_2 * Q_t[, , (irep - 1)] + (1 - l_2) * Gf_t[1:r, 1:r]
    }

    if (irep <= (nlag + 1)) {
      beta_pred[, irep] <- beta_0.mean
      beta_update[, irep] <- beta_pred[, irep]
      Rb_t[, , irep] <- beta_0.var
    } else {
      beta_pred[, irep] <- beta_update[, (irep - 1)]
      Rb_t[, , irep] <- (1 / l_4) * Sb_t[, , (irep - 1)]
    }

    if (irep >= (nlag + 1)) {
      H <- Flag[((irep - 1) * r + 1):(irep * r), ]
      Rx <- Rb_t[, , irep] %*% t(H)
      KV_b <- Q_t[, , irep] + H %*% Rx
      KG <- Rx %*% ginv(KV_b)
      beta_update[, irep] <- matrix(beta_pred[, irep], ncol = 1) +
        KG %*% (t(matrix(x[irep, ], nrow = 1)) - H %*% matrix(beta_pred[, irep], ncol = 1))
      Sb_t[, , irep] <- Rb_t[, , irep] - KG %*% (H %*% Rb_t[, , irep])
      if (guardar_kalman) {
        K_t[, , irep] <- KG
        e_t[, irep] <- as.numeric(x[irep, ]) - as.numeric(H %*% beta_pred[, irep])
      }
    }

    bb <- matrix(beta_update[, irep], ncol = 1)
    splace <- 0
    biga <- matrix(0, r, r * nlag)
    for (ii in 1:nlag) {
      for (iii in 1:r) {
        biga[iii, ((ii - 1) * r + 1):(ii * r)] <- t(bb[(splace + 1):(splace + r), 1])
        splace <- splace + r
      }
    }
    B <- rbind(biga, cbind(diag(r * (nlag - 1)), matrix(0, nrow = r * (nlag - 1), ncol = r)))

    if ((max(abs(eigen(B)$values)) <= 1) || (irep == 1)) {
      beta_t[, , irep] <- B
    } else {
      beta_t[, , irep] <- beta_t[, , (irep - 1)]
      beta_update[, irep] <- 0.99 * beta_update[, (irep - 1)]
    }
    if (verbose) setTxtProgressBar(pb, irep)
  }
  if (verbose) close(pb)

  B_t <- beta_t[1:ncol(Q_t), , ]
  out <- list(B_t = B_t, Q_t = Q_t)
  if (guardar_kalman) {

    out$kalman <- list(beta_pred = beta_pred, K = K_t, e = e_t, xc = matrix(as.numeric(x), nrow = t),
                       media = media_x, nlag = nlag, W = NULL, fechas = zoo::index(x),
                       nombres = NAMES)
  }
  out
}

IRF <- function(Phi, Sigma, nfore = 10, orth = TRUE) {
  if (nfore <= 0) stop("nfore must be a positive integer")
  if (length(dim(Sigma)) > 2) Sigma <- Sigma[, , 1]

  k <- dim(Phi)[1]
  p <- floor(dim(Phi)[2] / k)
  if (orth) {
    m1 <- eigen(Sigma)
    Sh <- m1$vectors %*% diag(sqrt(m1$values)) %*% t(m1$vectors)
  }

  PSI <- diag(rep(1, k))
  WGT <- if (orth) c(PSI %*% Sh) else c(PSI)
  for (il in 1:nfore) {
    tmp <- matrix(0, k, k)
    iend <- min(il, p)
    for (j in 1:iend) {
      jdx <- (il - j)
      kdx <- j * k
      tmp <- tmp + Phi[, (kdx - k + 1):kdx] %*% PSI[, (jdx * k + 1):(jdx * k + k)]
    }
    PSI <- cbind(PSI, tmp)
    WGT <- cbind(WGT, if (orth) c(tmp %*% Sh) else c(tmp))
  }

  irf <- list()
  for (i in 1:k) {
    irf[[i]] <- t(WGT[((i - 1) * k + 1):(i * k), ])
  }
  list(irf = irf)
}

FEVD <- function(Phi, Sigma, nfore = 100, type = c("time", "frequency"),
                 generalized = TRUE, range = NULL) {
  type <- match.arg(type)
  if (nfore <= 0) stop("nfore must be a positive integer")
  if (length(dim(Sigma)) > 2) Sigma <- Sigma[, , 1]
  k <- ncol(Sigma)
  irf <- IRF(Phi = Phi, Sigma = Sigma, nfore = nfore, orth = FALSE)$irf

  if (type == "time") {

    Psi <- lapply(1:nfore, function(j) sapply(irf, function(i) i[j, ]))
    denom <- diag(Reduce("+", lapply(Psi, function(i) i %*% Sigma %*% t(i))))
    if (generalized) {

      irf_ <- lapply(Psi, function(i) t(i %*% Sigma %*% solve(diag(sqrt(diag(Sigma))))))
      IRFg <- array(unlist(irf_), c(k, k, nfore))
      enum <- apply(IRFg^2, 1:2, sum)
      nfevd <- t(enum) / c(apply(enum, 1, sum))
      FEVDm <- nfevd / apply(nfevd, 1, sum)
    } else {
      irf_ <- lapply(Psi, function(i) chol(Sigma) %*% t(i))
      IRFg <- array(unlist(irf_), c(k, k, nfore))
      fevd <- IRFg^2
      tab_ <- array(0, c(k, k, nfore))
      for (i in 1:nfore) tab_[, , i] <- t(fevd[, , i] %*% solve(diag(denom)))
      FEVDm <- apply(tab_, 1:2, sum)
    }

  } else {
    fftir <- lapply(irf, function(i) apply(i, 2, fft))
    fftir <- lapply(1:(nfore + 1), function(j) sapply(fftir, function(i) i[j, ]))
    denom <- diag(Re(Reduce("+", lapply(fftir, function(i)
      i %*% Sigma %*% t(Conj(i)) / nfore)[range])))
    if (generalized) {
      irf_ <- lapply(fftir, function(i) (abs(i %*% Sigma)))
      IRF_ <- array(unlist(irf_), c(k, k, nfore + 1))
      enum_ <- IRF_^2 / (nfore + 1)
      tab <- list()
      for (i in 1:dim(enum_)[3]) {
        fevd_ <- enum_[, , 1]
        for (j in 1:dim(enum_)[2]) {
          fevd_[j, ] <- enum_[j, , i] / (denom[j] * diag(Sigma))
        }
        tab[[i]] <- t(fevd_)
      }
      tot <- apply(Reduce("+", tab[range]), 2, sum)
      FEVDm <- lapply(tab, function(i) t(i) / tot)
      IRFg <- NULL
    } else {
      irf_ <- lapply(fftir, function(i) (abs(i %*% t(chol(Sigma)))))
      IRFg <- array(unlist(irf_), c(k, k, nfore + 1))
      enum <- IRFg^2 / (nfore + 1)
      FEVDm <- list()
      for (i in 1:dim(enum)[3]) {
        fevd_ <- enum[, , 1]
        for (j in 1:dim(enum)[2]) fevd_[j, ] <- enum[j, , i] / denom[j]
        FEVDm[[i]] <- fevd_
      }
    }
  }
  list(IRF = IRFg, FEVD = FEVDm)
}

ConnectednessTable <- function(FEVD, digit = 2) {
  if (length(dim(FEVD)) <= 1) stop("FEVD must be at least a matrix")
  NAMES <- colnames(FEVD)
  k <- dim(FEVD)[1]
  if (is.null(NAMES)) NAMES <- 1:k

  CT <- apply(FEVD, 1:2, mean) * 100
  OWN <- diag(diag(CT))
  TO <- colSums(CT - OWN)
  FROM <- rowSums(CT - OWN)
  NET <- TO - FROM
  TCI <- mean(TO)
  cTCI <- TCI * k / (k - 1)
  NPDC <- CT - t(CT)
  NPT <- rowSums(NPDC < 0)
  INFLUENCE <- 100 * abs(NPDC / t(t(CT) + CT))

  fmt <- function(v) format(round(v, digit), nsmall = digit)
  table <- fmt(cbind(CT, FROM))
  to <- c(fmt(c(TO, sum(TO))))
  inc <- c(fmt(colSums(CT)), "cTCI/TCI")
  tci <- paste0(fmt(cTCI), "/", fmt(TCI))
  net <- c(fmt(NET), tci)
  npt <- c(fmt(NPT), "")
  TABLE <- rbind(table, to, inc, net, npt)
  colnames(TABLE) <- c(NAMES, "FROM")
  rownames(TABLE) <- c(NAMES, "TO", "Inc.Own", "NET", "NPT")

  PCI <- matrix(NA, k, k)
  for (i in 1:k) {
    for (j in 1:k) {
      PCI[i, j] <- 200 * (CT[i, j] + CT[j, i]) / (CT[i, i] + CT[i, j] + CT[j, i] + CT[j, j])
    }
  }
  list(FEVD = CT, TCI = TCI, cTCI = cTCI, PCI = PCI, TO = TO, FROM = FROM,
       NET = NET, NPDC = NPDC, TABLE = TABLE, NPT = NPT, INFLUENCE = INFLUENCE)
}

TimeConnectedness <- function(Phi = NULL, Sigma = NULL, nfore = 10, generalized = TRUE,
                              corrected = FALSE, FEVD = NULL, verbose = TRUE) {
  if ((is.null(Phi) || is.null(Sigma)) && is.null(FEVD)) {
    stop("Provide Phi and Sigma, or FEVD")
  }
  if (is.null(FEVD)) {
    NAMES <- colnames(Sigma)
    k <- ncol(Sigma)
    if (length(dim(Phi)) == 2) Phi <- array(Phi, c(nrow(Phi), ncol(Phi), 1))
    if (length(dim(Sigma)) == 2) Sigma <- array(Sigma, c(nrow(Sigma), ncol(Sigma), 1))
    t <- dim(Sigma)[3]
    if (is.null(NAMES)) NAMES <- 1:k
    date <- as.character(as.Date(dimnames(Sigma)[[3]]))
    CT <- array(NA, c(k, k, t), dimnames = list(NAMES, NAMES, date))
    for (i in 1:t) {
      CT[, , i] <- get("FEVD", mode = "function")(Phi = Phi[, , i], Sigma = Sigma[, , i],
                                                 nfore = nfore, generalized = generalized,
                                                 type = "time")$FEVD
    }
  } else {
    CT <- FEVD
    t <- dim(CT)[3]
    NAMES <- dimnames(CT)[[1]]
    k <- length(NAMES)
    date <- as.character(as.Date(dimnames(CT)[[3]]))
  }

  TCI <- array(NA, c(t, 1), dimnames = list(date, "TCI"))
  NPT <- NET <- FROM <- TO <- array(NA, c(t, k), dimnames = list(date, NAMES))
  PCI <- NPDC <- INFLUENCE <- array(NA, c(k, k, t), dimnames = list(NAMES, NAMES, date))
  if (verbose) pb <- txtProgressBar(min = 0, max = t, style = 3)
  for (i in 1:t) {
    ct <- ConnectednessTable(CT[, , i])
    TO[i, ] <- ct$TO
    FROM[i, ] <- ct$FROM
    NET[i, ] <- ct$NET
    NPT[i, ] <- ct$NPT
    NPDC[, , i] <- ct$NPDC
    INFLUENCE[, , i] <- ct$INFLUENCE
    PCI[, , i] <- ct$PCI
    TCI[i, ] <- if (corrected) ct$cTCI else ct$TCI
    if (verbose) setTxtProgressBar(pb, i)
  }
  if (verbose) close(pb)

  TABLE <- ConnectednessTable(CT)$TABLE
  config <- list(nfore = nfore, approach = "Time", generalized = generalized, corrected = corrected)
  list(TABLE = TABLE, CT = CT, TCI = TCI, TO = TO, FROM = FROM, NET = NET, NPT = NPT,
       NPDC = NPDC, PCI = PCI, INFLUENCE = INFLUENCE, config = config)
}

getIndeces <- function(n.ahead, up, down) {
  space <- (0:floor(n.ahead / 2)) / ((floor(n.ahead / 2))) * pi
  output <- ((space >= down) & (space < up)) * 1
  if (n.ahead %% 2 == 0) {
    output <- c(output, rev(output[2:(length(output) - 1)]))
  } else {
    output <- c(output, rev(output[2:length(output)]))
  }
  which(output == 1)
}

getPartition <- function(partition, n.ahead) {
  if (!all(sort(partition, decreasing = TRUE) == partition)) {
    stop("Limits must be in decreasing order.")
  }
  part <- lapply(1:(length(partition) - 1), function(i)
    getIndeces(n.ahead + 1, partition[i], partition[i + 1]))
  if (length(unique(do.call(c, part))) != (n.ahead + 1)) {
    warning("Partition does not cover the full range.")
  }
  if (any(sapply(part, length) == 0)) {
    stop("Empty band for this nfore. Change partition or increase nfore.")
  }
  part
}

FrequencyConnectedness <- function(Phi, Sigma, nfore = 100, partition = c(pi, pi / 2, 0),
                                   generalized = TRUE, orth = FALSE, scenario = "ABS",
                                   corrected = FALSE, verbose = TRUE) {
  NAMES <- colnames(Sigma)
  if (length(dim(Phi)) == 2) Phi <- array(Phi, c(nrow(Phi), ncol(Phi), 1))
  if (length(dim(Sigma)) == 2) Sigma <- array(Sigma, c(nrow(Sigma), ncol(Sigma), 1))
  k <- dim(Sigma)[1]
  t <- dim(Sigma)[3]
  if (is.null(NAMES)) NAMES <- 1:k

  periods <- round(pi / partition)
  period_names <- NULL
  for (i in 1:(length(periods) - 1)) {
    period_names <- c(period_names, paste0(periods[i], "-", periods[i + 1]))
  }
  period_names <- c("Total", period_names)
  date <- as.character(dimnames(Sigma)[[3]])
  interval <- length(period_names)
  new_p <- getPartition(partition, nfore)
  range <- sort(unique(do.call(c, new_p)))

  TCI <- array(0, c(t, interval), dimnames = list(date, period_names))
  CT <- NPDC <- array(0, c(k, k, t, interval), dimnames = list(NAMES, NAMES, date, period_names))
  NET <- FROM <- TO <- NPT <- array(0, c(t, k, interval), dimnames = list(date, NAMES, period_names))
  PCI <- INFLUENCE <- array(0, c(k, k, t, interval), dimnames = list(NAMES, NAMES, date, period_names))

  if (verbose) pb <- txtProgressBar(min = 0, max = t, style = 3)
  for (i in 1:t) {
    decomp <- FEVD(Phi = Phi[, , i], Sigma = Sigma[, , i], nfore = nfore,
                   generalized = generalized, type = "frequency", range = range)$FEVD
    for (ij in 1:length(decomp)) {
      rownames(decomp[[ij]]) <- colnames(decomp[[ij]]) <- 1:ncol(Sigma)
    }
    tables <- lapply(new_p, function(j) Reduce("+", decomp[j]))
    for (j in 2:interval) {
      tab_j <- if (scenario == "ABS") tables[[j - 1]] else tables[[j - 1]] / sum(sum(tables[[j - 1]])) * k
      dca <- ConnectednessTable(tab_j)
      CT[, , i, j] <- dca$FEVD
      TO[i, , j] <- dca$TO
      FROM[i, , j] <- dca$FROM
      NET[i, , j] <- dca$NET
      NPDC[, , i, j] <- dca$NPDC
      INFLUENCE[, , i, j] <- dca$INFLUENCE
      NPT[i, , j] <- dca$NPT
      TCI[i, j] <- if (corrected) dca$cTCI else dca$TCI
    }
    if (verbose) setTxtProgressBar(pb, i)
  }
  if (verbose) close(pb)

  CT[, , , 1] <- apply(CT, 1:3, sum)
  TCI[, 1] <- apply(TCI, 1, sum)
  TO[, , 1] <- apply(TO, 1:2, sum)
  FROM[, , 1] <- apply(FROM, 1:2, sum)
  NET[, , 1] <- apply(NET, 1:2, sum)
  NPDC[, , , 1] <- apply(NPDC, 1:3, sum)
  for (ij in 1:t) {
    for (jl in interval:1) {
      for (i in 1:k) {
        for (j in 1:k) {
          PCI[i, j, ij, jl] <- 200 * (CT[i, j, ij, jl] + CT[j, i, ij, jl]) /
            (CT[i, i, ij, 1] + CT[i, j, ij, 1] + CT[j, i, ij, 1] + CT[j, j, ij, 1])
        }
      }
      INFLUENCE[, , ij, jl] <- 100 * abs(NPDC[, , ij, jl] / t(t(CT[, , ij, 1]) + CT[, , ij, 1]))
    }
    NPT[ij, , 1] <- rowSums(NPDC[, , ij, 1] < 0)
  }

  TABLE <- array(NA, c(k + 4, k + 1, interval),
                 dimnames = list(c(NAMES, "TO", "Inc.Own", "Net", "NPDC"), c(NAMES, "FROM"), period_names))
  for (i in 1:interval) TABLE[, , i] <- ConnectednessTable(CT[, , , i] / 100)$TABLE

  config <- list(partition = partition, nfore = nfore, generalized = generalized, orth = orth,
                 scenario = scenario, corrected = corrected, approach = "Frequency")
  list(TABLE = TABLE, CT = CT / 100, TCI = TCI, TO = TO, FROM = FROM, NET = NET, NPT = NPT,
       NPDC = NPDC, PCI = PCI, INFLUENCE = INFLUENCE, config = config)
}

ConnectednessApproach <- function(x,
                                  nlag = 1,
                                  nfore = 10,
                                  window.size = NULL,
                                  corrected = FALSE,
                                  model = "TVP-VAR",
                                  connectedness = c("Time", "Frequency"),
                                  VAR_config = list(TVPVAR = list(kappa1 = 0.99, kappa2 = 0.99,
                                                                  prior = "BayesPrior", gamma = 0.01)),
                                  Connectedness_config = list(
                                    TimeConnectedness = list(generalized = TRUE),
                                    FrequencyConnectedness = list(partition = c(pi, pi / 2, 0),
                                                                  generalized = TRUE, scenario = "ABS")),
                                  verbose = TRUE) {
  if (!inherits(x, "zoo")) stop("Data must be a 'zoo' object")
  if (model != "TVP-VAR") stop("Only model 'TVP-VAR' is supported")
  connectedness <- if (length(connectedness) > 1) "Time" else match.arg(connectedness, c("Time", "Frequency"))
  if (nlag <= 0) stop("nlag must be a positive integer")

  k <- ncol(x)
  NAMES <- colnames(x)
  if (is.null(NAMES)) NAMES <- 1:k
  if (is.null(window.size)) window.size <- nrow(x) else window.size <- window.size - nlag

  tv <- VAR_config$TVPVAR
  prior <- switch(tv$prior,
                  MinnesotaPrior = MinnesotaPrior(gamma = tv$gamma, k = k, nlag = nlag),
                  UninformativePrior = UninformativePrior(k = k, nlag = nlag),
                  BayesPrior = BayesPrior(x = x, size = window.size, nlag = nlag),
                  stop("Unknown prior"))

  if (verbose) message("Estimating TVP-VAR")
  fit <- TVPVAR(x, configuration = list(l = c(tv$kappa1, tv$kappa2), nlag = nlag, prior = prior),
                verbose = verbose)
  B_t <- fit$B_t
  Q_t <- fit$Q_t
  DATE <- as.character(zoo::index(x))
  dimnames(Q_t)[[1]] <- dimnames(Q_t)[[2]] <- NAMES
  dimnames(Q_t)[[3]] <- DATE[(length(DATE) - dim(Q_t)[3] + 1):length(DATE)]

  if (verbose) message("Computing connectedness measures")
  if (connectedness == "Time") {
    TimeConnectedness(Phi = B_t, Sigma = Q_t, nfore = nfore,
                      generalized = Connectedness_config$TimeConnectedness$generalized,
                      corrected = corrected, verbose = verbose)
  } else {
    fc <- Connectedness_config$FrequencyConnectedness
    FrequencyConnectedness(Phi = B_t, Sigma = Q_t, nfore = nfore, partition = fc$partition,
                           generalized = fc$generalized, scenario = fc$scenario,
                           corrected = corrected, verbose = verbose)
  }
}


build_W <- function(z, zlag) {
  z <- as.matrix(z)
  t <- nrow(z)
  K <- ncol(z)
  nz <- colnames(z)
  if (is.null(nz)) nz <- paste0("z", 1:K)
  W <- matrix(NA_real_, t, K * (zlag + 1))
  cn <- character(0)
  for (j in 0:zlag) {
    cols <- (j * K + 1):((j + 1) * K)
    if (j == 0) {
      W[, cols] <- z
    } else {
      W[(j + 1):t, cols] <- z[1:(t - j), , drop = FALSE]
    }
    cn <- c(cn, if (j == 0) nz else paste0(nz, "_L", j))
  }
  colnames(W) <- cn
  W
}

BayesPriorX <- function(x, W = NULL, size = NULL, nlag, sin_z = integer(0), lambda5 = 0.001) {
  if (!inherits(x, "zoo")) stop("Data must be a 'zoo' object")
  if (nlag <= 0) stop("nlag must be a positive integer")
  if (is.null(size)) size <- nrow(x)
  if ((size - nlag) <= 0) stop("size must be larger than nlag")

  p <- ncol(x)
  q <- if (is.null(W)) 0 else ncol(W)
  if (q > 0 && nrow(W) != nrow(x)) stop("W and x must have the same number of rows")
  if (q > 0 && anyNA(W[(nlag + 1):size, ])) stop("W has NA within the prior sample")

  size <- size - nlag
  yt <- t(x[(nlag + 1):(size + nlag), ])
  m <- p + nlag * (p^2) + p * q

  Zt <- NULL
  for (i in (nlag + 1):(size + nlag)) {
    ztemp <- diag(p)
    for (j in 1:nlag) {
      xlag <- x[(i - j), 1:p]
      xtemp <- matrix(0, p, p * p)
      for (jj in 1:p) xtemp[jj, ((jj - 1) * p + 1):(jj * p)] <- xlag
      ztemp <- cbind(ztemp, xtemp)
    }
    if (q > 0) {
      wtemp <- matrix(0, p, p * q)
      for (jj in 1:p) wtemp[jj, ((jj - 1) * q + 1):(jj * q)] <- W[i, ]
      ztemp <- cbind(ztemp, wtemp)
    }
    Zt <- rbind(Zt, ztemp)
  }

  vbar <- matrix(0, m, m)
  xhy <- matrix(0, m, 1)
  for (i in 1:size) {
    zhat1 <- Zt[((i - 1) * p + 1):(i * p), ]
    vbar <- vbar + t(zhat1) %*% zhat1
    xhy <- xhy + t(zhat1) %*% as.matrix(yt[, i])
  }
  vbar <- ginv(vbar)
  aols <- vbar %*% xhy

  idx_r <- integer(0)
  if (q > 0 && length(sin_z) > 0) {
    idx_r <- unlist(lapply(sin_z, function(i) p + nlag * p^2 + (i - 1) * q + 1:q))
    libres <- setdiff(seq_len(ncol(Zt)), idx_r)
    vr <- matrix(0, length(libres), length(libres)); xr <- matrix(0, length(libres), 1)
    for (i in 1:size) {
      zr <- Zt[((i - 1) * p + 1):(i * p), libres, drop = FALSE]
      vr <- vr + t(zr) %*% zr
      xr <- xr + t(zr) %*% as.matrix(yt[, i])
    }
    aols <- matrix(0, ncol(Zt), 1)
    aols[libres, ] <- ginv(vr) %*% xr
    d <- rep(1, ncol(Zt)); d[idx_r] <- lambda5
    vbar <- vbar * outer(d, d)
  }

  sse2 <- matrix(0, p, p)
  for (i in 1:size) {
    zhat1 <- Zt[((i - 1) * p + 1):(i * p), ]
    sse2 <- sse2 + (yt[, i] - zhat1 %*% aols) %*% t(yt[, i] - zhat1 %*% aols)
  }
  hbar <- sse2 / size

  list(bprior = aols[-c(1:p)], Vprior = vbar[-c(1:p), -c(1:p)], Q = hbar,
       nphi = nlag * p^2, q = q, idx_restringidos = if (length(idx_r)) idx_r - p else integer(0))
}

ExtendPriorX <- function(prior, k, q, var_gamma, sin_z = integer(0), lambda5 = 0.001) {
  if (q == 0) return(prior)
  prior$bprior <- c(prior$bprior, rep(0, k * q))
  vg <- rep(var_gamma, k * q)
  for (i in sin_z) vg[(i - 1) * q + 1:q] <- var_gamma * lambda5^2
  prior$Vprior <- as.matrix(Matrix_blockdiag(prior$Vprior, diag(vg, k * q)))
  prior
}

Matrix_blockdiag <- function(A, B) {
  out <- matrix(0, nrow(A) + nrow(B), ncol(A) + ncol(B))
  out[1:nrow(A), 1:ncol(A)] <- A
  out[nrow(A) + 1:nrow(B), ncol(A) + 1:ncol(B)] <- B
  out
}

TVPVARX <- function(x, W = NULL,
                    configuration = list(l = c(0.99, 0.99), nlag = 1, prior = NULL),
                    verbose = TRUE, guardar_kalman = FALSE) {
  l <- as.numeric(configuration$l)
  nlag <- configuration$nlag
  prior <- configuration$prior
  if (nlag <= 0) stop("nlag must be a positive integer")
  if (!inherits(x, "zoo")) stop("Data must be a 'zoo' object")
  if (l[1] <= 0 || l[1] >= 1) stop("kappa1 must be in (0, 1)")
  if (l[2] <= 0 || l[2] >= 1) stop("kappa2 must be in (0, 1)")

  r <- ncol(x)
  t <- nrow(x)
  q <- if (is.null(W)) 0 else ncol(W)
  if (q > 0) {
    W <- as.matrix(W)
    if (nrow(W) != t) stop("W and x must have the same number of rows")
    if (anyNA(W[(nlag + 1):t, , drop = FALSE])) {
      stop("W has NA after the first nlag months")
    }
  }
  if (is.null(prior)) {
    prior <- ExtendPriorX(UninformativePrior(r, nlag), r, q, 4)
  }
  NAMES <- colnames(x)
  if (is.null(NAMES)) NAMES <- 1:r
  WNAMES <- if (q > 0) colnames(W) else NULL
  if (q > 0 && is.null(WNAMES)) WNAMES <- paste0("w", 1:q)

  mphi <- nlag * (r^2)
  m <- mphi + r * q
  if (length(prior$bprior) != m) stop("Prior has wrong dimension (", m, ")")

  beta_0.mean <- prior$bprior
  beta_0.var <- prior$Vprior
  Q <- prior$Q
  if (is.null(Q)) Q <- cov(x)

  DATES <- as.character(zoo::index(x))

  media_x <- colMeans(x)
  x <- scale(x, TRUE, FALSE)
  if (q > 0) {
    W <- sweep(W, 2, colMeans(W[(nlag + 1):t, , drop = FALSE]), "-")
  }

  k <- nlag * r
  beta_pred <- matrix(0, m, t)
  beta_update <- matrix(0, m, t)
  Rb_t <- array(0, c(m, m, t))
  Sb_t <- array(0, c(m, m, t))
  beta_t <- array(0, c(k, k, t))
  Q_t <- array(0, c(r, r, t), dimnames = list(NAMES, NAMES, DATES))
  Gamma_t <- if (q > 0) array(0, c(r, q, t), dimnames = list(NAMES, WNAMES, DATES)) else NULL

  l_2 <- l[2]
  l_4 <- l[1]

  yy <- x[(nlag + 1):t, ]
  xx <- embed(x, nlag + 1)[, -c(1:r)]

  Hrow <- function(irep) {
    xl <- xx[irep - nlag, ]
    H <- NULL
    for (j in 1:nlag) {
      H <- cbind(H, kronecker(diag(r), t(xl[((j - 1) * r + 1):(j * r)])))
    }
    if (q > 0) H <- cbind(H, kronecker(diag(r), t(W[irep, ])))
    H
  }

  B <- NULL
  G <- NULL
  if (guardar_kalman) {
    K_t <- array(0, c(m, r, t)); e_t <- matrix(0, r, t)
  }
  if (verbose) pb <- txtProgressBar(min = 0, max = t, style = 3)
  for (irep in 1:t) {

    if (irep == 1) {
      Q_t[, , irep] <- Q
    } else {
      if (irep <= (nlag + 1)) {
        Gf_t <- 0.1 * (t(matrix(x[irep, ], nrow = 1)) %*% (x[irep, ]))
      } else {
        res <- yy[(irep - nlag), ] - xx[(irep - nlag), ] %*% t(B[1:r, 1:k])
        if (q > 0) res <- res - W[irep, , drop = FALSE] %*% t(G)
        Gf_t <- t(res) %*% res
      }
      Q_t[, , irep] <- l_2 * Q_t[, , (irep - 1)] + (1 - l_2) * Gf_t[1:r, 1:r]
    }

    if (irep <= (nlag + 1)) {
      beta_pred[, irep] <- beta_0.mean
      beta_update[, irep] <- beta_pred[, irep]
      Rb_t[, , irep] <- beta_0.var
    } else {
      beta_pred[, irep] <- beta_update[, (irep - 1)]
      Rb_t[, , irep] <- (1 / l_4) * Sb_t[, , (irep - 1)]
    }

    if (irep >= (nlag + 1)) {
      H <- Hrow(irep)
      Rx <- Rb_t[, , irep] %*% t(H)
      KV_b <- Q_t[, , irep] + H %*% Rx
      KG <- Rx %*% ginv(KV_b)
      beta_update[, irep] <- matrix(beta_pred[, irep], ncol = 1) +
        KG %*% (t(matrix(x[irep, ], nrow = 1)) - H %*% matrix(beta_pred[, irep], ncol = 1))
      Sb_t[, , irep] <- Rb_t[, , irep] - KG %*% (H %*% Rb_t[, , irep])
      if (guardar_kalman) {
        K_t[, , irep] <- KG
        e_t[, irep] <- as.numeric(x[irep, ]) - as.numeric(H %*% beta_pred[, irep])
      }
    }

    bb <- matrix(beta_update[, irep], ncol = 1)
    splace <- 0
    biga <- matrix(0, r, r * nlag)
    for (ii in 1:nlag) {
      for (iii in 1:r) {
        biga[iii, ((ii - 1) * r + 1):(ii * r)] <- t(bb[(splace + 1):(splace + r), 1])
        splace <- splace + r
      }
    }
    B <- rbind(biga, cbind(diag(r * (nlag - 1)), matrix(0, nrow = r * (nlag - 1), ncol = r)))
    if (q > 0) G <- matrix(bb[(mphi + 1):m, 1], nrow = r, ncol = q, byrow = TRUE)

    if ((max(abs(eigen(B)$values)) <= 1) || (irep == 1)) {
      beta_t[, , irep] <- B
      if (q > 0) Gamma_t[, , irep] <- G
    } else {
      beta_t[, , irep] <- beta_t[, , (irep - 1)]
      if (q > 0) Gamma_t[, , irep] <- Gamma_t[, , (irep - 1)]
      beta_update[, irep] <- 0.99 * beta_update[, (irep - 1)]
    }
    if (verbose) setTxtProgressBar(pb, irep)
  }
  if (verbose) close(pb)

  B_t <- beta_t[1:ncol(Q_t), , ]
  out <- list(B_t = B_t, Q_t = Q_t, Gamma_t = Gamma_t)
  if (guardar_kalman) {
    out$kalman <- list(beta_pred = beta_pred, K = K_t, e = e_t, xc = matrix(as.numeric(x), nrow = t),
                       media = media_x, nlag = nlag, W = W, fechas = as.Date(DATES), nombres = NAMES)
  }
  out
}

ConnectednessApproachX <- function(x, z = NULL,
                                   nlag = 1,
                                   zlag = 0,
                                   nfore = 10,
                                   window.size = NULL,
                                   corrected = FALSE,
                                   connectedness = c("Time", "Frequency"),
                                   VAR_config = list(TVPVAR = list(kappa1 = 0.99, kappa2 = 0.99,
                                                                   prior = "BayesPrior", gamma = 0.01)),
                                   Connectedness_config = list(
                                     TimeConnectedness = list(generalized = TRUE),
                                     FrequencyConnectedness = list(partition = c(pi, pi / 2, 0),
                                                                   generalized = TRUE, scenario = "ABS")),
                                   x_sin_z = NULL,
                                   lambda5 = 0.001,
                                   verbose = TRUE, guardar_kalman = FALSE) {
  if (!inherits(x, "zoo")) stop("x must be a 'zoo' object")
  connectedness <- if (length(connectedness) > 1) "Time" else match.arg(connectedness, c("Time", "Frequency"))
  if (nlag <= 0) stop("nlag must be a positive integer")
  if (zlag < 0) stop("zlag must be 0 or positive")
  if (anyNA(x)) stop("x has NA: use na.omit() o complete los datos antes")

  if (!is.null(z)) {
    if (!inherits(z, "zoo")) stop("z must be a 'zoo' object")
    if (is.null(colnames(z))) colnames(z) <- paste0("z", 1:ncol(z))
    if (is.null(colnames(x))) colnames(x) <- paste0("x", 1:ncol(x))
    nx <- colnames(x)
    n_antes <- nrow(x)

    Wz <- zoo(build_W(coredata(z), zlag), order.by = zoo::index(z))
    nw <- colnames(Wz)

    pos <- match(zoo::index(x), zoo::index(Wz))
    W <- matrix(NA_real_, nrow(x), ncol(Wz), dimnames = list(NULL, nw))
    W[!is.na(pos), ] <- coredata(Wz)[pos[!is.na(pos)], , drop = FALSE]
    completo <- !apply(is.na(W), 1, any)
    if (!any(completo)) stop("No common dates between x and z (with lags).")

    f <- which(completo)[1]
    u <- max(which(completo))
    ini <- max(1, f - nlag)
    filas <- ini:u
    if (!all(completo[(ini + nlag):u])) {
      stop("z has missing months within the sample of x.")
    }
    x <- x[filas, ]
    W <- W[filas, , drop = FALSE]
    z <- zoo(coredata(z)[match(zoo::index(x), zoo::index(z)), , drop = FALSE], zoo::index(x))
    colnames(x) <- nx
    if (verbose && nrow(x) < n_antes) {
      message("Common sample with controls: ", format(start(x)), " to ", format(end(x)),
              " (", nrow(x), " of ", n_antes, " months).")
    }
  } else {
    W <- NULL
  }

  k <- ncol(x)
  q <- if (is.null(W)) 0 else ncol(W)
  NAMES <- colnames(x)
  if (is.null(NAMES)) NAMES <- 1:k
  ws <- if (is.null(window.size)) nrow(x) else window.size - nlag

  sin_z <- integer(0)
  if (!is.null(x_sin_z) && length(x_sin_z) > 0) {
    if (q == 0) stop("x_sin_z requires controls z")
    if (!all(x_sin_z %in% NAMES)) stop("x_sin_z not in columns of x: ", paste(setdiff(x_sin_z, NAMES), collapse = ", "))
    sin_z <- match(x_sin_z, NAMES)
  }

  tv <- VAR_config$TVPVAR
  prior <- switch(tv$prior,
                  BayesPrior = BayesPriorX(x = x, W = W, size = ws, nlag = nlag,
                                           sin_z = sin_z, lambda5 = lambda5),
                  MinnesotaPrior = ExtendPriorX(MinnesotaPrior(gamma = tv$gamma, k = k, nlag = nlag),
                                                k, q, tv$gamma, sin_z, lambda5),
                  UninformativePrior = ExtendPriorX(UninformativePrior(k = k, nlag = nlag), k, q, 4,
                                                    sin_z, lambda5),
                  stop("Unknown prior"))

  if (verbose) message("Estimating TVP-VARX (", q, " exogenous regressors)")
  fit <- TVPVARX(x, W = W,
                 configuration = list(l = c(tv$kappa1, tv$kappa2), nlag = nlag, prior = prior),
                 verbose = verbose, guardar_kalman = guardar_kalman)
  B_t <- fit$B_t
  Q_t <- fit$Q_t
  DATE <- as.character(zoo::index(x))
  dimnames(Q_t)[[1]] <- dimnames(Q_t)[[2]] <- NAMES
  dimnames(Q_t)[[3]] <- DATE

  if (verbose) message("Computing connectedness measures")
  if (connectedness == "Time") {
    dca <- TimeConnectedness(Phi = B_t, Sigma = Q_t, nfore = nfore,
                             generalized = Connectedness_config$TimeConnectedness$generalized,
                             corrected = corrected, verbose = verbose)
  } else {
    fc <- Connectedness_config$FrequencyConnectedness
    dca <- FrequencyConnectedness(Phi = B_t, Sigma = Q_t, nfore = nfore, partition = fc$partition,
                                  generalized = fc$generalized, scenario = fc$scenario,
                                  corrected = corrected, verbose = verbose)
  }

  dca$Gamma_t <- fit$Gamma_t
  dca$kalman <- fit$kalman
  dca$B_t <- B_t
  dca$Q_t <- Q_t
  dca$x_usado <- x
  dca$z_usado <- z
  dca$config_x <- list(nlag = nlag, zlag = zlag, controles = if (q > 0) colnames(W) else NULL,
                       prior = tv$prior, kappa = c(tv$kappa1, tv$kappa2),
                       x_sin_z = if (length(sin_z)) NAMES[sin_z] else NULL, lambda5 = lambda5)

  dca$max_abs_restringido <- if (length(sin_z) && q > 0) max(abs(fit$Gamma_t[sin_z, , ])) else NA
  dca
}

seleccion_rezagos_X <- function(x, z = NULL, p_max = 6, s_max = 3, s_igual_p = FALSE,
                                p_fijo = NULL, s_fijo = NULL, x_sin_z = NULL, verbose = TRUE) {
  if (!is.null(p_fijo) && is.na(p_fijo)) p_fijo <- NULL
  if (!is.null(s_fijo) && is.na(s_fijo)) s_fijo <- NULL
  if (!is.null(p_fijo)) p_max <- max(p_max, p_fijo)
  if (!is.null(s_fijo)) s_max <- max(s_max, s_fijo)
  if (!inherits(x, "zoo")) stop("x must be a 'zoo' object")
  if (anyNA(x)) stop("x has NA")
  N <- ncol(x)
  X <- coredata(x)
  fx <- zoo::index(x)

  if (!is.null(z)) {
    if (!inherits(z, "zoo")) stop("z must be a 'zoo' object")
    if (is.null(colnames(z))) colnames(z) <- paste0("z", 1:ncol(z))
    K <- ncol(z)
    smax_eff <- if (s_igual_p) p_max else s_max

    Wz_full <- build_W(coredata(z), smax_eff)
    pos <- match(fx, zoo::index(z))
    W_all <- matrix(NA_real_, nrow(X), ncol(Wz_full))
    W_all[!is.na(pos), ] <- Wz_full[pos[!is.na(pos)], , drop = FALSE]
  } else {
    K <- 0
    smax_eff <- 0
    W_all <- NULL
  }

  lagmax <- p_max
  ok <- rep(TRUE, nrow(X))
  ok[seq_len(lagmax)] <- FALSE
  if (!is.null(W_all)) ok <- ok & !apply(is.na(W_all), 1, any)
  filas <- which(ok)
  if (length(filas) < 2) stop("Insufficient sample")
  if (any(diff(filas) != 1)) stop("Common sample is not contiguous: check missing values in z")
  Tn <- length(filas)

  combos <- if (is.null(z)) data.frame(p = 1:p_max, s = NA) else
    if (s_igual_p) data.frame(p = 1:p_max, s = 1:p_max) else
      expand.grid(p = 1:p_max, s = 0:s_max)
  if (!is.null(p_fijo)) combos <- combos[combos$p == p_fijo, , drop = FALSE]
  if (!is.null(s_fijo) && !is.null(z)) combos <- combos[combos$s == s_fijo, , drop = FALSE]
  if (nrow(combos) == 0) stop("No combinations with these fixed values (check s_igual_p)")
  rownames(combos) <- NULL

  res <- t(sapply(seq_len(nrow(combos)), function(r) {
    p <- combos$p[r]; s <- combos$s[r]
    Y <- X[filas, , drop = FALSE]
    R <- matrix(1, Tn, 1)
    for (j in 1:p) R <- cbind(R, X[filas - j, , drop = FALSE])
    nR <- ncol(R)
    if (!is.null(z)) R <- cbind(R, W_all[filas, 1:(K * (s + 1)), drop = FALSE])

    U <- matrix(0, Tn, N); k <- 0
    for (i in 1:N) {
      cols <- if (!is.null(z) && isTRUE(colnames(x)[i] %in% x_sin_z)) 1:nR else seq_len(ncol(R))
      b <- qr.solve(R[, cols, drop = FALSE], Y[, i])
      U[, i] <- Y[, i] - R[, cols, drop = FALSE] %*% b
      k <- k + length(cols)
    }
    S <- crossprod(U) / Tn
    ld <- as.numeric(determinant(S, logarithm = TRUE)$modulus)
    c(logdet = ld, k = k,
      BIC = ld + log(Tn) / Tn * k,
      AIC = ld + 2 / Tn * k,
      HQ  = ld + 2 * log(log(Tn)) / Tn * k)
  }))
  tabla <- cbind(combos, res)
  elegido <- lapply(c(BIC = "BIC", AIC = "AIC", HQ = "HQ"), function(cr) {
    i <- which.min(tabla[[cr]]); c(nlag = tabla$p[i], zlag = tabla$s[i])
  })

  if (verbose) {
    cat("Lag selection (OLS VARX, common sample of", Tn, "months:",
        format(fx[filas[1]]), "to", format(fx[filas[Tn]]), ")\n")
    print(tabla[, c("p", "s", "BIC", "AIC", "HQ")], digits = 5, row.names = FALSE)
    for (cr in names(elegido)) {
      cat(sprintf("  %-3s: nlag = %d%s\n", cr, elegido[[cr]]["nlag"],
                  if (is.null(z)) "" else sprintf(", zlag = %d", elegido[[cr]]["zlag"])))
    }
  }
  list(tabla = tabla, elegido = elegido, muestra = fx[filas], T = Tn)
}

seleccion_rezagos_VAR <- function(x, p_max = 6, verbose = TRUE) {
  r <- seleccion_rezagos_X(x, z = NULL, p_max = p_max, verbose = verbose)
  r$tabla$s <- NULL
  r$elegido <- sapply(r$elegido, function(e) e["nlag"])
  names(r$elegido) <- c("BIC", "AIC", "HQ")
  r
}


config_bootstrap <- function(nlag0 = 1, nlag1 = 1, zlag1 = 0, nlag2 = 1, nfore = 12,
                             nfore_frec = 100, window_size = 200, kappa1 = 0.99,
                             kappa2 = 0.96, lambda5 = 0.001,
                             particion = c(pi + 0.00001, pi / 3, 0),
                             modelos = c("M0", "M1")) {
  list(nlag0 = nlag0, nlag1 = nlag1, zlag1 = zlag1, nlag2 = nlag2, nfore = nfore,
       nfore_frec = nfore_frec, window_size = window_size, kappa1 = kappa1,
       kappa2 = kappa2, lambda5 = lambda5, particion = particion, modelos = modelos)
}
modelos_cfg <- function(cfg) if (is.null(cfg$modelos)) c("M0", "M1") else cfg$modelos

resumir_t <- function(d) {
  list(TCI = as.numeric(d$TCI[, 1]), TO = unclass(d$TO), FROM = unclass(d$FROM),
       NET = unclass(d$NET), CT = apply(d$CT, 1:2, mean) * 100)
}
resumir_f <- function(d) {
  list(TCI = unclass(d$TCI), TO = unclass(d$TO), FROM = unclass(d$FROM), NET = unclass(d$NET),
       CT = apply(d$CT, c(1, 2, 4), mean) * 100)
}

estimar_modelos <- function(x, z, cfg, guardar_kalman = FALSE) {
  mods <- modelos_cfg(cfg)
  out <- list(); kal <- list()
  for (m in mods) {
    r <- estimar_un_modelo(m, x, z, cfg, guardar_kalman = guardar_kalman)
    kal[[m]] <- r$kalman; r$kalman <- NULL
    out[[m]] <- r
  }
  out$fechas <- out[[mods[1]]]$fechas; out$paises <- colnames(x); out$bandas <- out[[mods[1]]]$bandas
  for (m in mods) out[[m]]$fechas <- out[[m]]$bandas <- NULL
  if (guardar_kalman) out$kalman <- kal
  out
}

generar_pseudo <- function(kal, eta, trayectoria = c("fija", "estado")) {
  trayectoria <- match.arg(trayectoria)
  xc <- kal$xc; Tn <- nrow(xc); r <- ncol(xc); p <- kal$nlag
  m <- nrow(kal$beta_pred); mphi <- p * r^2
  W <- kal$W
  xs <- xc
  d <- rep(0, m)
  for (t in (p + 1):Tn) {
    bp <- kal$beta_pred[, t] + d

    A <- matrix(0, r, r * p); sp <- 0
    for (l in 1:p) for (i in 1:r) { A[i, (l - 1) * r + 1:r] <- bp[sp + 1:r]; sp <- sp + r }
    C <- rbind(A, cbind(diag(r * (p - 1)), matrix(0, r * (p - 1), r)))
    if (max(abs(eigen(C, only.values = TRUE)$values)) > 1) return(NULL)
    H <- NULL
    for (l in 1:p) H <- cbind(H, kronecker(diag(r), t(xs[t - l, ])))
    if (!is.null(W)) H <- cbind(H, kronecker(diag(r), t(W[t, ])))
    es <- eta[t] * kal$e[, t]
    xs[t, ] <- H %*% bp + es
    if (trayectoria == "estado") d <- d + kal$K[, , t] %*% (es - kal$e[, t])
  }
  sweep(xs, 2, kal$media, "+")
}

replica_bootstrap <- function(b, est, x, z, cfg, semilla = 12345, max_intentos = 50) {
  mods <- modelos_cfg(cfg)
  fechas <- est$kalman[[mods[1]]]$fechas
  for (k in 0:max_intentos) {
    set.seed(semilla + b + 1000003 * k)
    eta <- sample(c(-1, 1), length(fechas), replace = TRUE)
    names(eta) <- as.character(fechas)
    pseudo <- list()
    for (m in mods) {
      em <- eta[as.character(est$kalman[[m]]$fechas)]
      if (anyNA(em)) stop("Model dates do not match")
      pseudo[[m]] <- generar_pseudo(est$kalman[[m]], em)
    }
    if (any(sapply(pseudo, is.null))) next
    res <- tryCatch({
      out <- list()
      for (m in mods) {
        fe <- est$kalman[[m]]$fechas
        Xm <- zoo(pseudo[[m]], fe); colnames(Xm) <- colnames(x); Zm <- z
        r <- estimar_un_modelo(m, Xm, Zm, cfg)
        r$fechas <- r$bandas <- NULL
        out[[m]] <- r
      }
      out
    }, error = function(e) NULL)
    if (is.null(res)) next
    res$intentos <- k
    return(res)
  }
  stop("Replication ", b, ": no stable sample after ", max_intentos, " attempts")
}

estimar_un_modelo <- function(modelo, x, z, cfg, guardar_kalman = FALSE) {
  vc <- list(TVPVAR = list(kappa1 = cfg$kappa1, kappa2 = cfg$kappa2, prior = "BayesPrior", gamma = 0.01))
  if (modelo == "M0") {
    p0 <- cfg$nlag0
    pr0 <- BayesPrior(x = x, size = if (is.null(cfg$window_size)) nrow(x) else cfg$window_size - p0, nlag = p0)
    f0 <- TVPVAR(x, configuration = list(l = c(cfg$kappa1, cfg$kappa2), nlag = p0, prior = pr0),
                 verbose = FALSE, guardar_kalman = guardar_kalman)
    Q0 <- f0$Q_t
    dimnames(Q0) <- list(colnames(x), colnames(x), as.character(zoo::index(x)))
    t0 <- TimeConnectedness(Phi = f0$B_t, Sigma = Q0, nfore = cfg$nfore, generalized = TRUE,
                            corrected = FALSE, verbose = FALSE)
    q0 <- FrequencyConnectedness(Phi = f0$B_t, Sigma = Q0, nfore = cfg$nfore_frec, partition = cfg$particion,
                                 generalized = TRUE, scenario = "ABS", corrected = FALSE, verbose = FALSE)
    return(list(t = resumir_t(t0), f = resumir_f(q0), kalman = f0$kalman,
                fechas = as.Date(rownames(t0$TCI)), bandas = colnames(q0$TCI)))
  }
  if (modelo == "M1") {
    d1 <- ConnectednessApproachX(x, z, nlag = cfg$nlag1, zlag = cfg$zlag1, nfore = cfg$nfore,
                                 window.size = cfg$window_size, corrected = FALSE, connectedness = "Time",
                                 VAR_config = vc, verbose = FALSE, guardar_kalman = guardar_kalman)
    q1 <- FrequencyConnectedness(Phi = d1$B_t, Sigma = d1$Q_t, nfore = cfg$nfore_frec, partition = cfg$particion,
                                 generalized = TRUE, scenario = "ABS", corrected = FALSE, verbose = FALSE)
    return(list(t = resumir_t(d1), f = resumir_f(q1), kalman = d1$kalman,
                fechas = as.Date(rownames(d1$TCI)), bandas = colnames(q1$TCI)))
  }
  stop("unknown model")
}

juntar_replicas <- function(reps) {
  pila <- function(f) { a <- lapply(reps, f); simplify2array(a) }
  out <- list()
  for (m in intersect(c("M0", "M1"), names(reps[[1]]))) {
    o <- list(
      TCI = pila(function(r) r[[m]]$t$TCI),
      TO = pila(function(r) r[[m]]$t$TO), FROM = pila(function(r) r[[m]]$t$FROM),
      NET = pila(function(r) r[[m]]$t$NET), CT = pila(function(r) r[[m]]$t$CT),
      TCI_f = pila(function(r) r[[m]]$f$TCI),
      TO_f = pila(function(r) r[[m]]$f$TO), FROM_f = pila(function(r) r[[m]]$f$FROM),
      NET_f = pila(function(r) r[[m]]$f$NET), CT_f = pila(function(r) r[[m]]$f$CT))
    out[[m]] <- o
  }
  out$B <- length(reps)
  out$intentos <- sapply(reps, function(r) r$intentos)
  out
}
