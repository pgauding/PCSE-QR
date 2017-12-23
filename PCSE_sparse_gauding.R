## Patrick Gauding - POLS 904 - PCSE with QR
## 2017-12-05

pcse_sparse <- function(object, groupN, groupT, pairwise = FALSE){
  library(pcse)
  library(Matrix)
  mc <- match.call()
  check <- class(object)
  if (!("lm" %in% check)) {
    stop("Formula object must be of class 'lm'.")
  }
  extractLMvar <- function(lm.result, groupvar) {
    isNamed <- try(inherits(groupvar, "character"), silent = TRUE)
    if (inherits(isNamed, "try-error")) {
      groupvar <- as.character(match.call()$groupvar)
      isNamed <- TRUE
    }
    if (isNamed & length(groupvar) == 1) {
      groupvar <- model.frame(lm.result)[[groupvar]]
    }
    groupvar
  }
  groupT <- extractLMvar(object, groupT)
  groupN <- extractLMvar(object, groupN)
  check <- !is.null(groupT)
  if (!check) {
    stop(paste(deparse(mc$groupT), "is not found."))
  }
  check <- !is.null(groupN)
  if (!check) {
    stop(paste(deparse(mc$groupN), "is not found."))
  }
  check <- length(groupN) == length(groupT)
  if (!check) {
    stop("Length of groupT and groupN be of equal length.")
  }
  check <- length(groupN) == dim(model.matrix(object))[1]
  if (!check) {
    stop("Length of groupN and groupT must equal nrows of using data.")
  }
  check <- is.na(groupN)
  if (any(check)) {
    stop("There must not be any missing values in the CS groupN!")
  }
  check <- is.na(groupT)
  if (any(check)) {
    stop("There must not be any missing values in the TS groupT!")
  }
  nCS <- length(na.omit(unique(groupN)))
  nTS <- length(na.omit(unique(groupT)))
  check <- nCS * nTS >= dim(model.matrix(object))[1]
  if (!check) {
    stop("There cannot be more than nCS*nTS rows in the using data!")
  }
  if ("factor" %in% class(groupN)) {
    groupN <- as.numeric(groupN)
  }
  if ("factor" %in% class(groupT)) {
    groupT <- as.numeric(groupT)
  }
  units <- unique(groupN)
  units <- na.omit(units)
  nCS <- length(units)
  time <- unique(groupT)
  time <- na.omit(time)
  nTS <- length(time)
  check <- 0
  for (i in 1:nCS) {
    if (sum(groupN == units[i]) == nTS) {
      check <- check + 1
    }
    else {
      check <- check
    }
  }
  flag <- ifelse(check == nCS, TRUE, FALSE)
  using <- data.frame(groupN = groupN, groupT = groupT)
  using <- na.omit(using)
  using$resid <- resid(object)
  using <- data.frame(using, model.matrix(object))
  ord <- order(using$groupN, using$groupT)
  using <- using[ord, ]
  units <- unique(using$groupN)
  nCS <- length(units)
  time <- unique(using$groupT)
  nTS <- length(time)
  avgN <- dim(using)[1]/nCS
  brows <- c()
  for (i in 1:nTS) {
    br <- which(using$groupT == time[i])
    check <- length(br) == nCS
    if (check) {
      brows <- c(brows, br)
    }
  }
  balanced <- using[brows, ]
  ord <- order(balanced$groupN, balanced$groupT)
  balanced <- balanced[ord, ]
  Bunits <- unique(balanced$groupN)
  BnCS <- length(Bunits)
  Btime <- unique(balanced$groupT)
  BnTS <- length(Btime)
  rect <- using
  rect$groupN <- as.numeric(rect$groupN)
  Runits <- unique(rect$groupN)
  Runits <- na.omit(Runits)
  missN <- matrix(NA, nCS, 2)
  for (i in 1:nCS) {
    if (sum(rect$groupN == Runits[i]) != nTS) {
      missN[i, 1] <- Runits[i]
      missN[i, 2] <- nTS - (sum(rect$groupN == Runits[i]))
    }
  }
  missN <- na.omit(missN)
  missT <- c()
  tmp <- c()
  if (dim(missN)[1] != 0 & dim(missN)[2] != 0) {
    for (i in 1:dim(missN)[1]) {
      tt <- time %in% rect$groupT[rect$groupN == missN[i,
                                                       1]]
      missT <- c(missT, time[!tt])
      tmp <- c(tmp, rep(missN[i, 1], missN[i, 2]))
    }
    missN <- tmp
    nM <- length(missN)
    R <- dim(rect)[1]
    C <- dim(rect)[2]
    if (R != nCS * nTS) {
      for (i in (R + 1):(nCS * nTS)) {
        rect[i, ] <- rep(NA, C)
      }
    }
    rect[c((R + 1):(R + nM)), 1] <- missN
    rect[c((R + 1):(R + nM)), 2] <- missT
  }
  ord <- order(rect$groupN, rect$groupT)
  rect <- rect[ord, ]
  Runits <- unique(rect$groupN)
  Runits <- na.omit(Runits)
  RnCS <- length(Runits)
  Rtime <- unique(rect$groupT)
  Rtime <- na.omit(Rtime)
  RnTS <- length(Rtime)
  if (flag) {
    e <- using$resid
    E <- Matrix(e, nCS, nTS, byrow = TRUE, sparse = TRUE)
    E <- t(E)
    Sigma.hat <- crossprod(E)/nTS
    X <- as.matrix(using[, 4:dim(using)[2]])
    omega <- kronecker(Sigma.hat, diag(1, nTS))
    middle <- t(X) %*% omega %*% X
    nobs <- length(e)
    dataX <- X
  }
  if (!flag) {
    if (!pairwise) {
      if (BnCS == 0 | BnTS == 0) {
        stop("Either the number of CS observations per panel ",
             "or the number of TS observations per panel ",
             "used to compute the vcov matrix is zero. You must use",
             " pairwise selection.")
      }
      e <- balanced$resid
      E <- matrix(e, BnCS, BnTS, byrow = TRUE)
      E <- qr(E)
      Sigma.hat <- crossprod(E)/BnTS
      if (avgN/2 > BnTS) {
        warning("Caution! The number of CS observations per panel, ",
                BnTS, ", used to compute the vcov matrix is less than half the",
                "average number of obs per panel in the original data.",
                "You should consider using pairwise selection.")
      }
      X <- as.matrix(rect[, 4:dim(rect)[2]])
      X[is.na(X)] <- 0
      omega <- kronecker(Sigma.hat, diag(1, nTS))
      middle <- t(X) %*% omega %*% X
      nobs <- length(resid(object))
    }
    if (pairwise) {
      V <- rect[, 4:dim(rect)[2]]
      valid <- apply(!is.na.data.frame(V), 1, prod)
      nobs <- sum(valid)
      e <- rect$resid
      e[is.na(e)] <- 0
      E <- matrix(e, RnCS, RnTS, byrow = TRUE)
      E <- t(E)
      V <- matrix(valid, RnCS, RnTS, byrow = TRUE)
      V <- t(V)
      numer <- crossprod(E)
      denom <- crossprod(V)
      denom[denom == 0] <- NA
      check <- is.na(denom)
      if (sum(check) != 0) {
        stop("Error! A CS-unit exists without any obs or without any obs in\n              common with another CS-unit. You must remove that unit from the\n              data passed to pcse().")
      }
      Sigma.hat <- numer/denom
      X <- as.matrix(rect[, 4:dim(rect)[2]])
      X[is.na(X)] <- 0
      omega <- kronecker(Sigma.hat, diag(1, nTS))
      middle <- t(X) %*% omega %*% X
    }
  }
  XX <- t(X) %*% X
  XXinv <- solve(XX)
  vcov <- XXinv %*% middle %*% XXinv
  pcse <- sqrt(diag(vcov))
  b <- summary(object)$coef[, 1]
  tstats <- b/pcse
  df <- nobs - ncol(X)
  pval <- 2 * pt(abs(tstats), df, lower.tail = FALSE)
  res <- list(vcov = vcov, pcse = pcse, b = b, tstats = tstats,
              df = df, pval = pval, pairwise = pairwise, nobs = nobs,
              nmiss = (nCS * nTS) - nobs, call = mc)
  class(res) <- "pcse"
  return(res)
}

## Generate X
gen.x <- function(t, x.bar, sd.x){
  e.x <- rnorm(t, mean = 0, sd = 1)
  x <- rnorm(t, mean = x.bar, sd = sd.x) + e.x
  x
}

## Generate Y
gen.y <- function(t, b1 = 5, x){
  e.y <- rnorm(t, mean = 0, sd = 1)
  y <- double(t)
  for(i in 1:t){
    y[i] <- b1*x[i] + e.y[i]
  }
  y
}

## Generate Data
data.gen <- function(t, b1, x.bar, sd.x){
  x <- gen.x(t = t, x.bar = x.bar, sd.x = sd.x)
  y <- gen.y(t = t, b1 = b1, x)
  return(df <- data.frame(cbind(x=x,y=y)))
}

## Generate Data Set
data.set <- function(cs, t, b1, x.bar, sd.x){
  
  out <- do.call(rbind, lapply(1:cs, function(x) data.gen(t = t, b1 = b1, x.bar = x.bar, sd.x = sd.x)))
  
  mat.x <- matrix(unlist(out[,"x"]), ncol=cs)
  mat.y <- matrix(unlist(out[,"y"]), ncol=cs)
  
  x.col <- matrix(mat.x, ncol=1)
  y.col <- matrix(mat.y, ncol=1)
  
  ind.mat <- sapply(c(1:cs), rep.int, t)
  ind.col <- matrix(ind.mat, ncol=1)
  
  time.mat <- matrix(rep(1:t,cs), ncol=cs)
  time.col <- matrix(time.mat, ncol=1)
  
  data.mat<- data.frame(cbind(ind.col, time.col, y.col, x.col))
  colnames(data.mat) <- c("cs", "t", "y", "x")
  data.mat
  
}

## Generate TSCS dataset
set.seed(281330)
cs <- 830
t <- 830
b1 <- 25
x.bar <- 50
sd.x <- 5

panel.one <- data.set(cs, t, b1, x.bar, sd.x)
panel.one

panel.two <- data.set(cs, t, b1, x.bar, sd.x)
panel.two

panel.three <- data.set(cs, t, b1, x.bar, sd.x)
panel.three

panel.four <- data.set(cs, t, b1, x.bar, sd.x)
panel.four

## PCSE Test
library(pcse)

p1_lm <- lm(y ~ x, data = panel.one)
p1_pcse <- pcse(p1_lm, groupN = panel.one$cs, groupT = panel.one$t)
p1_pcse_sparse <- pcse_sparse(p1_lm, groupN = panel.one$cs, groupT = panel.one$t)
summary(p1_pcse)
summary(p1_pcse_sparse)

p2_lm <- lm(y ~ x, data = panel.two)
p2_pcse <- pcse(p2_lm, groupN = panel.two$cs, groupT = panel.two$t)
p2_pcse_sparse <- pcse_sparse(p2_lm, groupN = panel.two$cs, groupT = panel.two$t)
summary(p2_pcse)
summary(p2_pcse_sparse)

p3_lm <- lm(y ~ x, data = panel.three)
p3_pcse <- pcse(p3_lm, groupN = panel.three$cs, groupT = panel.three$t)
p3_pcse_sparse <- pcse_sparse(p3_lm, groupN = panel.three$cs, groupT = panel.three$t)
summary(p3_pcse)
summary(p3_pcse_sparse)

p4_lm <- lm(y ~ x, data = panel.four)
#p3_pcse <- pcse(p3_lm, groupN = panel.three$cs, groupT = panel.three$t)
p4_pcse_sparse <- pcse_sparse(p4_lm, groupN = panel.four$cs, groupT = panel.four$t)
#summary(p3_pcse)
summary(p4_pcse_sparse)

library(rbenchmark)

rbenchmark::benchmark(
    #"pcse" = {
   # p1_pcse <- pcse(p1_lm, groupN = panel.one$cs, groupT = panel.one$t)
  #},
    "pcse_sparse" = {
    p1_pcse_sparse <- pcse_sparse(p1_lm, groupN = panel.one$cs, groupT = panel.one$t)
  },
    replications = 1,
    columns = c("test", "replications", "elapsed", "relative",
            "user.self", "sys.self"))

rbenchmark::benchmark(
    "pcse" = {
    p1_pcse <- pcse(p2_lm, groupN = panel.two$cs, groupT = panel.two$t)
  },
    "pcse_sparse" = {
    p1_pcse_sparse <- pcse_sparse(p2_lm, groupN = panel.two$cs, groupT = panel.two$t)
  },
    replications = 1000,
    columns = c("test", "replications", "elapsed", "relative",
              "user.self", "sys.self"))

rbenchmark::benchmark(
    "pcse" = {
    p3_pcse_sparse <- pcse_sparse(p3_lm, groupN = panel.three$cs, groupT = panel.three$t)
  },
    replications = 1,
    columns = c("test", "replications", "elapsed", "relative",
              "user.self", "sys.self"))

library(profvis)

profvis({
  p1_pcse <- pcse(p1_lm, groupN = panel.one$cs, groupT = panel.one$t)
  p1_pcse_sparse <- pcse_sparse(p1_lm, groupN = panel.one$cs, groupT = panel.one$t)
})

