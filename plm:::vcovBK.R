## Patrick Gauding <patrick.gauding@ku.edu>
## Last edited: 2019-10-07

## plm:::vcovBK.R

## > plm:::vcovBK.plm
function (x, type = c("HC0", "HC1", "HC2", "HC3", "HC4"), cluster = c("group", 
    "time"), diagonal = FALSE, ...) 
{
    if (!is.null(x$weights)) 
        stop("vcovXX functions not implemented for weighted panel regressions")
    type <- match.arg(type)
    model <- describe(x, "model")
    if (!model %in% c("random", "within", "pooling", "fd")) {
        stop("Model has to be either random, within, pooling or fd model")
    }
    demy <- pmodel.response(x, model = model)
    demX <- model.matrix(x, model = model, rhs = 1, cstcovar.rm = "all")
    if (!is.null(x$aliased) && any(x$aliased, na.rm = TRUE)) 
        demX <- demX[, !x$aliased, drop = FALSE]
    if (length(formula(x))[2] > 1) {
        demZ <- model.matrix(x, model = model, rhs = 2, cstcovar.rm = "all")
        demX <- fitted(lm.fit(demZ, demX))
    }
    pdim <- pdim(x)
    nT <- pdim$nT$N
    Ti <- pdim$Tint$Ti
    k <- dim(demX)[[2]]
    n0 <- pdim$nT$n
    t0 <- pdim$nT$T
    uhat <- x$residuals
    groupind <- as.numeric(attr(x$model, "index")[, 1])
    timeind <- as.numeric(attr(x$model, "index")[, 2])
    if (model == "fd") {
        groupind <- groupind[timeind > 1]
        timeind <- timeind[timeind > 1]
        nT <- nT - n0
        Ti <- Ti - 1
        t0 <- t0 - 1
    }
    switch(match.arg(cluster), group = {
        n <- n0
        t <- t0
        relevant.ind <- groupind
        lab <- timeind
    }, time = {
        n <- t0
        t <- n0
        relevant.ind <- timeind
        lab <- groupind
    })
    tind <- vector("list", n)
    tlab <- vector("list", n)
    for (i in 1:length(unique(relevant.ind))) {
        tind[[i]] <- which(relevant.ind == i)
        tlab[[i]] <- lab[which(relevant.ind == i)]
    }
    dhat <- function(x) {
        tx <- t(x)
        diag(crossprod(tx, solve(crossprod(x), tx)))
    }
    switch(match.arg(type), HC0 = {
        diaghat <- NULL
    }, HC1 = {
        diaghat <- NULL
    }, HC2 = {
        diaghat <- try(dhat(demX), silent = TRUE)
    }, HC3 = {
        diaghat <- try(dhat(demX), silent = TRUE)
    }, HC4 = {
        diaghat <- try(dhat(demX), silent = TRUE)
    })
    df <- nT - k
    switch(match.arg(type), HC0 = {
        omega <- function(residuals, diaghat, df) residuals
    }, HC1 = {
        omega <- function(residuals, diaghat, df) residuals * 
            sqrt(length(residuals)/df)
    }, HC2 = {
        omega <- function(residuals, diaghat, df) residuals/sqrt(1 - 
            diaghat)
    }, HC3 = {
        omega <- function(residuals, diaghat, df) residuals/(1 - 
            diaghat)
    }, HC4 = {
        omega <- function(residuals, diaghat, df) residuals/sqrt(1 - 
            diaghat)^pmin(4, length(residuals) * diaghat/as.integer(round(sum(diaghat), 
            digits = 0)))
    })
    uhat <- omega(uhat, diaghat, df)
    tres <- array(dim = c(t, t, n))
    for (i in 1:n) {
        ut <- uhat[tind[[i]]]
        tpos <- (1:t)[unique(lab) %in% tlab[[i]]]
        if (diagonal) {
            tres[tpos, tpos, i] <- diag(diag(ut %o% ut))
        }
        else {
            tres[tpos, tpos, i] <- ut %o% ut
        }
    }
    OmegaT <- apply(tres, 1:2, mean, na.rm = TRUE)
    unlabs <- unique(lab)
    salame <- array(dim = c(k, k, n))
    for (i in 1:n) {
        groupinds <- tind[[i]]
        grouplabs <- tlab[[i]]
        xi <- demX[groupinds, , drop = FALSE]
        tpos <- unlabs %in% grouplabs
        OmegaTi <- OmegaT[tpos, tpos, drop = FALSE]
        salame[, , i] <- crossprod(xi, OmegaTi) %*% xi
    }
    salame <- apply(salame, 1:2, sum)
    pane <- solve(crossprod(demX))
    mycov <- pane %*% salame %*% pane
    attr(mycov, which = "cluster") <- match.arg(cluster)
    return(mycov)
}
<bytecode: 0x7fb80b008238>
<environment: namespace:plm>
