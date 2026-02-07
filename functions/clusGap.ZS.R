clusGap.ZS <- function (x, FUNcluster, K.max, B = 100, d.power = 1, spaceH0 = c("scaledPCA", "original"), verbose = interactive(), cores=getOption("mc.cores", 2L),...) 
{
    stopifnot(is.function(FUNcluster), 
              length(dim(x)) == 2, 
              K.max >= 1, #JCV changed from 2 to 1
              (n <- nrow(x)) >= 1, ncol(x) >= 1)

    if (B != (B. <- as.integer(B)) || (B <- B.) <= 0)stop("'B' has to be a positive integer")

    cl. <- match.call()
    if (is.data.frame(x)) x <- as.matrix(x)
    ii <- seq_len(n)
    W.k <- function(X, kk) {
        if (kk > 1){ 
            clus <- FUNcluster(X, kk)
        }else{
            clus <- rep.int(1L, nrow(X))
        }
        0.5 * sum(vapply(split(ii, clus), function(I) {xs <- X[I, , drop = FALSE]; sum(dist(xs)^d.power/nrow(xs))}, 0))
    }
    logW <- E.logW <- SE.sim <- numeric(K.max)
    if (verbose)cat("Clustering k = 1,2,..., K.max (= ", K.max, "): .. ", sep = "")
    # for (k in 1:K.max) logW[k] <- log(W.k(x, k))
    par.fun <- function(j) log(W.k(X=x, kk=j))
    logW <- unlist(parallel::mclapply(as.list(seq(1,K.max)), FUN=par.fun))
    if (verbose) cat("done\n")
    spaceH0 <- match.arg(spaceH0)
    xs <- scale(x, center = TRUE, scale = FALSE)
    m.x <- rep(attr(xs, "scaled:center"), each = n)
    switch(spaceH0, 
           scaledPCA = {V.sx <- svd(xs, nu = 0)$v; xs <- xs %*% V.sx}, 
           original = {}, 
           stop("invalid 'spaceH0':", spaceH0))
    rng.x1 <- apply(xs, 2L, range)
    # logWks <- matrix(0, B, K.max)
    # if(verbose) cat("Bootstrapping, b = 1,2,..., B (= ", B, ")  [one \".\" per sample]:\n", sep = "")
    par.fun <- function(b){
        z1 <- apply(rng.x1, 2, function(M, nn) runif(nn, min = M[1], max = M[2]), nn = n)
        z <- switch(spaceH0, scaledPCA = tcrossprod(z1, V.sx), 
            original = z1) + m.x
        # for (k in 1:K.max) {
        sapply(1:K.max,FUN=function(k) log(W.k(z, k)))
        # }
        # tmplogWks <- unlist(pbmcapply::pbmclapply(1:K.max, FUN=function(k) log(W.k(z, k))))
        # logWks[b,1:K.max] <- tmplogWks
        # if (verbose) cat(".", if (b%%50 == 0) paste(b, "\n"))
    }
    logWks <- t(simplify2array(pbmcapply::pbmclapply(1:B, FUN=function(b) par.fun(b),
                                                     mc.cores = pmin(length(seq(1,K.max)),cores))))

    # if (verbose && (B%%50 != 0)) cat("", B, "\n")
    E.logW <- colMeans(logWks)
    SE.sim <- sqrt((1 + 1/B) * apply(logWks, 2, var))
    structure(class = "clusGap", list(Tab = cbind(logW, E.logW, 
        gap = E.logW - logW, SE.sim), call = cl., spaceH0 = spaceH0, 
        n = n, B = B, FUNcluster = FUNcluster))
}
plot.clusGap <- function (gap_stat, linecolor = "steelblue", maxSE = list(method = "firstSEmax", 
    SE.factor = 1)) 
{
    if (!inherits(gap_stat, "clusGap")) 
        stop("Only an object of class clusGap is allowed. (cluster package)")
    if (is.list(maxSE)) {
        if (is.null(maxSE$method)) 
            maxSE$method = "firstmax"
        if (is.null(maxSE$SE.factor)) 
            maxSE$SE.factor = 1
    }else stop("The argument maxSE must be a list containing the parameters method and SE.factor")
    gap <- gap_stat$Tab[, "gap"]
    se <- gap_stat$Tab[, "SE.sim"]
    decr <- diff(gap) <= 0
    k <- .maxSE(gap_stat$Tab[, "gap"], gap_stat$Tab[, "SE.sim"],
                method = maxSE$method, SE.factor = maxSE$SE.factor)
    df <- as.data.frame(gap_stat$Tab, stringsAsFactors = TRUE)
    df$clusters <- 1:nrow(df)
    df$ymin <- df$gap - df$SE.sim
    df$ymax <- df$gap + df$SE.sim

    plot(df$clusters, df$gap,
         ylim = range(c(df$ymin,df$ymax)),
         type='l', col = linecolor,
         ylab = "Gap statistic (k)",
         xlab = "Number of clusters k",
         main = paste("Optimal number of clusters is", k),
         xaxs='i', xlim = c(0.5,max(df$clusters)+0.5))
    points(df$clusters, df$gap, col = linecolor, pch=16)
    segments(df$clusters,
             df$ymin,
             df$clusters,
             df$ymax,
             col = linecolor)
    abline(v = k, col = linecolor, lty = 3)
    k
}
.maxSE <- function (f, SE.f, method = c("firstSEmax", "Tibs2001SEmax", 
                                        "globalSEmax", "firstmax", "globalmax"), SE.factor = 1) 
{
  method <- match.arg(method)
  stopifnot((K <- length(f)) >= 1, K == length(SE.f), SE.f >= 
              0, SE.factor >= 0)
  fSE <- SE.factor * SE.f
  switch(method, firstmax = {
    decr <- diff(f) <= 0
    if (any(decr)) which.max(decr) else K
  }, globalmax = {
    which.max(f)
  }, Tibs2001SEmax = {
    g.s <- f - fSE
    if (any(mp <- f[-K] >= g.s[-1])) which.max(mp) else K
  }, firstSEmax = {
    decr <- diff(f) <= 0
    nc <- if (any(decr)) which.max(decr) else K
    if (any(mp <- f[seq_len(nc - 1)] >= f[nc] - fSE[nc])) which(mp)[1] else nc
  }, globalSEmax = {
    nc <- which.max(f)
    if (any(mp <- f[seq_len(nc - 1)] >= f[nc] - fSE[nc])) which(mp)[1] else nc
  })
}
