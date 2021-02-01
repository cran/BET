uni <- function(x){
  Fx <- ecdf(x)
  Fx(x)
}

frac.2 <- function(a,d){
  a0 <- a
  b <- rep(NA,d)
  for (i in 1:d){
    temp <- 2^i*a0
    if (temp>1){
      b[i] <- 1
      a0 <- a0-1/2^i
    }else{
      b[i] <- 0
    }
  }
  b
}


frac2 <- function(a,d){
  if (d==1){
    t(t(sapply(a,function(x){frac.2(x,d)})))
  }else{
    t(sapply(a,function(x){frac.2(x,d)}))
  }
}


bex.centers <- function(depth){			# depth >=1
  cbind(rep((1:2^depth)/2^depth,2^depth),rep( (1:2^depth)/2^depth, rep(2^depth,2^depth)   ))-1/2^(depth+1)
}


plot.bid <- function(depth, be.ind1, be.ind2){
  xyc <- bex.centers(depth)
  BEx <- frac2(xyc[,1], depth)
  BEy <- frac2(xyc[,2], depth)

  RDx <- 2*BEx-1
  RDy <- 2*BEy-1

  be.ind1.num <- as.numeric(unlist(strsplit(be.ind1,":")))
  x.prod <- apply(RDx[,be.ind1.num,drop=F],1,prod) #1:row

  be.ind2.num <- as.numeric(unlist(strsplit(be.ind2,":")))
  y.prod <- apply(RDy[,be.ind2.num,drop=F],1,prod)

  col.ind <- x.prod*y.prod

  for (i.col in 1: nrow(xyc)){
    if (col.ind[i.col]<0){
      xycc <- xyc[i.col,]
      xp <- c(xycc[1]-1/2^(depth+1),xycc[1]+1/2^(depth+1),xycc[1]+1/2^(depth+1),xycc[1]-1/2^(depth+1))
      yp <- c(xycc[2]-1/2^(depth+1),xycc[2]-1/2^(depth+1),xycc[2]+1/2^(depth+1),xycc[2]+1/2^(depth+1))
      polygon(xp,yp,border=NA,col=rgb(0,0,1,1/4))
    }
  }
}

bet <- function(X, d, unif.margin = FALSE, cex=0.5, ...) UseMethod("bet")


bet.plot <- function(X, d, unif.margin = FALSE, cex=0.5, ...){
  if(ncol(X) != 2) stop("X does not have two columns.")
  bet.res <- BETCpp(X, d, unif.margin, asymptotic = T, test_uniformity = T, test_independence = F, independence_index = list())
  i1 <- as.character(bet.res$Interaction$X1)
  i2 <- as.character(bet.res$Interaction$X2)
  be.ind1 <- unlist(strsplit(i1, " "))[1]
  be.ind2 <- unlist(strsplit(i2, " "))[1]
  x <- uni(X[,1])
  y <- uni(X[,2])

  # par(mgp = c(1.8, 0.5, 0),mar=c(3,3,3,1))
  plot(c(0,1), c(0,1), xlab=expression(U[x]),ylab=expression(U[y]),type = "n")
  points(x,y,mgp = c(1.8, 0.5, 0),xlim=c(0,1),ylim=c(0,1),cex=cex,col=2, pch=16)
  plot.bid(d, be.ind1, be.ind2)
}

BET <- function(X, d, unif.margin = FALSE, asymptotic = TRUE, plot = FALSE, test.independence = FALSE, index = NULL, test.uniformity = TRUE){
  n <- nrow(X)
  p <- ncol(X)
  if (p == 1){
    for (i in 1:n){
      if (X[n][1] > 1 || X[n][1] < 0) stop("Data out of range [0, 1]")
    }
  }
  if (plot && (p == 2))
    bet.plot(X, d, unif.margin)
  if (plot && (p != 2)) warning("plot not available: X does not have two columns.")

  if (test.uniformity){
    test.independence <- FALSE
  }else if ((!test.independence)){
    stop("Choose uniformity or independence to test.")
  }else if (is.null(index)){
    stop("Need a list of index!")
  }

  BETCpp(X, d, unif.margin, asymptotic, test.uniformity, test.independence, index)
}

symm <- function(X, d, unif.margin = FALSE){
  n <- nrow(X)
  p <- ncol(X)
  if (p == 1){
    for (i in 1:n){
      if (X[n][1] > 1 || X[n][1] < 0) stop("Data out of range [0, 1]")
    }
  }
  symmCpp(X, d, unif.margin)
}

BETs <- function(X, d.max=4, unif.margin = FALSE, asymptotic = TRUE, plot = FALSE, test.independence = FALSE, index = NULL, test.uniformity = TRUE){
  n <- nrow(X)
  p <- ncol(X)
  temp <- BET(X, 1, unif.margin, asymptotic, FALSE, test.independence, index, test.uniformity) #BET
  bet.adj.pvalues <- rep(NA,d.max)
  bet.extreme.asymmetry <- rep(NA,d.max)

  max.abs.count.interaction <- abs(temp$Extreme.Asymmetry)
  bet.extreme.asymmetry[1] <- temp$Extreme.Asymmetry
  # table22 <- matrix(c(max.abs.count.interaction/2+n/4, -max.abs.count.interaction/2+n/4, -max.abs.count.interaction/2+n/4, max.abs.count.interaction/2+n/4), 2, 2)
  # FE22 <- fisher.test(table22,conf.int=FALSE)$p.value- dhyper(table22[1,1],n/2,n/2,n/2)/2
  FE.pvalue0 <- min(temp$p.value.bonf, 1)

  bet.adj.pvalues[1] <- FE.pvalue0
  bet.s.interaction <- temp$Interaction

  if (d.max==1){
    return(list(bet.s.pvalue=temp$p.value.bonf,bet.s.extreme.asymmetry=temp$Extreme.Asymmetry, bet.s.index=temp$Interaction, bet.s.zstatistic=temp$z.statistic))
  }else{
    for (id in 2:d.max){
      tempa <- BET(X, id, unif.margin, asymptotic, FALSE, test.independence, index, test.uniformity) #BET

      max.abs.count.interaction <- abs(tempa$Extreme.Asymmetry)
      bet.extreme.asymmetry[id] <- tempa$Extreme.Asymmetry

      if (p == 1){
        FE.pvalue <- min(tempa$p.value.bonf/(2^id-1) * ((2^id-1) - (2^(id-1)-1)), 1)
      }else{
        FE.pvalue <- min((tempa$p.value.bonf/(2^(p*id)-p*(2^id-1)-1)) * ((2^(p*id)-p*(2^id-1)-1) - (2^(p*(id-1))-p*(2^(id-1)-1)-1)), 1)
      }

      bet.adj.pvalues[id] <- FE.pvalue
      if (FE.pvalue < FE.pvalue0){
        bet.s.interaction <- tempa$Interaction
        FE.pvalue0 <- FE.pvalue
      }
    }
    bet.s.pvalue <- min(min(bet.adj.pvalues)*d.max,1)
    dp = which(bet.adj.pvalues==min(bet.adj.pvalues),arr.ind=TRUE)[1]
    if (plot && p == 2)
      bet.plot(X, dp, unif.margin)
    if (plot && p != 2) warning('plot not available: X does not have two columns.')
    bet.s.extreme.asymmetry <- bet.extreme.asymmetry[which(bet.adj.pvalues==min(bet.adj.pvalues))]
    bet.s.zstat <- abs(bet.s.extreme.asymmetry)/sqrt(n)
    return(list(bet.s.pvalue.bonf=bet.s.pvalue,bet.s.extreme.asymmetry=bet.s.extreme.asymmetry, bet.s.index=bet.s.interaction, bet.s.zstatistic=bet.s.zstat))
  }
}

BEAST <- function(X, d, subsample.percent = 1/2, B = 100, unif.margin = FALSE, lambda = NULL, test.independence = FALSE, index = NULL, test.uniformity = TRUE, p.value.method = "p", num.permutations = 100){
  n <- nrow(X)
  p <- ncol(X)
  if (p == 1){
    for (i in 1:n){
      if (X[n][1] > 1 || X[n][1] < 0) stop("Data out of range [0, 1]")
    }
  }

  if(is.null(lambda)){
    lambda <- sqrt(2 * log(2^(p * d)) / n) / 2
  }

  if (test.uniformity){
    test.independence <- FALSE
  }else if ((!test.independence)){
    stop("Choose uniformity or independence to test.")
  }else if (is.null(index)){
    stop("Need a list of index!")
  }

  if(!p.value.method %in% c("p", "s")){
    p.value.method <- "NA"
  }

  m <- n * subsample.percent

  BeastCpp(X, d, m, B, unif.margin, lambda, test.uniformity, test.independence, index, p.value.method, num.permutations)
}
