
ca_jo_jus06_hr3_fun <- function(z,
                                H1, H2, H3,
                                r1, r2, r3,
                                nor.id1, nor.id2, nor.id3,
                                conv.val = 1e-04, max.iter = 50, df = 1) {

  # # inputs
  # ca.jo.res.02$P <- 6
  # z <- ca.jo.res.02
  #
  # H1 <- matrix(c( 1, 0, 0,
  #                -1, 0, 0,
  #                 0, 0, 0,
  #                 0, 1, 0,
  #                 0,-1, 0,
  #                 0, 0, 1), nrow=6, ncol=3, byrow=T)
  # H2 <- matrix(c( 0, 0, 0,
  #                 1, 0, 0,
  #                 0, 1, 0,
  #                 0, 0, 0,
  #                 0, 0, 0,
  #                 0, 0, 1), nrow=6, ncol=3, byrow=T)
  # H3 <- matrix(c( 0, 0, 0,
  #                 0, 0, 0,
  #                 0, 0, 0,
  #                 1, 0, 0,
  #                 0, 1, 0,
  #                 0, 0, 1), nrow=6, ncol=3, byrow=T)
  #
  # r1 <- 1
  # r2 <- 1
  # r3 <- 1
  #
  # nor.id1 <- 1
  # nor.id2 <- 3
  # nor.id3 <- 4
  #
  # conv.val <- 1e-04
  # max.iter <- 50
  #
  # df <- 1





  P <- z$P

  # # if (!(class(z) == "ca.jo")) {
  # #   stop("\nPlease, provide object of class 'ca.jo' as 'z'.\n")
  # # }
  # if (r >= z$P || r < 1) {
  #   stop("\nCount of cointegrating relationships is out of allowable range.\n")
  # }
  # if (z$ecdet == "none") {
  #   P <- z$P
  # } else {
  #   P <- z$P + 1
  # }
  # r <- as.integer(r)
  H1 <- as.matrix(H1)
  H2 <- as.matrix(H2)
  H3 <- as.matrix(H3)
  # if (!(nrow(H1) == P)) {
  #   stop("\nRow number of 'H' is unequal to VAR order.\n")
  # }
  s1 <- ncol(H1)
  s2 <- ncol(H2)
  s3 <- ncol(H3)
  r1 <- 1
  r2 <- 1
  r3 <- 1
  # r2 <- r - r1
  r <- r1 + r2 + r3

  lambda <- z$lambda
  type <- "Estimation and testing under partly known beta"
  N <- nrow(z$Z0)
  M00 <- crossprod(z$Z0)/N
  M11 <- crossprod(z$Z1)/N
  MKK <- crossprod(z$ZK)/N
  M01 <- crossprod(z$Z0, z$Z1)/N
  M0K <- crossprod(z$Z0, z$ZK)/N
  MK0 <- crossprod(z$ZK, z$Z0)/N
  M10 <- crossprod(z$Z1, z$Z0)/N
  M1K <- crossprod(z$Z1, z$ZK)/N
  MK1 <- crossprod(z$ZK, z$Z1)/N
  M11inv <- solve(M11)
  S00 <- M00 - M01 %*% M11inv %*% M10
  S0K <- M0K - M01 %*% M11inv %*% M1K
  SK0 <- MK0 - MK1 %*% M11inv %*% M10
  SKK <- MKK - MK1 %*% M11inv %*% M1K

  bet.i1 <- as.matrix(z$Vorg[,1])
  bet.i2 <- as.matrix(z$Vorg[,2])
  bet.i3 <- as.matrix(z$Vorg[,3])

  i <- 0
  last <- 1
  diff <- 1
  while (diff > conv.val) {
    # while (i < max.iter) {

    # 1) Concentrate out given bet_1 and bet_2 and solve for bet_3 (see Johansen, p. 109, th. 7.3)
    beta1 <- cbind(bet.i1, bet.i2)
    S00.b1 <- S00 - S0K %*% beta1 %*% solve(t(beta1) %*% SKK %*% beta1) %*% t(beta1) %*% SK0
    S0K.b1 <- S0K - S0K %*% beta1 %*% solve(t(beta1) %*% SKK %*% beta1) %*% t(beta1) %*% SKK
    SK0.b1 <- SK0 - SKK %*% beta1 %*% solve(t(beta1) %*% SKK %*% beta1) %*% t(beta1) %*% SK0
    SKK.b1 <- SKK - SKK %*% beta1 %*% solve(t(beta1) %*% SKK %*% beta1) %*% t(beta1) %*% SKK

    Ctemp <- chol(t(H3) %*% SKK.b1 %*% H3, pivot = TRUE)
    pivot <- attr(Ctemp, "pivot")
    oo <- order(pivot)
    C <- t(Ctemp[, oo])
    Cinv <- solve(C)
    S00inv.b1 <- solve(S00.b1)
    valeigen <- eigen(Cinv %*% t(H3) %*% SK0.b1 %*% S00inv.b1 %*% S0K.b1 %*% H3 %*% t(Cinv))
    lam  <- valeigen$values # MM: added
    e <- as.matrix(valeigen$vectors[, 1:r3])
    bet.i3 <- H3 %*% t(Cinv) %*% e

    # L max for given bet_1 and bet_2 (see Johansen, p. 109, th. 7.3)
    Dtemp <- chol(S00, pivot = TRUE)
    pivot <- attr(Dtemp, "pivot")
    oo <- order(pivot)
    D <- t(Dtemp[, oo])
    Dinv <- solve(D)
    valeigen <- eigen(Dinv %*% S0K %*% beta1 %*% solve(t(beta1) %*% SKK %*% beta1) %*% t(beta1) %*% SK0 %*% t(Dinv))
    rho <- valeigen$values
    teststat <- N * (log(1 - rho[1]) + log(1 - rho[2]) + log(1 - lam[1]) - sum(log(1 - lambda[1:r])))
    lik.max.01 <- det(S00) * (1 - rho[1]) * (1 - rho[2]) * (1 - lam[1])

    # for convergence (Pfaff)
    lambda.res <- valeigen$values
    lambda.res <- Re(valeigen$values)
    diff <- t(lambda.res - last) %*% (lambda.res - last)
    last <- lambda.res

    # 2) Concentrate out given bet_2 and bet_3 and solve for bet_1
    beta2 <- cbind(bet.i2, bet.i3)
    S00.b2 <- S00 - S0K %*% beta2 %*% solve(t(beta2) %*% SKK %*% beta2) %*% t(beta2) %*% SK0
    S0K.b2 <- S0K - S0K %*% beta2 %*% solve(t(beta2) %*% SKK %*% beta2) %*% t(beta2) %*% SKK
    SK0.b2 <- SK0 - SKK %*% beta2 %*% solve(t(beta2) %*% SKK %*% beta2) %*% t(beta2) %*% SK0
    SKK.b2 <- SKK - SKK %*% beta2 %*% solve(t(beta2) %*% SKK %*% beta2) %*% t(beta2) %*% SKK

    Ctemp <- chol(t(H1) %*% SKK.b2 %*% H1, pivot = TRUE)
    pivot <- attr(Ctemp, "pivot")
    oo <- order(pivot)
    C <- t(Ctemp[, oo])
    Cinv <- solve(C)
    S00inv.b2 <- solve(S00.b2)
    valeigen <- eigen(Cinv %*% t(H1) %*% SK0.b2 %*% S00inv.b2 %*% S0K.b2 %*% H1 %*% t(Cinv))
    lam  <- valeigen$values # MM: added
    e <- as.matrix(valeigen$vectors[, 1:r1])
    bet.i1 <- H1 %*% t(Cinv) %*% e

    # L max for given bet_2 and bet_3
    # Dtemp <- chol(S00, pivot = TRUE)
    # pivot <- attr(Dtemp, "pivot")
    # oo <- order(pivot)
    # D <- t(Dtemp[, oo])
    # Dinv <- solve(D)
    valeigen <- eigen(Dinv %*% S0K %*% beta2 %*% solve(t(beta2) %*% SKK %*% beta2) %*% t(beta2) %*% SK0 %*% t(Dinv))
    rho <- valeigen$values
    teststat <- N * (log(1 - rho[1]) + log(1 - rho[2]) + log(1 - lam[1]) - sum(log(1 - lambda[1:r])))
    lik.max.02 <- det(S00) * (1 - rho[1]) * (1 - rho[2]) * (1 - lam[1])

    # for convergence (CATS)
    diff.cats <- abs(lik.max.02 - lik.max.01)

    # 3) Concentrate out given bet_1 and bet_3 and solve for bet_2
    beta3 <- cbind(bet.i1, bet.i3)
    S00.b3 <- S00 - S0K %*% beta3 %*% solve(t(beta3) %*% SKK %*% beta3) %*% t(beta3) %*% SK0
    S0K.b3 <- S0K - S0K %*% beta3 %*% solve(t(beta3) %*% SKK %*% beta3) %*% t(beta3) %*% SKK
    SK0.b3 <- SK0 - SKK %*% beta3 %*% solve(t(beta3) %*% SKK %*% beta3) %*% t(beta3) %*% SK0
    SKK.b3 <- SKK - SKK %*% beta3 %*% solve(t(beta3) %*% SKK %*% beta3) %*% t(beta3) %*% SKK

    Ctemp <- chol(t(H2) %*% SKK.b3 %*% H2, pivot = TRUE)
    pivot <- attr(Ctemp, "pivot")
    oo <- order(pivot)
    C <- t(Ctemp[, oo])
    Cinv <- solve(C)
    S00inv.b3 <- solve(S00.b3)
    valeigen <- eigen(Cinv %*% t(H2) %*% SK0.b3 %*% S00inv.b3 %*% S0K.b3 %*% H2 %*% t(Cinv))
    lam  <- valeigen$values # MM: added
    e <- as.matrix(valeigen$vectors[, 1:r3])
    bet.i2 <- H2 %*% t(Cinv) %*% e

    # L max for given bet_1 and bet_3
    # Dtemp <- chol(S00, pivot = TRUE)
    # pivot <- attr(Dtemp, "pivot")
    # oo <- order(pivot)
    # D <- t(Dtemp[, oo])
    # Dinv <- solve(D)
    valeigen <- eigen(Dinv %*% S0K %*% beta3 %*% solve(t(beta3) %*% SKK %*% beta3) %*% t(beta3) %*% SK0 %*% t(Dinv))
    rho <- valeigen$values
    teststat <- N * (log(1 - rho[1]) + log(1 - rho[2]) + log(1 - lam[1]) - sum(log(1 - lambda[1:r])))
    lik.max.02 <- det(S00) * (1 - rho[1]) * (1 - rho[2]) * (1 - lam[1])

    i <- i + 1
    if (i > max.iter) {
      warning("\nNo convergence, used last iterations values.\n")
      break
    }
  }

  bet.i1.nor <- bet.i1/bet.i1[nor.id1,1]; bet.i1.nor
  bet.i2.nor <- bet.i2/bet.i2[nor.id2,1]; bet.i2.nor
  bet.i3.nor <- bet.i3/bet.i3[nor.id3,1]; bet.i3.nor
  bet.nor <- cbind(bet.i1.nor, bet.i2.nor, bet.i3.nor)

  V <- cbind(bet.i1.nor, bet.i2.nor, bet.i3.nor)
  W <- S0K %*% V %*% solve(t(V) %*% SKK %*% V)

  PI <- W %*% t(V)
  DELTA <- S00 - S0K %*% V %*% solve(t(V) %*% SKK %*% V) %*% t(V) %*% SK0
  GAMMA <- M01 %*% M11inv - PI %*% MK1 %*% M11inv





  # df <- (P - s1 - r2) * r1
  pval <- c(1 - pchisq(Re(teststat), df), df)

  # ret.lis <- list("cajo.test", Z0 = z$Z0, Z1 = z$Z1, ZK = z$ZK, ecdet = z$ecdet,
  #                 H = H, A = NULL, B = NULL, type = type, teststat = teststat,
  #                 pval = pval, lambda = lambda.res, Vorg = Vorg, V = V,
  #                 W = W, PI = PI, DELTA = DELTA, DELTA.bb = NULL, DELTA.ab = NULL,
  #                 DELTA.aa.b = NULL, GAMMA = GAMMA, test.name = "Johansen-Procedure")


  psi.01 <- bet.i1
  psi.02 <- bet.i2
  psi.03 <- bet.i3
  Psi <- cbind(psi.01, psi.02, psi.02)


  ret.lis <- list(lambda = lambda.res,
                  V = V, W = W, PI = PI, DELTA = DELTA,
                  DELTA.bb = NULL, DELTA.ab = NULL, DELTA.aa.b = NULL,
                  GAMMA = GAMMA,
                  test.name = "Johansen-Procedure",
                  z = z,
                  H1 = H1, H2 = H2, H3 = H3,
                  r1 = r1, r2 = r2, H3 = H3,
                  nor.id1 = nor.id1, nor.id2 = nor.id2, nor.id3 = nor.id3,
                  lag = z$lag, ecdet = z$ecdet, x = z$x,
                  bet.nor = bet.nor,
                  teststat = teststat, pval = pval)

}
