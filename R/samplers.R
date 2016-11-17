
UpdatePhi <- function(IndividualData_all, M_all, FF, SS, p, d, maxd,individual_variable_index) {
  phi <- array(0,dim = c(maxd,p, FF*SS))
  data = IndividualData_all[individual_variable_index,]
  groupIndex <- SS*(M_all[1,]-1)+M_all[2,]
  for (j in 1:p) {
    phicount <- groupcount(groupIndex, data[j,], FF*SS, d[j])
    phi_j <- apply(phicount, c(1,2), function(x) rgamma(1,x+1,1))
    phi[1:d[j], j,] <- apply(phi_j, 1, function(x) x / sum(x))
  }
  dim(phi) <- c(maxd*p,FF * SS) #reshape to a 2D matrix
  return(phi)
}

UpdatePhiWeighted <- function(IndividualData_all, M_all, FF, SS, p, d, maxd,individual_variable_index, struc_weight) {
  phi <- array(0,dim = c(maxd,p, FF*SS))
  data <- lapply(IndividualData_all,function(x) x[individual_variable_index,])
  groupIndex <- lapply(M_all,function(x) SS*(x[1,]-1)+x[2,])
  for (j in 1:p) {
    phicount <- 0
    for(w_i in 1:length(struc_weight)){
      data_w_i <- data[[w_i]]
      phicount <- phicount + (groupcount(groupIndex[[w_i]], data_w_i[j,], FF*SS, d[j]) / struc_weight[w_i])
    }
    phi_j <- apply(phicount, c(1,2), function(x) rgamma(1,x+1,1))
    phi[1:d[j], j,] <- apply(phi_j, 1, function(x) x / sum(x))
  }
  dim(phi) <- c(maxd*p,FF * SS) #reshape to a 2D matrix
  return(phi)
}

UpdateOmega <- function(beta,M_all, FF, SS) {
  phicountcluster <- groupcount(M_all[1,],M_all[2,],FF,SS)

  cum <- t(apply(phicountcluster[,seq(SS,1)],1, cumsum))[,seq(SS,1)]
  v <- mapply(function(x,y) rbeta(1,x,y), 1 + phicountcluster[,1:(SS-1)], beta + cum[,2:SS])
  dim(v) <- c(FF, SS-1)
  v[v>1-1e-5] <- 1-1e-5
  v = cbind(v,1)
  omega <- v * t(apply(cbind(1, 1-v[,1:(SS-1)]), 1, cumprod))

  return(list(omega = omega, v = v))
}

UpdateOmegaWeighted <- function(beta,M_all, FF, SS, struc_weight) {
  phicountcluster <- 0
  for(w_i in 1:length(struc_weight)){
    M_all_w_i <- M_all[[w_i]]
    phicountcluster <- phicountcluster + (groupcount(M_all_w_i[1,],M_all_w_i[2,],FF,SS) / struc_weight[w_i])
  }

  cum <- t(apply(phicountcluster[,seq(SS,1)],1, cumsum))[,seq(SS,1)]
  v <- mapply(function(x,y) rbeta(1,x,y), 1 + phicountcluster[,1:(SS-1)], beta + cum[,2:SS])
  dim(v) <- c(FF, SS-1)
  v[v>1-1e-5] <- 1-1e-5
  v = cbind(v,1)
  omega <- v * t(apply(cbind(1, 1-v[,1:(SS-1)]), 1, cumprod))

  return(list(omega = omega, v = v))
}

UpdateLambda <- function(dHH,FF,G_all,HHdata_all) {
  lambda <- list()
  for (i in 1:length(dHH)) {
    lambdacount <- groupcount(G_all,HHdata_all[i,],FF, dHH[i])
    lam <- apply(lambdacount, c(1,2), function(x) rgamma(1,x+1,1))
    lam <- t(apply(lam, 1, function(x) x / sum(x)))
    lambda[[i]] = lam;
  }

  return(lambda)
}

UpdateLambdaWeighted <- function(dHH,FF,G_all,HHdata_all,struc_weight) {
  lambda <- list()
  for (i in 1:length(dHH)) {
    lambdacount <- 0
    for(w_i in 1:length(struc_weight)){
      HHdata_all_w_i <- HHdata_all[[w_i]]
      lambdacount <- lambdacount + (groupcount(G_all[[w_i]],HHdata_all_w_i[i,],FF, dHH[i])/ struc_weight[w_i])
    }
    lam <- apply(lambdacount, c(1,2), function(x) rgamma(1,x+1,1))
    lam <- t(apply(lam, 1, function(x) x / sum(x)))
    lambda[[i]] = lam;
  }

  return(lambda)
}

UpdatePi <- function(alpha,G_all,FF) {
  kcount <- groupcount1D(G_all, FF)
  s <- seq(FF,1)
  cum <- cumsum(kcount[s])[s]
  u <- mapply(function(x,y) rbeta(1,x,y), 1 + kcount[1:FF-1], alpha + cum[2:FF])
  u[u > 1-1e-5] <- 1-1e-5
  u <- c(u,1)
  u[FF] <- 1

  pi  <- u* cumprod(c(1,1-u[1:FF-1]))

  return(list(pi = pi, u = u))
}

UpdatePiWeighted <- function(alpha,G_all,FF,struc_weight) {
  kcount <- 0
  for(w_i in 1:length(struc_weight)){
    kcount <- kcount + (groupcount1D(G_all[[w_i]], FF)/ struc_weight[w_i])
  }
  s <- seq(FF,1)
  cum <- cumsum(kcount[s])[s]
  u <- mapply(function(x,y) rbeta(1,x,y), 1 + kcount[1:FF-1], alpha + cum[2:FF])
  u[u > 1-1e-5] <- 1-1e-5
  u <- c(u,1)
  u[FF] <- 1

  pi  <- u* cumprod(c(1,1-u[1:FF-1]))

  return(list(pi = pi, u = u))
}

UpdateAlpha <- function(aa,ab,u) {
  FF <- length(u)
  alpha <- rgamma(1,aa + FF - 1,scale = 1/(ab - sum(log(1-u[1:FF-1]))))
  return(alpha)
}

UpdateBeta <- function(ba,bb,v) {
  FF <- dim(v)[1]
  SS <- dim(v)[2]
  beta <- rgamma(1,ba + FF*(SS-1), scale = 1/(bb - sum(log(1-v[,1:SS-1])) ))
}



