
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
  return(beta)
}

SampleMissing <- function(MissData,para,orig, household_variable_index,individual_variable_index,
                          G_household,M,hyper){
  nonstruc_zero_variables <- MissData$nonstruc_zero_variables
  nonstruc_zero_variables_house <-
    nonstruc_zero_variables[is.element(nonstruc_zero_variables,household_variable_index)]
  nonstruc_zero_variables_indiv <-
    nonstruc_zero_variables[is.element(nonstruc_zero_variables,individual_variable_index)]
  struc_zero_variables <- MissData$struc_zero_variables
  struc_zero_variables_house <-
    struc_zero_variables[is.element(struc_zero_variables,household_variable_index)]
  struc_zero_variables_indiv <-
    struc_zero_variables[is.element(struc_zero_variables,individual_variable_index)]

  #sample non structural zeros variables for everyone at once
  if(sum(is.na(MissData$household_with_miss[,nonstruc_zero_variables_indiv])) > 0){
    phi_m_g <- t(para$phi[,(M + (G_household$G_Individuals-1)*hyper$SS)])
    for(k in nonstruc_zero_variables_indiv){
      if(sum(is.na(MissData$household_with_miss[,k]))>0){
        real_k <- which(individual_variable_index==k)
        pr_X_miss_p <- phi_m_g[which(is.na(MissData$household_with_miss[,k])==TRUE),
                               ((1:orig$d[real_k]) + (real_k-1)*orig$maxd)]
        Ran_unif_miss_p <- runif(nrow(pr_X_miss_p))
        cumul_miss_p <- pr_X_miss_p%*%upper.tri(diag(ncol(pr_X_miss_p)),diag=TRUE)
        MissData$household[is.na(MissData$household_with_miss[,k]),k] <- rowSums(Ran_unif_miss_p>cumul_miss_p) + 1L
      }
    }
  }
  if(sum(is.na(MissData$household_with_miss[,nonstruc_zero_variables_house])) > 0){
    lambda_g <- lapply(para$lambda,function(x) x[G_household$G,])
    for(kk in nonstruc_zero_variables_house){
      if(sum(is.na(MissData$household_with_miss[,kk]))>0){
        real_kk <- which(household_variable_index==kk)
        lambda_g_kk <- lambda_g[[real_kk]]
        pr_X_miss_p <-
          lambda_g_kk[which(is.na(MissData$household_with_miss[c(1,cumsum(orig$n_i[-orig$n])+1),kk])==TRUE),]
        Ran_unif_miss_p <- runif(nrow(pr_X_miss_p))
        cumul_miss_p <- pr_X_miss_p%*%upper.tri(diag(ncol(pr_X_miss_p)),diag=TRUE)
        sampled_values <- rowSums(Ran_unif_miss_p>cumul_miss_p) + 1L
        MissData$household[is.na(MissData$household_with_miss[,kk]),kk] <- rep(sampled_values,orig$n_i[
          which(is.na(MissData$household_with_miss[c(1,cumsum(orig$n_i[-orig$n])+1),kk])==TRUE)])
      }
    }
  }

  #sample structural zeros variables one household at a time
  for(sss in MissData$miss_Hhindex){
    another_index <- which(MissData$household_with_miss$Hhindex==sss)
    X_house_sss_prop <- MissData$household[another_index[1],household_variable_index]
    X_house_sss_prop <- matrix(rep(t(X_house_sss_prop),MissData$n_batch_imp[sss]),byrow=T,
                               ncol=length(household_variable_index))
    X_indiv_sss_prop <- MissData$household[another_index,individual_variable_index]
    relate_index <- which(colnames(X_indiv_sss_prop)=="relate")
    X_indiv_sss_prop <- matrix(rep(t(X_indiv_sss_prop),MissData$n_batch_imp[sss]),byrow=T,
                               ncol=length(individual_variable_index))
    NA_error_house_sss <- MissData$household_with_miss[another_index[1],household_variable_index]
    NA_error_house_sss <-  matrix(rep(t(NA_error_house_sss ),MissData$n_batch_imp[sss]),byrow=T,
                                  ncol=length(household_variable_index))
    NA_error_indiv_sss <- MissData$household_with_miss[another_index,individual_variable_index]
    NA_error_indiv_sss <- matrix(rep(t(NA_error_indiv_sss),MissData$n_batch_imp[sss]),byrow=T,
                                 ncol=length(individual_variable_index))
    G_prop <- rep(G_household$G[sss],MissData$n_batch_imp[sss])
    M_prop <- rep(M[another_index],MissData$n_batch_imp[sss])
    lambda_g <- lapply(para$lambda,function(x) x[G_prop,])
    phi_m_g <- t(para$phi[,(M_prop + (G_household$G[sss]-1)*hyper$SS)])
    check_counter_sss <- 0;
    while(check_counter_sss < 1){
      for(kkk in struc_zero_variables_house){
        real_kkk <- which(household_variable_index==kkk)
        if(sum(is.na(NA_error_house_sss[,real_kkk]))>0){
          pr_X_house_k <- lambda_g[[real_kkk]]
          Ran_unif_X_house_k <- runif(nrow(pr_X_house_k))
          cumul_X_house_k <- pr_X_house_k%*%upper.tri(diag(ncol(pr_X_house_k)),diag=TRUE)
          X_house_sss_prop[,real_kkk] <- rowSums(Ran_unif_X_house_k>cumul_X_house_k) + 1L
        }
      }
      for(kkkk in struc_zero_variables_indiv){
        real_kkkk <- which(individual_variable_index==kkkk)
        if(sum(is.na(NA_error_indiv_sss[,real_kkkk]))>0){
          pr_X_indiv_k <- phi_m_g[is.na(NA_error_indiv_sss[,real_kkkk]),
                                  ((1:orig$d[real_kkkk]) + (real_kkkk-1)*orig$maxd)]
          Ran_unif_X_indiv_k <- runif(nrow(pr_X_indiv_k))
          cumul_X_indiv_k <- pr_X_indiv_k%*%upper.tri(diag(ncol(pr_X_indiv_k)),diag=TRUE)
          X_indiv_sss_prop[is.na(NA_error_indiv_sss[,real_kkkk]),real_kkkk] <-
            rowSums(Ran_unif_X_indiv_k>cumul_X_indiv_k) + 1L
        }
      }
      #Check edit rules; Need to make this part more general, very specific for this data and assumes head is
      #at the household level
      X_indiv_sss_prop_orig <- X_indiv_sss_prop
      X_indiv_sss_prop_orig[,relate_index] <- X_indiv_sss_prop_orig[,relate_index] + 1 #recode relate
      comb_to_check <- X_house_sss_prop[,-1]
      comb_to_check[,relate_index] <- 1 #Set relate to 1
      comb_to_check <- cbind(comb_to_check,matrix(t(X_indiv_sss_prop_orig),nrow=MissData$n_batch_imp[sss],byrow=TRUE))
      check_counter <- checkSZ(comb_to_check,(length(another_index) + 1))
      check_counter_sss <- check_counter_sss + sum(check_counter)
      if(length(which(check_counter==1))>0){
        MissData$n_0_reject[sss] <- MissData$n_0_reject[sss] +
          length(which(check_counter[1:which(check_counter==1)[1]]==0))
      } else{
        MissData$n_0_reject[sss] <- MissData$n_0_reject[sss] + MissData$n_batch_imp[sss]
      }
    }
    X_house <- X_house_sss_prop[which(check_counter==1)[1],]
    X_indiv <- matrix(comb_to_check[which(check_counter==1)[1],-c(1:length(individual_variable_index))],
                      byrow=TRUE,nrow=length(another_index)) #remove household head
    X_indiv[,relate_index] <- X_indiv[,relate_index] - 1 #recode relate back
    MissData$household[another_index,household_variable_index] <- rep(X_house,each=length(another_index))
    MissData$household[another_index,individual_variable_index] <- X_indiv
  }


return(MissData)
}




