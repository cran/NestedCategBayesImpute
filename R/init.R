initData <- function(household, individual_variable_index, household_variable_index) {
  orig <- list()
  orig$origdata <- household
  if (which(names(household) == "Hhindex") != 1) {
    stop("The first column of input has to be 'Hhindex'")
  }
  orig$n_i <- as.data.frame(table(household[,'Hhindex']))$Freq
  orig$n <- length(orig$n_i)

  HHrowIndex <- c(1, cumsum(orig[["n_i"]])+1)

  #household level data is in household_variable_index
  orig$HHdataorigT <- t(household[HHrowIndex[1:orig$n],household_variable_index])
  orig$HHserial <- household[,"Hhindex"]

  orig$n_individuals <- dim(household)[1]
  orig$n_individuals_real <- dim(household)[1] + orig$n #the real number of individuals

  #orig$p holds the number of individual level variables
  orig$p <- length(individual_variable_index)
  #levels for each variable in the model. This might be different from the sub-sampled data in use
  orig$d <- rep(0,orig$p)
  for (i in 1:length(orig$d)) {
    orig$d[i] <- max(household[,individual_variable_index[i]])
  }
  orig$dataT <- t(household[,individual_variable_index])
  orig$maxd <- max(orig$d)

  counts <- as.data.frame(table(orig$n_i))
  orig$n_star_h <- counts[order(counts[,1]),2]
  return(orig)
}


initParameters <- function(data,hyper,HHhead_at_group_level) {
  para <- list()
  para$alpha <- 1 #hyperparameters for stick-breaking weights
  para$beta <- 1

  #intilize phi
  para$phi <- matrix(0, nrow  = data$maxd*data$p, ncol = hyper$FF*hyper$SS) #cell probabilities
  phi_1 <- matrix(0, nrow = data$maxd, ncol = data$p)
  for (i in 1:data$p) {
    for (j in 1:data$d[i]) {
      phi_1[j,i] <-sum(data$dataT[i,]==j)/data$n_individuals
    }
  }
  phi_1 <- as.vector(phi_1)
  for (i in 1:dim(para$phi)[2]) {
    para$phi[,i] <- phi_1
  }

  para$HHdata_all <- data$HHdataorigT
  if (!HHhead_at_group_level) {
    para$HHdata_all[2,] <- para$HHdata_all[2,] - 1 #household size
  }
  #initialize lambda
  para$lambda <- list()
  for (i in 1:length(hyper$dHH)) {
    lambda <- matrix(0, nrow = hyper$FF,ncol = hyper$dHH[i])
    for (j in 1:hyper$dHH[i]) {
      lambda[,j] <- sum(para$HHdata_all[i,]==j) / data$n
    }
    para$lambda[[i]] <- lambda
  }

  para$u <- c(rbeta(hyper$FF-1, 1,para$alpha),1)
  para$pi <- para$u * cumprod(c(1,1.0-para$u[1:hyper$FF-1]))

  ones <- matrix(1.0, hyper$FF,1)
  para$v <- c(rbeta(hyper$FF * (hyper$SS-1), 1, para$beta),ones)
  dim(para$v) <- c(hyper$FF, hyper$SS)

  para$omega <- matrix(0, nrow = hyper$FF,ncol = hyper$SS)
  for (i in 1:hyper$FF) {
    v1 <- para$v[i,]
    para$omega[i,]  <- v1 * cumprod(c(1, 1-v1[1:hyper$SS-1]))
  }

  return(para)
}

initOutput <- function(data,hyper,mc) {
  output <- list()
  output$alphaout <- matrix(0,nrow = mc$eff.sam,ncol = 1)
  output$betaout <- matrix(0, nrow = mc$eff.sam,ncol = 1)
  output$piout <- matrix(0, mc$eff.sam,hyper$FF)
  output$omegaout <- array(0, dim=c(mc$eff.sam,hyper$FF,hyper$SS))
  output$nout <- matrix(0,nrow = mc$nrun,ncol = 1)
  output$extrasize <- matrix(0,nrow = mc$nrun,ncol = length(data$n_star_h))
  output$F_occupied <- matrix(0, nrow = mc$eff.sam, ncol = 1)
  output$S_occupied_max <- matrix(0, nrow = mc$eff.sam, ncol = 1)
  output$elapsed_time <-  matrix(0,nrow = mc$nrun,ncol = 1)
  output$newphiout <- array(0, dim=c(mc$eff.sam,data$maxd*data$p,hyper$FF*hyper$SS))
  output$lambdaout = list()
  for (i in 1:length(hyper$dHH)) {
    output$lambdaout[[i]] = array(0, dim=c(mc$eff.sam,hyper$FF, hyper$dHH[i]))
  }
  return(output)
}
