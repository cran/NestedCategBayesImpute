

RunModel <- function(orig,mc,hyper,para,output,synindex,individual_variable_index,household_variable_index,
                     HHhead_at_group_level,weight_option,struc_weight){
  synData <- list()
  if(weight_option){
    struc_weight <- c(1,struc_weight) #add 1 for the weight of the observed data
    G_all_weighted <- vector("list",length(struc_weight))
    HHdata_all_weighted <- vector("list",length(struc_weight))
    IndividualData_all_weighted <- vector("list",length(struc_weight))
    M_all_weighted <- vector("list",length(struc_weight))
  }

  for (i in 1:mc$nrun) {
    cat(paste("iteration ", i,"\n", sep = ""))
    t <- proc.time()

    G_household <- sampleG(para$phi,orig$dataT,para$omega,para$pi,orig$n_i,t(para$HHdata_all[,1:orig$n]),para$lambda)

    M <- sampleM(para$phi,orig$dataT,para$omega,G_household$G,orig$HHserial)

    if(weight_option){
      data.extra <- GetImpossibleHouseholds(orig$d,ceiling(orig$n_star_h*struc_weight[-1]),para$lambda,para$omega,para$phi,
                                            para$pi,hyper$blocksize,orig$n,is.element(i,synindex),HHhead_at_group_level)
      para$hh_size_new <- as.vector(data.extra$hh_size_new)
      DIM <- dim(data.extra$IndividualData_extra)[1]
      if (is.element(i,synindex)) {
        forsynData <- GetImpossibleHouseholds(orig$d,orig$n_star_h,para$lambda,para$omega,para$phi,
                                              para$pi,hyper$blocksize,orig$n,is.element(i,synindex),HHhead_at_group_level) #synthetic data
        synData[[which(synindex ==i)]] <- t(forsynData$synIndividuals_all[1:DIM,])
        colnames(synData[[which(synindex ==i)]]) <- colnames(orig$origdata)[-ncol(orig$origdata)]
      }

      #combine data and indicators -- use lists for weighting
      n_i_extra <- as.data.frame(table(data.extra$IndividualData_extra[1,]))$Freq
      n_i_extra_index <- rep(n_i_extra,n_i_extra)
      G_all_weighted[[1]] <- G_household$G
      temp_temp <- orig$HHdataorigT
      if (!HHhead_at_group_level) {
        temp_temp[2,] <- temp_temp[2,] -1
      }
      HHdata_all_weighted[[1]] <- temp_temp
      IndividualData_all_weighted[[1]] <- t(orig$origdata[,1:DIM])
      M_all_weighted[[1]] <- rbind(G_household$G_Individuals,M)
      M_all_vector <- cbind(M_all_weighted[[1]],data.extra$G_Individuals_and_M_extra)

      if (!HHhead_at_group_level){
        for(w_i in 2:length(struc_weight)){
          G_all_weighted[[w_i]] <- data.extra$G_extra[(n_i_extra == w_i)]
          HHdata_all_weighted[[w_i]] <- data.extra$HHdata_extra[,(n_i_extra == w_i)]
          IndividualData_all_weighted[[w_i]] <- data.extra$IndividualData_extra[,(n_i_extra_index == w_i)]
          M_all_weighted[[w_i]] <- data.extra$G_Individuals_and_M_extra[,(n_i_extra_index == w_i)]
        }
      } else{
        for(w_i in 2:length(struc_weight)){
          G_all_weighted[[w_i]] <- data.extra$G_extra[(n_i_extra == (w_i-1))]
          HHdata_all_weighted[[w_i]] <- data.extra$HHdata_extra[,(n_i_extra == (w_i-1))]
          IndividualData_all_weighted[[w_i]] <- data.extra$IndividualData_extra[,(n_i_extra_index == (w_i-1))]
          M_all_weighted[[w_i]] <- data.extra$G_Individuals_and_M_extra[,(n_i_extra_index == (w_i-1))]
        }
      }


      # update phi
      para$phi <- UpdatePhiWeighted(IndividualData_all_weighted,M_all_weighted,
                                    hyper$FF,hyper$SS,orig$p,orig$d,orig$maxd,individual_variable_index,struc_weight)
      #update omega
      Omega <- UpdateOmegaWeighted(para$beta,M_all_weighted, hyper$FF, hyper$SS,struc_weight)
      para$omega <- Omega$omega
      para$v <- Omega$v

      # update lambda
      para$lambda <- UpdateLambdaWeighted(hyper$dHH,hyper$FF,G_all_weighted,HHdata_all_weighted,struc_weight)

      # update pi
      Pi <- UpdatePiWeighted(para$alpha,G_all_weighted,hyper$FF,struc_weight)
      para$pi <- Pi$pi
      para$u <- Pi$u

    } else {
      data.extra <- GetImpossibleHouseholds(orig$d,orig$n_star_h,para$lambda,para$omega,para$phi,
                                            para$pi,hyper$blocksize,orig$n,is.element(i,synindex),HHhead_at_group_level)
      para$hh_size_new <- as.vector(data.extra$hh_size_new)
      DIM <- dim(data.extra$IndividualData_extra)[1]
      if (is.element(i,synindex)) {
        synData[[which(synindex ==i)]] <- t(data.extra$synIndividuals_all[1:DIM,])
      }

      #combine data and indicators
      para$G_all <- c(G_household$G, data.extra$G_extra)
      para$HHdata_all <- orig$HHdataorigT
      if (!HHhead_at_group_level) {
        para$HHdata_all[2,] <- para$HHdata_all[2,] -1
      }
      para$HHdata_all <- cbind(para$HHdata_all,data.extra$HHdata_extra)
      para$IndividualData_all <- cbind(t(orig$origdata[,1:DIM]),data.extra$IndividualData_extra)

      #row 1 for FF groups and row 2 for SS groups
      temp <- rbind(G_household$G_Individuals,M)
      para$M_all  <- cbind(temp,data.extra$G_Individuals_and_M_extra)

      # update phi
      para$phi <- UpdatePhi(para$IndividualData_all,para$M_all,
                            hyper$FF,hyper$SS,orig$p,orig$d,orig$maxd,individual_variable_index)

      #update Omega
      Omega <- UpdateOmega(para$beta,para$M_all, hyper$FF, hyper$SS)
      para$omega <- Omega$omega
      para$v <- Omega$v

      # update lambda
      para$lambda <- UpdateLambda(hyper$dHH,hyper$FF,para$G_all,para$HHdata_all)

      # update pi
      Pi <- UpdatePi(para$alpha,para$G_all,hyper$FF)
      para$pi <- Pi$pi
      para$u <- Pi$u
    }

    #update alpha
    para$alpha <- UpdateAlpha(hyper$aa,hyper$ab,para$u)

    #update beta
    para$beta <- UpdateBeta(hyper$ba,hyper$bb,para$v)



    #post save
    if(weight_option){
      G_all <- unlist(G_all_weighted)
      rep_G_all <- t(M_all_vector[1,])
      M_all <- t(M_all_vector[2,])
    } else {
      G_all <- para$G_all
      rep_G_all <- t(para$M_all[1,])
      M_all <- t(para$M_all[2,])
    }
    S_occup <- NULL
    for(occ in sort(unique(G_all))){
      S_occup <- rbind(S_occup,dim(table(rep_G_all[which(rep_G_all==occ)],M_all[which(rep_G_all==occ)]))[2])
    }
    cat(paste("number of occupied household classes is ", length(unique(G_all)), "\n", sep = ''))
    cat(paste("max number of occupied individual classes is ", max(S_occup), "\n", sep = ''))

    total_household <- sum(c(orig$n,para$hh_size_new))
    if(weight_option){
      cat(paste("total number of households (capped) sampled is ", total_household, "\n", sep = ''))
      est_total_household <- sum(c(orig$n,para$hh_size_new)/struc_weight)
      cat(paste("true (estimated) total number of households is ", est_total_household, "\n", sep = ''))
    } else {
      cat(paste("total number of households sampled is ", total_household, "\n", sep = ''))
    }
    cat(paste("elapsed time = ", (proc.time() - t)[["elapsed"]], "\n\n", sep = ' '))

    output$nout[i] <- total_household
    output$extrasize[i,] <- para$hh_size_new
    output$elapsed_time[i] <- (proc.time() - t)[["elapsed"]]


    if (i %% mc$thin == 0 && i > mc$burn)  {
      index <- (i-mc$burn)/mc$thin
      output$piout[index,] <- para$pi
      output$omegaout[index,,] <- para$omega
      output$newphiout[index,,] <- para$phi
      for (i in 1:length(hyper$dHH)) {
        output$lambdaout[[i]][index,,] = para$lambda[[i]]
      }
      output$F_occupied[index] <- length(unique(G_all))
      output$S_occupied_max[index] <- max(S_occup)
      output$alphaout[index] <- para$alpha
      output$betaout[index] <- para$beta
    }
  }

  return(list(synData=synData,output=output))
}
