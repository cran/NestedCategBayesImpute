.GenerateData <- function(hh_size,lambda, omega, phi,pi, d, total.batch,possiblehhcount,howmany,synindex,HHhead_at_group_level) {
  #save(hh_size,lambda, omega, phi,pi, d, total.batch,possiblehhcount,howmany,synindex, file = "last2.RData")
  #return(NULL)
  Individuals_extra <- list()
  G_extra <- list()
  HHData_extra <- list()
  batch.index <- 0
  n_possible_household <- 0
  p <- length(d)
  synIndividuals <- list()

  while (n_possible_household< possiblehhcount) {
    batch.index <- batch.index + 1
    #print(batch.index)
    #generate a batch of 10K household
    if (HHhead_at_group_level) {
      data_to_check <- samplehouseholds_HHhead_at_group_level(phi,omega, pi, d, lambda,batch.index+total.batch, howmany,hh_size)
    } else {
      data_to_check <- samplehouseholds(phi,omega, pi, d, lambda,batch.index+total.batch, howmany,hh_size)
    }

    #impossible household

    #data_to_check_old_format <- ConvertDataForward(data_to_check,hh_size,p,lambda)
    #print("do_check")
    if (HHhead_at_group_level) {
      checked.households <- checkconstraints_HHhead_at_group_level(data_to_check,possiblehhcount-n_possible_household, hh_size)
    } else {
      checked.households <- checkconstraints(data_to_check,possiblehhcount-n_possible_household, hh_size)
    }
    n_possible_household <- n_possible_household + checked.households$possible
    if (length(checked.households$Households) > 0) {
      #checked.households$Households <- t(ConvertDataBack(t(checked.households$Households),hh_size,p,lambda))
      #print("do convert")
      Individuals_extra[[batch.index]] <- households2individuals(checked.households$Households,hh_size)
      DIM <- p + length(lambda) + 1
      G_extra[[batch.index]] <- checked.households$Households[hh_size * DIM +1,]
      HHData_extra[[batch.index]] <- checked.households$Households[(p+3): DIM,]
      #print("done convert")
    } else {
      print("no batch")
    }


    if (synindex > 0) {
      if (length(checked.households$synHouseholds) > 0) {
        #checked.households$synHouseholds <- t(ConvertDataBack(t(checked.households$synHouseholds),hh_size,p,lambda))
        synIndividuals[[batch.index]]  <- households2individuals(checked.households$synHouseholds,hh_size)
      }
    }
  }

  Individuals_extra <- do.call(cbind, Individuals_extra)
  G_extra <- unlist(G_extra)
  if (HHhead_at_group_level) {
    HHData_extra <- do.call(cbind, HHData_extra)
  } else {
    HHData_extra <- unlist(HHData_extra)
  }

  if (synindex > 0) {
    synIndividuals <- do.call(cbind, synIndividuals)
  }
  batch.index <- batch.index + total.batch
  return(list(Individuals_extra = Individuals_extra,
              G_extra = G_extra,
              HHData_extra = HHData_extra,
              synIndividuals = synIndividuals,
              batch.index = batch.index))
}

GetImpossibleHouseholds <- function(d,n_star_h,lambda,omega,phi,pi,howmany,n,synindex,HHhead_at_group_level) {
  #save(d,n_star_h,lambda,omega,phi,pi,howmany,n,synindex, file = "last.RData")
  #return(NULL)
  cumsize <- 0
  hh_size_new <-  matrix(0,nrow = length(n_star_h), ncol = 1)
  hh_index <- list()
  ImpossibleIndividuals <- list()
  G_extra <- list()
  HHdata_extra <- list()
  synIndividuals_all <- list()

  ##
  total.batch <- 0
  for (hh_size in  1:(length(n_star_h))) {
    if (HHhead_at_group_level) {
      hh_size_real <- hh_size
    } else {
      hh_size_real <- hh_size + 1
    }
    batch <- .GenerateData(hh_size_real,lambda, omega, phi,pi, d, total.batch,n_star_h[hh_size],howmany,synindex,HHhead_at_group_level)

    hh_size_new[hh_size] <- length(batch$G_extra)
    hh_index[[hh_size]] <- cumsize + rep(1:hh_size_new[hh_size], each = hh_size_real)
    cumsize <- cumsize + hh_size_new[hh_size]
    ImpossibleIndividuals[[hh_size]] <- batch$Individuals_extra
    G_extra[[hh_size]] <-  batch$G_extra
    HHdata_extra[[hh_size]] <- rbind(batch$HHData_extra,rep(hh_size, times = hh_size_new[hh_size])) # length(lambda) by ...
    if (synindex > 0) {
      synIndividuals_all[[hh_size]] <- batch$synIndividuals
    }
    total.batch <- batch$batch.index
  }

  ##
  hh_index <- unlist(hh_index)
  ImpossibleIndividuals <- do.call(cbind, ImpossibleIndividuals)
  G_extra <- unlist(G_extra)
  HHdata_extra <- do.call(cbind, HHdata_extra)
  if (synindex > 0) {
    synIndividuals_all <- do.call(cbind, synIndividuals_all)
  }

  ImpossibleIndividuals[1,] <- n + hh_index
  DIM <- dim(ImpossibleIndividuals)[1] - 2
  IndividualData_extra <- ImpossibleIndividuals[1:DIM,]
  G_Individuals_and_M_extra <- ImpossibleIndividuals[(DIM+1):(DIM+2),]
  return(list(G_Individuals_and_M_extra = G_Individuals_and_M_extra,
              G_extra = G_extra,
              IndividualData_extra = IndividualData_extra,
              HHdata_extra = HHdata_extra,
              hh_size_new = hh_size_new,
              synIndividuals_all = synIndividuals_all))
}
