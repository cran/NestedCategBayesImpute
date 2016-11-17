
###################################################################################
###################################### START ######################################
###################################################################################

### Empty environment and load required libraries
rm(list = ls())
library(NestedCategBayesImpute)
library(dplyr)


### Set indicator for whether of not to move the household head
### Also set indicator for the weighting/capping option
HHhead_at_group_level <- TRUE #set to TRUE to move household head to the group level
weight_option <- TRUE #set to TRUE for weighting/capping option. If TRUE, must supply weights


### Use data included in package; prepare data and specify variable indexes
if (HHhead_at_group_level) {
  orig.file <- system.file("extdata","origdata_newFormat.txt",package="NestedCategBayesImpute")
  orig.data <- read.table(orig.file,header = TRUE, sep = " ")
  orig.data$relate <- orig.data$relate - 1L #recode relate variable to 11 levels
  household.size <- as.data.frame(table(orig.data$Hhindex))
  household.size[,1] <- as.numeric(household.size[,1])
  names(household.size) <- c("Hhindex", 'householdsize')
  household <- orig.data %>% inner_join(household.size)

  individual_variable_index = c(3:7)
  household_variable_index = c(8:13) #make sure the last column represents household size
} else {
  orig.file <- system.file("extdata","origdata_oldFormat.txt",package="NestedCategBayesImpute")
  orig.data <- read.table(orig.file,header = TRUE, sep = " ")
  orig.data$Hhindex
  household.size <- as.data.frame(table(orig.data$Hhindex))
  household.size[,1] <- as.numeric(household.size[,1])
  names(household.size) <- c("Hhindex", 'householdsize')
  household <- orig.data %>% inner_join(household.size)

  individual_variable_index = c(3:7)
  household_variable_index = c(8,9) #make sure the last column represents household size
}


### Initialize the input data structure
orig <- initData(household,individual_variable_index,household_variable_index)


### Check first few lines of data
#head(household)


### Supply weights; one for each household size
if(weight_option){
  struc_weight <- c(1/2,1/2,1/3,1/3,1/3) #must be ordered & no household size must be excluded
} else {
  struc_weight <- rep(1,length(orig$n_star_h)) #just a dummy column of ones if struc_weight=FALSE
}


### Set mcmc parameters
mc <- list(nrun = 10000, burn = 5000, thin = 5)
mc$eff.sam <- (mc$nrun-mc$burn)/mc$thin


### Set number of categories for each household level variable
dHH <- rep(0,length(household_variable_index))
for (i in 1:length(dHH)) {
  dHH[i] <- max(household[,household_variable_index[i]])
  if (i == length(dHH) & !HHhead_at_group_level) {
    dHH[length(dHH)] <- dHH[length(dHH)] - 1 #Household head within household assumes that the household size starts from 2
  }
}


### Set hyper parameters
#aa & ab are gamma hyperparameters for alpha while ba & bb are gamma hyperparameters for beta
#blocksize is the number of impossible households to sample at once (we use batch sampling to speed up mcmc)
#FF is the max number of group-level latent classes
#SS is the max number of individual-level classes
hyper <- list(FF=20 , SS=15, aa=0.25, ab=0.25, ba=0.25,bb=0.25,dHH = dHH, blocksize = 10000)


### Initialize parameters and output
para <- initParameters(orig,hyper,HHhead_at_group_level)
output <- initOutput(orig,hyper,mc)


### Set number of synthetic data and the mcmc indexes for them
mm <- 5
synindex <- sort(sample(seq((mc$burn +1),mc$nrun,by=mc$thin),mm,replace=FALSE))


### Run model
ModelResults <- RunModel(orig,mc,hyper,para,output,synindex,individual_variable_index,household_variable_index,
                         HHhead_at_group_level,weight_option,struc_weight)


### View first few lines of the first synthetic data.
head((ModelResults$synData)[[1]]) # Remember that the relate variable has been recoded to 11 levels


### Some posterior summaries and plots
library(coda)
names(ModelResults$output)
dim(ModelResults$output$alphaout)
alpha_output <- mcmc(ModelResults$output$alphaout)
plot(alpha_output)
summary(alpha_output)

dim(ModelResults$output$betaout)
beta_output <- mcmc(ModelResults$output$betaout)
plot(beta_output)
summary(beta_output)

dim(ModelResults$output$nout)
total_households <-mcmc(ModelResults$output$nout)
plot(total_households)
summary(total_households)

dim(ModelResults$output$extrasize)
impossible_households <-mcmc(ModelResults$output$extrasize)
plot(impossible_households)
summary(impossible_households)

dim(ModelResults$output$elapsed_time)
time_per_iteration <-mcmc(ModelResults$output$elapsed_time)
plot(time_per_iteration)
summary(time_per_iteration)

dim(ModelResults$output$F_occupied)
F_occupied <-mcmc(ModelResults$output$F_occupied)
plot(F_occupied)
summary(F_occupied)

dim(ModelResults$output$S_occupied_max)
S_occupied_max <-mcmc(ModelResults$output$S_occupied_max)
plot(S_occupied_max)
summary(S_occupied_max)

###################################################################################
####################################### END #######################################
###################################################################################







