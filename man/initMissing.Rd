\name{initMissing}
\alias{initMissing}
\title{
Initilize the misising data structure from input data
}
\description{
Initilize the misising data structure from input data
}
\usage{
initMissing(household, individual_variable_index, household_variable_index,
            miss_batch, struc_zero_variables)
}
\arguments{
  \item{household}{A data matrix for the input household data.}
   \item{individual_variable_index}{Vector of column indexes for the individual-level variables.}
  \item{household_variable_index}{Vector of column indexes for the group-level variables.}
  \item{miss_batch}{initial number of batches to sample for each household with missing data.}
  \item{struc_zero_variables}{column indexes for the variables that define structural zeros like age and relate (including those for the household head).}
}
