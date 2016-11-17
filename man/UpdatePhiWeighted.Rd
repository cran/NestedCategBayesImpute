\name{UpdatePhiWeighted}
\alias{UpdatePhiWeighted}
\title{
Update phi.
}
\description{
Update phi -- the matrix of posterior samples for the individual-level probabilities for each individual-level variable by each pair of group-level and individual-level latent classes -- when the weighting/capping option is used. The weighting options allows capping the number of impossible households to sample and re-weight the multinomial counts within each latent class back to the expected truth.
}
\usage{
UpdatePhiWeighted(IndividualData_all, M_all, FF, SS, p, d, maxd,
                    individual_variable_index, struc_weight)
}
\arguments{
  \item{IndividualData_all}{
  Data matrix for the individual-level data from both the original data and the sampled impossible households.
}
  \item{M_all}{
 A vector of both the household-level and individual-level latent class indexes for all households both in the original data and the sampled impossible households.
}
  \item{FF}{
  Maximum number of household-level latent classes allowed.
}
  \item{SS}{
  Maximum number of individual-level latent classes allowed.
}
  \item{p}{
  Number of individual-level variables.
}
  \item{d}{
  A vector for the number of levels of each individual-level variable.
}
  \item{maxd}{
  Maximum value in \verb{d}.
}
\item{individual_variable_index}{Vector of column indexes for the individual-level variables.}
\item{struc_weight}{A vector of weights by household sizes used in capping the number of sampled impossible households.}
}
\details{
Function for obtaining a posterior sample of phi when the weighting/capping option is used.
}
\value{
Updated (posterior) value for phi based on the corresponding full conditional distribution.
}

\author{
Quanli Wang, Olanrewaju Akande
}

\keyword{ sampler }
