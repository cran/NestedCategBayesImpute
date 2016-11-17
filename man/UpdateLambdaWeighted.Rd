\name{UpdateLambdaWeighted}
\alias{UpdateLambdaWeighted}
\title{
Update lambda.
}
\description{
Update lambda -- the list of matrices of the group-level probabilities for each group-level variable by the group-level latent classes -- when the weighting/capping option is used. The weighting options allows capping the number of impossible households to sample and re-weight the multinomial counts within each latent class back to the expected truth. Each matrix in the list is for each group-level variable.
}
\usage{
UpdateLambdaWeighted(dHH, FF, G_all, HHdata_all,struc_weight)
}
\arguments{
  \item{dHH}{
  A vector containing the number of levels for each household-level variable.
}
  \item{FF}{
  Maximum number of household-level latent classes allowed.
}
  \item{G_all}{
  A vector of the household-level latent class indexes for all households both in the original data and the sampled impossible households.
}
  \item{HHdata_all}{
  Data matrix for the household-level data from both the original data and the sampled impossible households.
}
\item{struc_weight}{A vector of weights by household sizes used in capping the number of sampled impossible households.}
}
\details{
Function for obtaining a posterior sample of lambda when the weighting/capping option is used.
}
\value{
Updated (posterior) value for lambda based on the corresponding full conditional distribution.
}

\author{
Quanli Wang, Olanrewaju Akande
}

\keyword{ sampler }

