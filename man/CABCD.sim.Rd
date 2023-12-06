\name{CABCD.sim}
\alias{CABCD.sim}

\title{Simulations of the Covariate-Adjusted Biased Coin Design}

\description{
Implements the Covariate-adjusted Biased Coin Design by Baldi Antognini and Zagoraiou (2011) by simulating the covariate profile of each patient using an existing dataset or specifying number and levels of the covariates. The procedure works with qualitative covariates only.
}

\usage{
#With existing dataframe
CABCD.sim(data, covar = NULL, n = NULL, a = 3, nrep = 1000,
          print.results = TRUE)
#With covariates
CABCD.sim(data = NULL, covar, n, a = 3, nrep = 1000,
          print.results = TRUE)
}

\arguments{
   \item{data}{a data frame or a matrix. Each row of \code{data} corresponds to
              the covariate profile of a patient. To be specified if covariate
              profiles should be sampled from an existing dataset provided in
              \code{data}.}
   \item{covar}{either a vector or a list to be specified only if \code{data = NULL}.
              It could be a vector with length equal to the number
              of covariates and elements equal to the number of levels for each
              covariate. Otherwise it is a list containing the covariates with their
              levels (e.g. one covariate with two levels and one with three \code{covar
              = list(cov1 = c("lev1", "lev2"), cov2 = c("lev1", "lev2", "lev3"))}. }
  \item{n}{number of patients (to be specified only if \code{data = NULL}).}
  \item{a}{(non-negative) design parameter determining the degree of randomness:
          \eqn{a = 0} gives the completely randomized design;
          \eqn{a \rightarrow \infty} gives a deterministic design.
          The default value is set to 3.}
  \item{nrep}{number of trial replications.}
  \item{print.results}{logical. If TRUE a summary of the results is printed.}
}
\details{

This function simulates \code{nrep} times a clinical study assigning patients to treatments A and B with the Covariate-Adjusted Biased Coin Design (see \code{\link{CABCD}}).

When \code{covar} is provided, the function finds all the possible combination of the levels of the covariates, i.e., the strata and, at each trial replication, the patients' covariate profiles are uniformly sampled within those strata. The specification of \code{covat} requires the specification of the number of patients \code{n}.

When \code{data} is provided, at each trial replication, the patients' covariate profiles are sampled from the observed strata with uniform distribution. In this case the number of patients equals the number of rows of \code{data}.

The summary printed when \code{print.results = TRUE} reports the averages, in absolute value, of the imbalance measures, strata imbalances and within-covariate imbalances of the \code{nrep} trial replications. See also \code{\link{CABCD}}.


}



\value{
It returns an object of \code{\link{class}} \code{"covadapsim"}, which is a list containing the following elements:
  \item{summary.info}{

                      \code{Design} name of the design,

                      \code{Sample_size} number of patients,

                      \code{n_cov} number of covariates,

                      \code{n_levels} number of levels of each qualitative covariate,

                      \code{var_names} name of the covariates,

                      \code{n.rep} number of replications,

                      \code{parameter_a} design parameter (see above).

  }


\item{Imbalances}{a list with the imbalance measures at the end of each simulated trial:

            \code{Imb.measures} summary of overall imbalances  (\code{Loss} loss,
                                \code{Mahal} Mahalanobis distance, \code{overall.imb}
                                difference in the total number of patients
                                assigned to A and B),

            \code{within.imb}  within-covariate imbalance: difference in the number of
                                patients assigned to A and B for each level of each
                                qualitative covariate,

            \code{strata.imb}  the within-stratum imbalance (i.e. difference in the
                                total number of patients assigned to A and B within
                                the stratum),

            \code{strata.A}    total number of patients assigned to A within the
                                stratum,

            \code{strata.N}  total number of patients assigned to each stratum,

            \code{obs.strata} matrix of the possible strata.
     }


\item{out}{For each replication returns a list of the data provided in input (\code{data}) and the resulting assignments (\code{Assignment}). }

}
\references{
Baldi Antognini A and Zagoraiou M. \emph{The covariate-adaptive biased coin design for balancing clinical trials in the presence of prognostic factors}. Biometrika, 2011, 98(3): 519-535.
}

\seealso{
See Also  \code{\link{CABCD}}.
}

\examples{
require(covadap)
# Here we set nrep = 100 for illustrative purposes,
# Set it equal to at least 5000 for more reliable Monte Carlo estimates.


### With existing dataframe
df1 <- data.frame("gender" = sample(c("female", "male"), 100, TRUE, c(1 / 3, 2 / 3)),
                 "age" = sample(c("18-35", "36-50", ">50"), 100, TRUE),
                 "bloodpressure" = sample(c("normal", "high", "hyper"), 100, TRUE),
                 stringsAsFactors = TRUE)
# Simulate the design
res1 <- CABCD.sim(data = df1, n = NULL, a = 3, nrep = 100)
### With covariates
# e.g. two binary covariates and one with three levels and 100 patients
res2 <- CABCD.sim(covar = c(2,2,3), n = 100, a = 3,  nrep = 100)
}

\keyword{Covariate-adjusted biased coin design}
\keyword{Covariate-Adaptive randomization}

