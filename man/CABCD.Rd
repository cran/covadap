\name{CABCD}
\alias{CABCD}

\title{Covariate-Adjusted Biased Coin Design}

\description{
Implements the Covariate-adjusted Biased Coin Design by Baldi Antognini and Zagoraiou (2011), a stratified randomization procedure for two treatments A and B. The procedure works with qualitative covariates only.
}

\usage{
CABCD(data, a = 3, print.results = TRUE)
}

\arguments{
  \item{data}{a data frame or a matrix. Each row of \code{data} corresponds to the covariate profile of a patient. }
  \item{a}{(non-negative) design parameter determining the degree of randomness: \eqn{a = 0} gives the completely
randomized design; \eqn{a \rightarrow \infty} gives a deterministic design. The default value is set to 3. }
  \item{print.results}{ logical. If TRUE a summary of the results is printed. }
}

\details{
The function assigns patients to treatments A or B as described in Baldi Antognini and Zagoraiou (2011).

The parameter \code{a} determines the degree of randomness of the procedure.

At the end of the study, the imbalance measures reported are the loss of estimation precision as described in Atkinson (1982), the Mahalanobis distance and the overall imbalance, defined as the difference in the total number of patients assigned to treatment A and B. The strata imbalances measures report, for each stratum, the total number of patients assigned (\code{N.strata}), the number of patients assigned to A (\code{A.strata}) and the within-stratum imbalance (\code{D.strata}), calculated as \code{2*A.strata}-\code{N.strata}. The within-covariate imbalances report, for each level of each qualitative covariate, the difference in the number of patients assigned to A and B. See also Value.

}

\value{
It returns an object of \code{\link{class}} \code{"covadap"}, which is a list containing the following elements:
  \item{summary.info}{

        \code{Design} name of the design,

        \code{Sample_size} number of patients,

        \code{n_cov} number of covariates,

        \code{n_levels} number of levels of each covariate,

        \code{var_names} name of covariates and levels,

        \code{parameter_a} design parameter (see above).
        }

  \item{Assignments }{a vector with the treatment assignments.}

  \item{Imbalances.summary}{summary of overall imbalance measures at the end of the
                            study (\code{Loss} loss, \code{Mahal} Mahalanobis
                            distance, \code{overall.imb} difference in the total
                            number of patients assigned to A and B).}

  \item{Strata.measures}{a data frame containing for each possiblue stratum the corresponding
                          imbalances:
                                  \code{N.strata} is the total number of
                                  patients assigned to the stratum;
                                  \code{A.strata} is the total number of patients
                                  assigned to A  within the stratum;
                                  \code{D.strata} is the within-stratum imbalance,
                                  i.e. difference in the total number of patients
                                  assigned to A and B within the stratum.
                        }

   \item{Imbalances}{a list containing all the imbalance measures:

                    \code{Imb.measures} (\code{Loss} loss,
                    \code{Mahal} Mahalanobis distance),

                    \code{Overall.imb} difference in the total number of patients
                    assigned to A and B,

                    \code{Within.strata} within-stratum imbalance for all strata,

                    \code{Within.cov} within-covariate imbalance: difference in the
                    number of patients assigned to A and B for each level of each
                    qualitative covariate.
                  }

  \item{data}{the data provided in input.}


  \item{observed.strata}{a data frame with all the observed strata.}


}

\references{
Baldi Antognini A and Zagoraiou M. \emph{The covariate-adaptive biased coin design for balancing clinical trials in the presence of prognostic factors}. Biometrika, 2011, 98(3): 519-535.


Atkinson A. C. \emph{Optimum biased coin designs for sequential clinical trials with prognostic factors}. Biometrika, 1982, 69(1): 61-67.
}

\seealso{
\code{\link{CABCD.sim}} for allocating patients by simulating their covariate profiles.
}


\examples{
require(covadap)

# Create a sample dataset
df1 <- data.frame("gender" = sample(c("female", "male"), 100, TRUE, c(1 / 3, 2 / 3)),
                  "age" = sample(c("18-35", "36-50", ">50"), 100, TRUE),
                  "bloodpressure" = sample(c("normal", "high", "hyper"), 100, TRUE),
                   stringsAsFactors = TRUE)
# To just view a summary of the metrics of the design
CABCD(data = df1, a = 3)
# To view a summary
# and create a list containing all the metrics of the design
res <- CABCD(data = df1, a = 3)
res
}


\keyword{Covariate-adjusted biased coin design}
\keyword{Covariate-Adaptive randomization}
