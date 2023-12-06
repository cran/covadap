\name{KER}
\alias{KER}

\title{Covariate-Adaptive randomization by Ma and Hu}

\description{
Implements the Covariate-Adaptive randomization by Ma and Hu (2013) for assigning patients to two treatments A and B in order to minimize the distance
between the covariate distribution in the two treatment groups. The procedure works with qualitative and quantitative covariates.
}

\usage{
KER(data, all.cat, p = 0.8, print.results = TRUE)
}

\arguments{
  \item{data}{a data frame or a matrix. It can be a matrix only when \code{all.cat =
                TRUE}. Each row of \code{data} corresponds to the covariate profile of a
                patient.
                }

  \item{all.cat}{logical. If all the covariates in \code{data} are qualitative must be set equal to                 \code{TRUE}, otherwise must be set equal to \code{FALSE}.  }

  \item{p}{biasing probability for the Efron's allocation function (\eqn{1/2 \leq p  \leq 1}). The default value is 0.8. }


  \item{print.results}{ logical. If TRUE a summary of the results is printed. }
}

\details{
The function assigns patients to treatments A or B with the Covariate-Adaptive randomization based on kernel density estimation as described in Ma and Hu (2013).

This randomization procedure can be used when \code{data} contains only qualitative covariate, in this case set \code{all.cat = TRUE}, when \code{data} contains only quantitative covariates or when covariates of mixed nature are present, in these two latter cases set \code{all.cat = FALSE}.  The function's output is slighly different according to these three scenarios as described in Value.

The assignment probability to A of each patient is based on the Efron's allocation function (Efron, 1971) with biasing probability equal to \code{p}.

At the end of the study the imbalance measures reported are the loss of estimation precision as described in Atkinson (1982), the Mahalanobis distance and the overall imbalance, defined as the difference in the total number of patients assigned to treatment A and B.

Only when \code{all.cat = TRUE}, the function returns the strata imbalances measures, that report, for each stratum, the total number of patients assigned (\code{N.strata}), the number of patients assigned to A (\code{A.strata}) and the within-stratum imbalance (\code{D.strata}), calculated as \code{2*A.strata-N.strata}.

If at least one qualitative covariate is present, the function returns the within-covariate imbalances reporting, for each level of each qualitative covariate, the difference in the number of patients assigned to A and B.

If at least one quantitative covariate is present, the function returns the difference in means. For each quantitative covariate, is reported the difference in the mean in group A and B.

See Value for more details.
  }

\value{
It returns an object of \code{\link{class}} \code{"covadap"}, which is a list containing the following elements:

  \item{summary.info}{

                      \code{Design} name of the design.

                      \code{Sample_size} number of patients.

                      \code{n_cov} number of covariates.

                      \code{n_categorical_variables}  number of levels of each
                      covariate. Is \code{NULL} if
                      \code{all.cat = TRUE} or only quantitative covariates are present).

                      \code{n_levels} number of levels of each qualitative covariate.
                      Is \code{NULL} if only quantitative covariates are present.

                      \code{var_names} name of the covariates.

                      \code{cov_levels_names} levels of each qualitative covariate.
                      Is \code{NULL} if only quantitative covariates are present.

                      \code{n_quantitative_variables} number of quantitative
                            covariates.  Is \code{NULL} if \code{all.cat = TRUE}.
  }

  \item{Assignments }{a vector with the treatment assignments.}

  \item{Imbalances.summary}{summary of overall imbalance measures at the end
                            of the study (\code{Loss} loss, \code{Mahal}
                            Mahalanobis distance,
                            \code{overall.imb} difference in the total number
                            of patients assigned to A and B ).}

  \item{Strata.measures}{(only if \code{all.cat = TRUE}) a data frame
                        containing for each possiblue stratum the
                        corresponding imbalances: \code{N.strata} is the total number of
                        patients assigned to the stratum; \code{A.strata} is the total
                        number of patients assigned to A within the stratum;
                        \code{D.strata} is the within-stratum imbalance,
                        i.e. difference in the total number of patients
                        assigned to A and B within the stratum).
  }

  \item{Imbalances}{a list containing all the imbalance measures.

                    \code{Imb.measures} summary of overall imbalances
                                          (\code{Loss} loss,
                                          \code{Mahal} Mahalanobis distance,
                                          \code{overall.imb} difference in the total
                                          number of patients assigned to A and B).

                    \code{Within.strata} (only if \code{all.cat = TRUE}) within-stratum
                                          imbalance for all strata.

                    \code{Within.cov} within-covariate imbalance: difference in the
                    number of patients assigned to A and B for each level of each
                    qualitative covariate (is \code{NULL} if only quantitative
                    covariates are present).
                    }

    \item{data}{the data provided in input.}

    \item{diff_mean}{(only if \code{all.cat = FALSE}) the difference in mean of the
                      quantitative covariates in group A and B.}

    \item{observed.strata}{(only if \code{all.cat = TRUE}) a data frame with all the
                            observed strata.}

}

\references{
Ma Z and Hu F. \emph{Balancing continuous covariates based on Kernel densities}. Contemporary Clinical Trials, 2013, 34(2): 262-269.

Atkinson A. C. \emph{Optimum biased coin designs for sequential clinical trials with prognostic factors}. Biometrika, 1982, 69(1): 61-67.

Efron B, \emph{Forcing a sequential experiment to be balanced}. Biometrika, 1971, 58(3): 403-418.
}

\seealso{
See Also as \code{\link{KER.sim}} for allocating patients by simulating their covariate profiles.
}

\examples{
require(covadap)

### Implement with qualitative covariates (set all.cat = TRUE)
# Create a sample dataset with qualitative covariates
df1 <- data.frame("gender" = sample(c("female", "male"), 100, TRUE, c(1 / 3, 2 / 3)),
                  "age" = sample(c("18-35", "36-50", ">50"), 100, TRUE),
                  "bloodpressure" = sample(c("normal", "high", "hyper"), 100, TRUE),
                   stringsAsFactors = TRUE)
# To just view a summary of the metrics of the design
KER(data = df1,  all.cat = TRUE, p = 0.8, print.results = TRUE)
# To view a summary
# and create a list containing all the metrics of the design
res1 <- KER(data = df1,  all.cat = TRUE, p = 0.8, print.results = TRUE)
res1


### Implement with quantitative or mixed covariates
# Create a sample dataset with covariates of mixed nature
ff1 <- data.frame("gender" = sample(c("female", "male"), 100, TRUE, c(1 / 3, 2 / 3)),
                 "age" = sample(c("0-30", "30-50", ">50"), 100, TRUE),
                 "bloodpressure" = sample(c("normal", "high", "hypertension"), 10,
                  TRUE),
                 "smoke" = sample(c("yes", "no"), 100, TRUE, c(2 / 3, 1 / 3)),
                 "cholesterol" = round(rnorm(100, 200, 8),1),
                 "height" = rpois(100,160),
                  stringsAsFactors = TRUE)

### With quantitative covariates only (set all.cat = FALSE)
# select only column 5 and 6 of the sample dataset
# To just view a summary of the metrics of the design
KER(data = ff1[,5:6],  all.cat = FALSE, p = 0.8, print.results = TRUE)
# To view a summary
# and create a list containing all the metrics of the design
res2 <- KER(data = ff1[,5:6], p = 0.8, all.cat = FALSE, print.results = TRUE)
res2

### With mixed covariates
# In this case the user must set all.cat = FALSE
# To just view a summary of the metrics of the design
KER(data = ff1, all.cat = FALSE, p = 0.8, print.results = TRUE)
# To view a summary
# and create a list containing all the metrics of the design
res3 <- KER(data = ff1, all.cat = FALSE, p = 0.8, print.results = TRUE)
res3

}


\keyword{Covariate-Adaptive randomization}
