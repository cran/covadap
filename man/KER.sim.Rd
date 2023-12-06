\name{KER.sim}
\alias{KER.sim}

\title{Simulations of the Covariate-Adaptive randomization by Ma and Hu}

\description{
Implements the Covariate-Adaptive randomization by Ma and Hu (2013) for assigning patients to two treatments A and B in order to minimize the distance
between the covariate distribution in the two treatment groups by simulating the covariate profile of each patient using an existing dataset or specifying number and levels of the covariates. The procedure works with qualitative and quantitative covariates.
}

\usage{
#With existing dataframe
KER.sim(data, covar = NULL, n = NULL, all.cat, nrep = 1000,
        p = 0.8, print.results = TRUE)
#With covariates
KER.sim(data = NULL, covar, n, all.cat, nrep = 1000,
        p = 0.8, print.results = TRUE)
}

\arguments{
   \item{data}{a data frame or a matrix. It can be a matrix only when \code{all.cat =
                TRUE}. Each row of \code{data} corresponds to the covariate profile of a
                patient.
                To be specified if covariate profiles should be sampled from an existing
                dataset provided in \code{data}.}
   \item{covar}{either a vector or a list to be specified only if \code{data = NULL}.
                If \code{all.cat = TRUE} is a vector with length equal to the number of
                covariates and elements equal to the number of levels for
                each covariate.
                Otherwise is a list containing \code{cat}, the list of the qualitative
                covariates with their level and \code{quant}, the list the quantitative
                covariates that are simulated from the normal distribution with
                the given mean
                and standard deviation.  }

   \item{n}{number of patients.}

    \item{all.cat}{logical. If all the covariates in \code{data} are qualitative
                    must be set equal to \code{TRUE}, otherwise must be set
                    equal to \code{FALSE}. }

  \item{nrep}{number of trial replications.}

  \item{p}{biasing probability for the Efron allocation function (\eqn{1/2 \leq p  \leq 1}). The default value is 0.8.}

  \item{print.results}{logical. If TRUE a summary of the results is printed. }
}
\details{
This function simulates \code{nrep} times a clinical study assigning patients to treatments A and B with the Efficient Covariate-Adaptive Design as described in Ma and Hu (see \code{\link{KER}}).

When \code{covar} is provided, the function finds all the possible combination of the levels of the covariates, i.e., the strata and, at each trial replication, the patients' covariate profiles are uniformly sampled within those strata. The specification of \code{covar} requires the specification of the number of patients \code{n}.

When \code{data} is provided, at each trial replication, the patients' covariate profiles are sampled from the observed strata with uniform distribution. In this case the number of patients equals the number of rows of \code{data}.

The summary printed when \code{print.results = TRUE} reports the averages, in absolute value, of the imbalance measures, strata imbalances and within-covariate imbalances of the \code{nrep} trial replications according to the nature of the covariates. See also \code{\link{KER}}.
}
\value{
It returns an object of \code{\link{class}} \code{"covadapsim"}, which is a list containing the following elements:
  \item{summary.info}{
                      \code{Design} name of the design.

                      \code{Sample_size} number of patients.

                      \code{n_cov} number of covariates.

                      \code{var_names} name of the covariates.

                      \code{n_quantitative_variables} number of
                      quantitative covariates. Is \code{NULL} if \code{all.cat = TRUE}.

                      \code{n_categorical_variables}  number of levels of each
                      covariate. Is \code{NULL} if \code{all.cat = TRUE}
                      or only quantitative covariates are present.

                      \code{n_levels} number of levels of each qualitative covariate.
                      Is \code{NULL} if only quantitative covariates are present.

                      \code{n.rep} number of replications.
                  }

  \item{Imbalances}{a list with the imbalance measures at the end of each simulated trial

            \code{Imb.measures} summary of overall imbalances
                    (\code{Loss} loss,
                     \code{Mahal} Mahalanobis distance,
                     \code{overall.imb} difference in the total number of patients
                     assigned to A and B).

            \code{within.imb}  within-covariate imbalance:
                difference in the number of patients assigned to A and B for
                each level of each qualitative covariate
                (is \code{NULL} if only quantitative covariates are present).

            \code{strata.imb} (only if \code{all.cat = TRUE}) the within-stratum
            imbalance (i.e. difference in the total number of patients assigned
            to A and B within the stratum).

            \code{strata.A} (only if \code{all.cat = TRUE}) is the total number of
            patients assigned to A within the stratum.

            \code{strata.N} (only if \code{all.cat = TRUE}) is the total number
            of patients assigned to each stratum.

            \code{diff_mean} (only if \code{all.cat = FALSE}) the difference in mean
                              in group A and B for each quantitative covariate.

            \code{obs.strata} (only if \code{all.cat = TRUE}) matrix of the possible
            strata.
     }

\item{out}{For each replication returns a list of the data provided in input (\code{data}) and the resulting assignments (\code{Assignment}). }


}
\references{
Ma Z and Hu F. \emph{Balancing continuous covariates based on Kernel densities}. Contemporary Clinical Trials, 2013, 34(2): 262-269.

}

\seealso{
See Also as \code{\link{KER}}.
}

\examples{
require(covadap)
# Here we set nrep = 50 for illustrative purposes,
# Set it equal to at least 5000 for more reliable Monte Carlo estimates.


### Implement with qualitative covariates (set all.cat = TRUE)
#### With an existing dataset
# Create a sample dataset with qualitative covariates
df1 <- data.frame("gender" = sample(c("female", "male"), 100, TRUE, c(1 / 3, 2 / 3)),
                  "age" = sample(c("18-35", "36-50", ">50"), 100, TRUE),
                  "bloodpressure" = sample(c("normal", "high", "hyper"), 100, TRUE),
                   stringsAsFactors = TRUE)
# To just view a summary of the metrics of the design
KER.sim(data = df1,  covar = NULL, n = NULL, all.cat = TRUE,
        p = 0.8, nrep = 50)
# To view a summary
# and create a list containing all the metrics of the design
res1 <- KER.sim(data = df1, covar = NULL, n = NULL, all.cat = TRUE,
                p = 0.8, nrep = 50)
#### By specifying the covariates
# e.g. two binary covariates and one with three levels and 100 patients
res2 <- KER.sim(data = NULL, covar = c(2,3,3), n = 100, all.cat = TRUE,
                p = 0.8, nrep = 50)



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
#### With an existing dataset
# select only column 5 and 6 of the sample dataset
# To just view a summary of the metrics of the design
KER.sim(data = ff1[,5:6], covar = NULL, n = NULL,  all.cat = FALSE, p = 0.8,
        nrep = 50)
# To view a summary
# and create a list containing all the metrics of the design
res3 <- KER.sim(data = ff1[,5:6],  covar = NULL, n = NULL, all.cat = FALSE,
                p = 0.8, nrep = 50)

#### By specifying the covariates
# BMI normally distributed with mean 26 and standard deviation 5
# cholesterol normally distributed with mean 200 and standard deviation 34
covar = list(quant = list(BMI = c(26, 5), cholesterol = c(200, 34)))
# To just view a summary of the metrics of the design
KER.sim(data = NULL, covar = covar, n = 100, all.cat = FALSE,
        p = 0.8, nrep = 50)
# To view a summary
# and create a list containing all the metrics of the design
res4 <- KER.sim(data = NULL, covar = covar, n = 100, all.cat = FALSE,
                  p = 0.8,  nrep = 50)

### With mixed covariates (set all.cat = FALSE)
#### With an existing dataset
# To just view a summary of the metrics of the design
KER.sim(data = ff1, covar = NULL, n = NULL, all.cat = FALSE, p = 0.8,
        nrep = 50)
# To view a summary
# and create a list containing all the metrics of the design
res5 <- KER.sim(data = ff1, covar = NULL, n = NULL, all.cat = FALSE,
                 p = 0.8, nrep = 50)
#### By specifying the covariates
# e.g. one qualitative covariate and 2 quantitative covariates:
# BMI normally distributed with mean 26 and standard deviation 5
# cholesterol normally distributed with mean 200 and standard deviation 34
# gender with levels M and F
covar = list(cat = list(gender = c("M", "F")),
             quant = list(BMI = c(26, 5), cholesterol = c(200, 34)))
# To just view a summary of the metrics of the design
KER.sim(data = NULL, covar = covar, n = 100, all.cat = FALSE,
        p = 0.8, nrep = 50)
# To view a summary and create a list containing all the metrics of the design
res6 <- KER.sim(data = NULL, covar = covar, n = 100, all.cat = FALSE,
                  p = 0.8, nrep = 50)

}

\keyword{Covariate-Adaptive randomization}
