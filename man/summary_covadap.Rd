\name{summary_covadap}
\alias{summary_covadap}

\title{
Summary of Covariate-Adaptive Designs
}

\description{
This function automatically recognizes the design implemented and provides a summary of the results.
}

\usage{
summary_covadap(res)
}

\arguments{
  \item{res}{
An object of class \code{"covadap"} or \code{"covadapsim"}  resulting from the application of a covariate-adaptive design.
}
}

\details{
When applied to an object of \code{\link{class}} \code{"covadap"}:
if at least one qualitative covariate is present, the function returns the within-covariate imbalances reporting, for each level of each qualitative covariate, the difference in the number of patients assigned to A and B.
If instead at least one quantitative covariate is present, the function returns the difference in means. For each quantitative covariate, the difference in the mean in group A and B is reported.

When applied to an object of \code{\link{class}} \code{"covadapsim"}, it reports the averages, in absolute value, of the imbalance measures, strata imbalances and within-covariate imbalances of the \code{nrep} trial replications according to the nature of the covariates.
}

\value{The form of the value returned by \code{summary_covadap} depends on the class of the argument provided (see Details). }

\examples{
#Create a sample dataset
df1 <- data.frame("gender" = sample(c("female", "male"), 100, TRUE, c(1 / 3, 2 / 3)),
                 "age" = sample(c("18-35", "36-50", ">50"), 100, TRUE),
                 "bloodpressure" = sample(c("normal", "high", "hyper"), 100, TRUE),
                 stringsAsFactors = TRUE)

res1 <- ECADE(data = df1, all.cat = TRUE,
                  alloc.function = "Efron", rho = 0.85)
summary_covadap(res1)


res2 <- ECADE.sim(data = df1, cov = NULL, n = NULL, all.cat = TRUE,
                  alloc.function = "Efron", rho = 0.85, nrep = 100)
summary_covadap(res2)

}
