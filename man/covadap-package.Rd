
\encoding{latin1}
\name{covadap-package}
\alias{covadap}
\docType{package}
\title{covadap: Implements Covariate-Adaptive Randomization procedures}

\description{
Implementing seven Covariate-Adaptive Randomization to assign patients to two treatments. Three of these procedures can also accommodate quantitative and mixed covariates. Given a set of covariates, the user can generate a single sequence of allocations or replicate the design multiple times by simulating the patients' covariate profiles. At the end, an extensive assessment of the performance of the randomization procedures is provided, calculating several imbalance measures.
}

\author{R. Frieri \email{rosamarie.frieri2@unibo.it}, M. Novelli \email{m.novelli@unibo.it}}

\section{Acknowledgement}{
This work was supported by the EU funding within the NextGenerationEU PRIN2022 \emph{Optimal and adaptive designs for modern medical experimentation} (2022TRB44L).}

\references{
Atkinson A C. \emph{Optimum biased coin designs for sequential clinical trials with prognostic factors}. Biometrika, 1982, 69(1): 61-67.

Baldi Antognini A, Frieri R, Zagoraiou M, Novelli M. \emph{The Efficient Covariate-Adaptive Design for high-order balancing of quantitative and qualitative covariates}. Statistical Papers, 2022.

Baldi Antognini A and Zagoraiou M. \emph{The covariate-adaptive biased coin design for balancing clinical trials in the presence of prognostic factors}. Biometrika, 2011, 98(3): 519-535.

Efron B, \emph{Forcing a sequential experiment to be balanced}. Biometrika, 1971, 58(3): 403-418.

Hu Y, Hu F. \emph{Asymptotic properties of covariate-adaptive randomization}. The Annals of Statistics, 2012, 40(3): 1974-1815.

Ma Z and Hu F. \emph{Balancing continuous covariates based on Kernel densities}. Contemporary Clinical Trials, 2013, 34(2): 262-269.

Pocock S J, Simon R. \emph{Sequential treatment assignment with balancing for prognostic factors in the controlled clinical trial}. Biometrics, 1975, 31(1): 103-115.

Soares F, Wu CFJ.  \emph{Some restricted randomization rules in sequential designs}. Communications in Statistics
Theory and Methods 1983, 12: 2017-2034.
}
\keyword{Covariate-Adaptive randomization}
\keyword{Quantitative covariates}
\keyword{Mixed covriates}
\keyword{covadap}

\seealso{
\code{\link{CABCD}, \link{HuHu}, \link{PocSim}, \link{BSD}, \link{DABCD}, \link{ECADE}, \link{KER}.}
}

