#Covariate-Adaptive randomization by Hu and Hu
HuHu <- function(data, p = 0.85, omega = NULL, print.results = TRUE){
  if(missing(data)){
    stop("data must be provided")
  }
  if(!is.numeric(p)){
    stop("p must be a postive number in [0.5, 1]")
  }
  if(p<0.5 | p >1){
    stop("p must be a postive number between in [0.5, 1]")
  }

  dr = data_preproc(data)
  n = dr$n
  if(n < 3){
    stop("at least 3 observations are needed")
  }
  delta = numeric(n)
  A.strata = N.strata = D.strata = numeric(dr$n.strata)
  strata = dr$strata
  obs.strata = dr$obs.strata
  n_cov = dr$n_cov
  n_levels = dr$n_levels
  if (is.null(omega)) {
    omega = rep(1/(n_cov + 2), n_cov + 2)
  }else if (length(omega) != (2 + n_cov)) {
    stop("Length of omega must equal to ncols(data) + 2 !")
  }else{
    omega = abs(omega)/sum(abs(omega))
  }

  delta[1] = rbinom(1,1,.5)
  sel = strata[1]
  N.strata[sel] = N.strata[sel] + 1
  A.strata[sel] = A.strata[sel] + delta[1]
  D.strata = 2 * A.strata - N.strata

  num_strata = dr$num_strata
  obs.num_strata = unique(dr$num_strata)
  map.strata = cbind(obs.num_strata[,1], obs.num_strata[,-1] + matrix(rep(cumsum(dr$n_levels[-dr$n_cov]),
                                                                          nrow(obs.num_strata)), nrow = nrow(obs.num_strata) ,byrow = T ) )

  poc.N = poc.A = poc.D = numeric(sum(dr$n_levels))
  sel.poc = map.strata[sel,]
  poc.N[sel.poc] = poc.N[sel.poc] + 1
  poc.A[sel.poc] = poc.A[sel.poc] + delta[1]

  for (i in 2:n) {
    sel = strata[i]
    sel.poc = map.strata[sel,]

    imb.H <- c(2 * sum(delta[1:(i - 1)]) - (i - 1), 2 * poc.A[sel.poc] - poc.N[sel.poc], D.strata[sel])
    imb.HA <- (omega%*%(imb.H + 1)^2)[1]
    imb.HB <- (omega%*%(imb.H - 1)^2)[1]

    if (imb.HA == imb.HB) {delta[i] <- rbinom(1, 1, .5)
    } else if (imb.HA > imb.HB) {delta[i] <- rbinom(1, 1, 1 - p)
    } else {delta[i] <- rbinom(1, 1, p)
    }
    N.strata[sel] = N.strata[sel] + 1
    A.strata[sel] = A.strata[sel] + delta[i]
    D.strata = 2 * A.strata - N.strata
    poc.N[sel.poc] = poc.N[sel.poc] + 1
    poc.A[sel.poc] = poc.A[sel.poc] + delta[i]

  }

  whit.cov = list()
  for (nc in 1:n_cov) {
    imb.temp = numeric()
    for (j in 1:n_levels[nc]) {
      N.tot = sum(num_strata[, nc] == j)
      A.tot = sum(delta[num_strata[, nc] == j])
      imb.temp[j] = 2*A.tot - N.tot
    }
    whit.cov[[nc]] = imb.temp
    names(whit.cov[[nc]]) = paste(dr$var_levels[[nc]][-1])
  }
  names(whit.cov) = colnames(dr$dm)
  Imb.measures = Imb.m(delta, dr$Fn)
  res = list(
    summary.info = list(
      Design = "Hu and Hu",
      Sample_size = n,
      n_cov = n_cov,
      n_levels = paste(n_levels, collapse = " "),
      var_names = paste(colnames(obs.strata), collapse = " "),
      cov_levels_names = dr$var_levels
    ),
    Assignments = factor(delta, labels = c("B", "A"), levels = 0:1),
    Imbalances.summary = Imb.measures,
    Strata.measures = cbind(obs.strata, N.strata, A.strata, D.strata),
    Imbalances = list(
      Imb.measures = Imb.measures[-3],
      Overall.imb = Imb.measures[[3]],
      Within.strata = cbind(obs.strata, D.strata),
      Within.cov = whit.cov
    ),
    data = data,
    observed.strata = obs.strata
  )

  class(res) = "covadap"
  if(print.results){
    print_covadap(res)
  }

  invisible(res)
}
