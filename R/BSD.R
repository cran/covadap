# Big Stick Design ####
BSD <- function(data, bound = 3, print.results = TRUE){
  if(missing(data)){
    stop("data must be provided")
  }
  if(!is.numeric(bound) | bound < 1){
    stop("bound must be a postive integer")
  }
  bound = round(bound)

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

  for (i in 1:n) {
    sel = strata[i]
    Dca = D.strata[sel]

    if (abs(Dca) == bound) {delta[i] =  c(1,0)[Dca == c(-bound,bound)]
    } else {delta[i] =  rbinom(1, 1, .5)}
    N.strata[sel] = N.strata[sel] + 1
    A.strata[sel] = A.strata[sel] + delta[i]
    D.strata = 2 * A.strata - N.strata

  }
  num_strata = dr$num_strata
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
      Design = "Big Stick Design",
      Sample_size = n,
      n_cov = n_cov,
      n_levels = paste(n_levels, collapse = " "),
      var_names = paste(colnames(obs.strata), collapse = " "),
      cov_levels_names = dr$var_levels,
      Maximum_tolerated_imbalance = bound
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
  if (print.results) {
    print_covadap(res)
  }
  invisible(res)
}
