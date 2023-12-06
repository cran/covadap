# Big Stick Design ####
BSD.sim <- function(data = NULL,
                    covar = NULL,
                    n = NULL,
                    bound = 3,
                    nrep = 1000,
                    print.results = TRUE) {
  if (is.null(data) & is.null(covar)) {
    stop("Either data or covar must be provided")
  }
  if (!is.null(data) & !is.null(covar)) {
    stop("Either data or covar must be provided")
  }

  if(!is.null(data) & !is.null(n)) {
    stop("Either data or n must be provided")
  }


  if (is.null(data) & is.null(n)) {
    stop("n must be provided")
  }
  if (!is.numeric(bound) | bound < 1) {
    stop("bound must be a postive integer")
  }
  bound = round(bound)

  if (is.null(data)) {
    if (is.numeric(covar)) {
      type = var_levels = list()
      var_names = character()
      for (nc in 1:length(covar)) {
        type[[nc]] = as.factor(1:covar[nc])
        var_levels[[nc]] = c(paste0("Var", nc), levels(type[[nc]]))
        var_names[nc] = var_levels[[nc]][1]
      }
      n_levels = covar
    } else{
      type = covar
      var_levels = list()
      n_levels = numeric()
      for (l in 1:length(type)) {
        var_levels[[l]] = c(names(type)[l], type[[l]])
        n_levels[l] = length(type[[l]])
      }
      var_names = names(type)
    }
    type1 = lapply(type, as.factor)
    obs.strata = as.data.frame(expand.grid(type1, KEEP.OUT.ATTRS = T))
    num.level = matrix(sapply(obs.strata, unclass), nrow = nrow(obs.strata))
    n_cov = ncol(obs.strata)

    if(n_cov==1){
      if(n_levels==1){
        stop("A covariate with at least two levels must be provided")
      }}
    if (sum(n_levels == 1) == length(n_levels)) {
      ##### controlla / cambia #####
      num.strata = matrix(1, 1, 1)
    } else if(n_cov==1){
      num.strata = model.matrix(~ obs.strata[, n_levels > 1])
    } else {
      num.strata = model.matrix( ~ ., data = obs.strata[, n_levels > 1])
    }
    n.strata = nrow(obs.strata)
    n_cov = ncol(obs.strata)

  } else{
    dr = data_preproc(data)
    n = dr$n
    n_cov = dr$n_cov
    n.strata = dr$n.strata
    n_levels = dr$n_levels
    var_levels = dr$var_levels
    obs.strata = dr$obs.strata
    num.level = matrix(sapply(obs.strata, unclass), nrow = nrow(obs.strata))
    num.strata = model.matrix( ~ ., data = obs.strata)
    var_names = colnames(obs.strata)
    prob.s = rep(1 / n.strata, n.strata)
  }

  res = res.data = list()
  within.imb = matrix(NA, nrow = nrep, ncol = sum(n_levels))
  var.l.names = paste(var_levels[[1]][1], var_levels[[1]][-1], sep = "_")
  if(n_cov>1) {
    for (nc in 2:n_cov) {
      var.l.names = c(var.l.names, paste(var_levels[[nc]][1], var_levels[[nc]][-1], sep = "_"))
    }}
  colnames(within.imb) = var.l.names
  Imb.measures = matrix(NA, nrow = nrep, ncol = 3)
  D.sim = A.sim = N.sim = matrix(NA, nrow = nrep, ncol = n.strata)

  for (r in 1:nrep) {
    delta = numeric(n)
    A.strata = N.strata = D.strata = numeric(n.strata)

    if (is.null(data)) {
      strata = sample(1:n.strata, n, T)
    } else{
      strata = sample(1:n.strata, n, T, prob = prob.s)
    }

    Fn = as.matrix(num.strata[strata, ])
    num_strata = num.level[strata, ]

    for (i in 1:n) {
      sel = strata[i]
      Dca = D.strata[sel]
      if (abs(Dca) == bound) {
        delta[i] =  c(1, 0)[Dca == c(-bound, bound)]
      } else {
        delta[i] =  rbinom(1, 1, .5)
      }
      N.strata[sel] = N.strata[sel] + 1
      A.strata[sel] = A.strata[sel] + delta[i]
      D.strata = 2 * A.strata - N.strata
    }

    ###
    start = stop = 0
    for (nc in 1:n_cov) {
      start =  stop + 1
      stop = start + n_levels[nc] - 1
      imb.temp = numeric()

      if( length(n_levels) ==1 ){
        for (j in 1:n_levels[nc]) {
          N.tot = sum(num_strata[nc] == j)
          A.tot = sum(delta[num_strata[nc] == j])
          imb.temp[j] = 2*A.tot - N.tot
        }
        within.imb[r, start:stop] = imb.temp
      } else {
        for (j in 1:n_levels[nc]) {
          N.tot = sum(num_strata[, nc] == j)
          A.tot = sum(delta[num_strata[, nc] == j])
          imb.temp[j] = 2*A.tot - N.tot
        }
        within.imb[r, start:stop] = imb.temp
      }

    }
    ###


    Imb.measures[r, ] = Imb.m.sim(delta, Fn)
    D.sim[r, ] = D.strata
    A.sim[r, ] = A.strata
    N.sim[r, ] = N.strata

    res.data[[r]] = list(data = obs.strata[strata, ],
                         Assignment = factor(delta, labels = c("B", "A"), levels = 0:1))

  } #end MC
  colnames(Imb.measures) <- c("Loss", "Mahal", "overall.imb")
  colnames(D.sim) <- colnames(A.sim) <- colnames(N.sim) <- 1:ncol(D.sim)
  res$summary.info = list(
    Design = "Big Stick Design",
    Sample_size = n,
    n_cov = n_cov,
    n_levels = paste(n_levels, collapse = " "),
    var_names = paste(var_names, collapse = " "),
    n.rep = nrep,
    Maximum_tolerated_imbalance = bound
  )
  res$Imbalances = list(
    Imb.measures = Imb.measures,
    within.imb = within.imb,
    strata.imb = D.sim,
    strata.A = A.sim,
    strata.N = N.sim,
    obs.strata = obs.strata
  )
  res$out = res.data
  class(res) = "covadapsim"
  if (print.results) {
    print_sim(res)
  }
  invisible(res)

}
