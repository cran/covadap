# Pocock and Simon's minimization method
PocSim.sim <- function(data = NULL,
           covar = NULL,
           n = NULL,
           p = 0.85,
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
    if(!is.numeric(p)){
      stop("p must be a postive number in [0.5, 1]")
    }
    if(p<0.5 | p >1){
      stop("p must be a postive number between in [0.5, 1]")
    }

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
      obs.strata = as.data.frame(expand.grid(type1, KEEP.OUT.ATTRS = T)) # strata with original labels
      num.level = matrix(sapply(obs.strata, unclass), nrow = nrow(obs.strata)) # strata with numeric labels
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

    } else{  #se vengono forniti i dati in input
      dr = data_preproc(data)
      n = dr$n
      n_cov = dr$n_cov
      n.strata = dr$n.strata
      n_levels = dr$n_levels
      var_levels = dr$var_levels
      obs.strata = dr$obs.strata # strata with original labels
      num.level = matrix(sapply(obs.strata, unclass), nrow = nrow(obs.strata)) # strata with numeric labels
      num.strata = model.matrix(~ ., data = obs.strata)
      var_names = colnames(obs.strata)
      prob.s = rep(1 / n.strata, n.strata)

    }
    if(n_cov==1){
      if(n_levels==1){
        stop("A covariate with at least two levels must be provided")
      }}
    # prepare empty boxes for results
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
      poc.N = poc.A = poc.D = numeric(sum(n_levels))

      if (is.null(data)) {
        strata = sample(1:n.strata, n, T)
      } else{
        strata = sample(1:n.strata, n, T, prob = prob.s)
      }

      Fn = as.matrix(num.strata[strata,])
      num_strata = num.level[strata,]

      map.strata = cbind(num.level[,1], num.level[,-1] + matrix(rep(cumsum(n_levels[-n_cov]),
                                                                    nrow(num.level)), nrow = nrow(num.level) ,byrow = T ) )

      delta[1] = rbinom(1,1,.5)
      sel = strata[1]
      N.strata[sel] = N.strata[sel] + 1
      A.strata[sel] = A.strata[sel] + delta[1]

      sel.poc = map.strata[sel,]
      poc.N[sel.poc] = poc.N[sel.poc] + 1
      poc.A[sel.poc] = poc.A[sel.poc] + delta[1]

      for (i in 2:n) {
        sel = strata[i]
        sel.poc = map.strata[sel,]

        t.imb = sum(2 * poc.A[sel.poc] - poc.N[sel.poc])

        if (t.imb == 0) {delta[i] <- rbinom(1, 1, .5)
        } else if (t.imb >= 1) {delta[i] <- rbinom(1, 1, 1 - p)
        } else {delta[i] <- rbinom(1, 1, p)
        }

        N.strata[sel] = N.strata[sel] + 1
        A.strata[sel] = A.strata[sel] + delta[i]

        poc.N[sel.poc] = poc.N[sel.poc] + 1
        poc.A[sel.poc] = poc.A[sel.poc] + delta[i]
      }
      within.imb[r, ] = 2 * poc.A - poc.N
      Imb.measures[r,] = Imb.m.sim(delta, Fn)
      D.sim[r,] = 2 * A.strata - N.strata
      A.sim[r,] = A.strata
      N.sim[r,] = N.strata

      res.data[[r]] = list(data = obs.strata[strata,],
                           Assignment = factor(delta, labels = c("B", "A"), levels = 0:1))

    }
    colnames(Imb.measures) <- c("Loss", "Mahal", "overall.imb")
    colnames(D.sim) <- colnames(A.sim) <- colnames(N.sim) <- 1:ncol(D.sim)

    res$summary.info = list(
      Design = "Pocock and Simon",
      Sample_size = n,
      n_cov = n_cov,
      n_levels = paste(n_levels, collapse = " "),
      var_names = paste(var_names, collapse = " "),
      n.rep = nrep
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

