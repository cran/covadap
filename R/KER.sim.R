# Covariate-Adaptive randomization by Ma and Hu
KER.sim <- function(data = NULL,
                    covar = NULL,
                    n = NULL,
                    all.cat,
                    nrep = 1000,
                    p = 0.8,
                    print.results = TRUE) {
  if(missing(all.cat)){
    stop("Specify whether all covariates are categorical")
  }
  if (all.cat & is.null(data) & is.null(covar)) {
    stop("Either data or covar must be provided")
  }
  if (all.cat & !is.null(data) & !is.null(covar)) {
    stop("Either data or covar must be provided")
  }
  if (is.null(data) & is.null(n)) {
    stop("n must be provided")
  }
  if (!all.cat & !is.list(covar) & is.null(data)) {
    stop("With mixed covariate profiles covar must be a list")
  }
  if(!is.numeric(p) | (p < 0.5 | p > 1)){
    stop("p must be a postive number in [0.5; 1]")
  }

  myKDE <- function(tA, tB, g1, g2, s1){
    h1 = s1 / length(g1)^.2
    h2 = s1 / length(g2)^.2
    kernelValues1 = dnorm((tA - g1) / h1) / h1
    kernelValues2 = dnorm((tB - g2) / h2) / h2
    return(sum(kernelValues1) - sum(kernelValues2))
  }

  alloc.fun = function(imb, p) {
    imb = round(imb, 10)
    if (imb < 0) {
      out =  rbinom(1, 1, p)

    } else if (imb > 0) {
      out = rbinom(1, 1, 1 - p)

    } else{
      out = rbinom(1, 1, .5)

    }

    return(out)
  }




  if(all.cat){

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

      if (sum(n_levels == 1) == length(n_levels)) {
        num.strata = matrix(1, 1, 1)
      } else{
        num.strata = model.matrix( ~ ., data = as.data.frame(obs.strata[, n_levels > 1]))
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
      num.level = matrix(sapply(obs.strata, unclass), nrow = nrow(obs.strata)) # strata with numeric labels
      num.strata = model.matrix(~ ., data = obs.strata)
      var_names = colnames(obs.strata)
      prob.s = rep(1 / n.strata, n.strata)
    }

    res = res.data = list()
    within.imb = matrix(NA, nrow = nrep, ncol = sum(n_levels))
    var.l.names = paste(var_levels[[1]][1], var_levels[[1]][-1], sep = "_")
    if(n_cov > 1){
    for (nc in 2:n_cov) {
      var.l.names = c(var.l.names, paste(var_levels[[nc]][1], var_levels[[nc]][-1], sep = "_"))
    }
    }
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

      Fn = as.matrix(num.strata[strata,])
      num_strata = as.matrix(num.level[strata,])

      delta[1] = rbinom(1,1,.5)
      sel = strata[1]
      N.strata[sel] = N.strata[sel] + 1
      A.strata[sel] = A.strata[sel] + delta[1]
      D.strata = 2 * A.strata - N.strata


      for (i in 2:n) {

        sel = strata[i]
        imb = numeric()

        for(uu in 2:ncol(Fn)){
          seld = Fn[1:(i-1),uu] == Fn[i,uu]
          imb[uu-1] = (sum(delta[1:(i-1)][seld] == 1) - sum(delta[1:(i-1)][seld] == 0) ) / (i - 1)
        }

        delta[i] = alloc.fun(mean(imb), p)
        N.strata[sel] = N.strata[sel] + 1
        A.strata[sel] = A.strata[sel] + delta[i]
        D.strata = 2 * A.strata - N.strata

      }

      start = stop = 0
      for (nc in 1:n_cov) {
        start =  stop + 1
        stop = start + n_levels[nc] - 1
        imb.temp = numeric()
        for (j in 1:n_levels[nc]) {
          N.tot = sum(num_strata[, nc] == j)
          A.tot = sum(delta[num_strata[, nc] == j])
          imb.temp[j] = 2*A.tot - N.tot
        }
        within.imb[r, start:stop] = imb.temp
      }

      Imb.measures[r,] = Imb.m.sim(delta, Fn)
      D.sim[r,] = D.strata
      A.sim[r,] = A.strata
      N.sim[r,] = N.strata

      res.data[[r]] = list(data = obs.strata[strata,],
                           Assignment = factor(delta, labels = c("B", "A"), levels = 0:1))

    }# end MC
    colnames(Imb.measures) <- c("Loss", "Mahal", "overall.imb")

    res$summary.info = list(
      Design =  "Kernel-Minimization",
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

  }else{

    if (is.null(data)) {

      type.cov = pmatch(tolower(names(covar)), c("categorical", "quantitative"))

      if(length(type.cov) == 2){
        type = covar[[which(type.cov == 1)]]
        var_levels = list()
        n_levels = numeric()
        for (l in 1:length(type)) {
          var_levels[[l]] = c(names(type)[l], type[[l]])
          n_levels[l] = length(type[[l]])
        }
        type1 = lapply(type, as.factor)
        obs.strata = as.data.frame(expand.grid(type1, KEEP.OUT.ATTRS = T)) # strata with original labels
        num.level = matrix(sapply(obs.strata, unclass), nrow = nrow(obs.strata)) # strata with numeric labels

        if (sum(n_levels == 1) == length(n_levels)) {
          num.strata = matrix(1, 1, 1)
        } else if(length(n_levels) ==1   ){
          num.strata= model.matrix( ~ data, data = obs.strata[, n_levels > 1])
          var.l.names =  paste(var_levels[[1]][1], var_levels[[1]][-1], sep = "_")
        } else {
          num.strata = model.matrix( ~ ., data = obs.strata[, n_levels > 1])
          var.l.names = paste(var_levels[[1]][1], var_levels[[1]][-1], sep = "_")
          for (nc in 2:n_cov) {
            var.l.names = c(var.l.names, paste(var_levels[[nc]][1], var_levels[[nc]][-1], sep = "_"))
          }}

        n.strata = nrow(obs.strata)
        n_cov = ncol(obs.strata)

        var_names = names(c(covar[[1]], covar[[2]]))
        type.quant = covar[[which(type.cov == 2)]]
        quant.par = matrix(unlist(type.quant), nrow = 2)
        n_cov_tot = length(var_names)
        n_cov_quant = length(type.quant)
        n_cov_cat = n_cov_tot - n_cov_quant
        res = res.data = list()
        within.imb = matrix(NA, nrow = nrep, ncol = sum(n_levels))


        colnames(within.imb) = var.l.names
        Imb.measures = matrix(NA, nrow = nrep, ncol = 3)
        dif.cont = matrix(NA, nrow = nrep, ncol = n_cov_quant)
        colnames(dif.cont) = names(type.quant)

        for (r in 1:nrep) {

          delta = numeric(n)
          strata = sample(1:n.strata, n, T)
          Fn1 = as.matrix(num.strata[strata,])
          Fn2 = matrix(rnorm(n_cov_quant*n, mean = quant.par[,1], sd = quant.par[,2]), nrow = n, byrow = T )
          Fn = cbind(Fn1, Fn2)
          num_strata = num.level[strata,]

          delta[1:10][sample(1:10, 5, F)] = 1


          if(n==10){
            stop("n must be greater than 10")
          }

          for (i in 11:n) {

            imb = numeric()
            for(uu in 2:ncol(Fn1)){
              seld = Fn1[1:(i-1),uu] == Fn1[i,uu]
              imb[uu-1] = (sum(delta[1:(i-1)][seld] == 1) - sum(delta[1:(i-1)][seld] == 0) ) / (i - 1)
            }

            sel = delta[1:(i-1)] == 1
            imb2 = numeric()
            for(u in 1:ncol(Fn2)){

              s1 = sd(Fn2[1:(i-1),u])
              tA = tB = (Fn2[i,u] ) / s1
              g1 = (Fn2[1:(i-1),u][sel]) / s1
              g2 = (Fn2[1:(i-1),u][!sel]) / s1
              imb2[u] = myKDE(tA, tB, g1, g2, s1) / (i - 1)

            }

            imb = c(imb, imb2)
            delta[i] = alloc.fun(mean(imb), p)

          }

          start = stop = 0
          start = stop = 0
          for (nc in 1:n_cov) {
            start =  stop + 1
            stop = start + n_levels[nc] - 1
            imb.temp = numeric()

            if(n_cov_cat==1){
              for (j in 1:n_levels[nc]) {
                N.tot = sum(num_strata[nc] == j)
                A.tot = sum(delta[num_strata[nc] == j])
                imb.temp[j] = 2*A.tot - N.tot
              }
            } else {
              for (j in 1:n_levels[nc]) {
                N.tot = sum(num_strata[, nc] == j)
                A.tot = sum(delta[num_strata[, nc] == j])
                imb.temp[j] = 2*A.tot - N.tot
              }
            }
            within.imb[r, start:stop] = imb.temp
          }

          temp.data = as.data.frame(cbind(obs.strata[strata, ],
                                          Fn[, (ncol(Fn) - n_cov_quant + 1):ncol(Fn)]))
          rownames(temp.data) = NULL
          colnames(temp.data) = var_names
          Imb.measures[r, ] = Imb.m.sim(delta, Fn)
          res.data[[r]] = list(data = temp.data,
                               Assignment = factor(delta,
                                                   labels = c("B", "A"),
                                                   levels = 0:1))

          dif.temp = numeric()
          for (vc in ((ncol(Fn) - n_cov_quant + 1):ncol(Fn))) {
            dif.temp = append(dif.temp,
                              diff(tapply(Fn[, vc],
                                          list(delta),
                                          mean, na.rm = T)))
          }
          dif.cont[r,] = dif.temp

        }
        colnames(Imb.measures) <- c("Loss", "Mahal", "overall.imb")

        res$summary.info = list(
          Design = "Kernel-Minimization",
          Sample_size = n,
          n_cov = n_cov_tot,
          var_names = paste(var_names, collapse = " "),
          n_quantitative_variables = n_cov_quant,
          n_categorical_variables = n_cov,
          n_levels = paste(n_levels, collapse = " "),
          n.rep = nrep
        )
        res$Imbalances = list(
          Imb.measures = Imb.measures,
          within.imb = within.imb,
          diff_mean = dif.cont
        )
        res$out = res.data
        case = 1

      }else if(length(type.cov) == 1 & type.cov == 2){

        var_names = names(covar[[1]])
        quant.par = matrix(unlist(covar), nrow = 2)
        n_cov_tot = length(var_names)
        type.quant = covar[[which(type.cov == 2)]]
        n_cov_quant = length(type.quant)

        res = res.data = list()
        Imb.measures = matrix(NA, nrow = nrep, ncol = 3)
        dif.cont = matrix(NA, nrow = nrep, ncol = n_cov_quant)
        colnames(dif.cont) = var_names

        for (r in 1:nrep) {
          delta = numeric(n)
          Fn = cbind(1, matrix(
            rnorm(n_cov_quant * n,
                  mean = quant.par[, 1],
                  sd = quant.par[, 2]),
            nrow = n,
            byrow = T
          ))

          delta[1:10][sample(1:10, 5, F)] = 1

          for (i in 11:n) {

            sel = delta[1:(i-1)] == 1
            imb = numeric()
            for(u in 2:ncol(Fn)){

              s1 = sd(Fn[1:(i-1),u])
              tA = (Fn[i,u] ) / s1
              tB = (Fn[i,u] ) / s1
              g1 = (Fn[1:(i-1),u][sel]) / s1
              g2 = (Fn[1:(i-1),u][!sel]) / s1
              imb[u-1] = myKDE(tA, tB, g1, g2, s1) / (i - 1)

            }

            delta[i] = alloc.fun(mean(imb), p)

          }

          temp.data = as.data.frame(Fn)
          rownames(temp.data) = NULL
          colnames(temp.data) = var_names
          Imb.measures[r, ] = Imb.m.sim(delta, Fn) # names FN
          res.data[[r]] = list(data = temp.data,
                               Assignment = factor(delta,
                                                   labels = c("B", "A"),
                                                   levels = 0:1))

          dif.temp = numeric()
          for (vc in 2:ncol(Fn)) {
            dif.temp = append(dif.temp,
                              diff(tapply(Fn[, vc],
                                          list(delta),
                                          mean, na.rm = T)))
          }
          dif.cont[r, ] = dif.temp

        }
        colnames(Imb.measures) <- c("Loss", "Mahal", "overall.imb")

        res$summary.info = list(
          Design = "Kernel-Minimization",
          Sample_size = n,
          n_cov = n_cov_tot,
          var_names = paste(var_names, collapse = " "),
          n.rep = nrep
        )
        res$Imbalances = list(
          Imb.measures = Imb.measures,
          diff_mean = dif.cont
        )
        res$out = res.data
        case = 2

      }else{
        stop("With only categorical covariates all.cat must be set to TRUE")
      }
      class(res) = "covadapsim"
      if (print.results) {
        print_mixed.sim(res, case)
      }
      invisible(res)

    } else{
      res = res.data = list()
      f.index = sapply(data, is.factor)
      n = nrow(data)
      var_names = colnames(data)

      if(sum(f.index) > 0){
        case = 1
        n_cov_cat = sum(f.index)
        n_cov_quant = ncol(data) - n_cov_cat
        var_levels = list()
        n_levels = numeric()

        if(n_cov_cat==1){
          Fn0 = as.matrix(cbind(model.matrix(~ data[,f.index]), data[,!f.index]))
          Fn1 = as.matrix(model.matrix(~  data[,f.index]))
          Fn2 = as.matrix(data[,!f.index])

          lev = unique(data[, f.index])
          var_levels = c(colnames(data[, f.index]), as.character(lev))
          n_levels = length(lev)
          within.imb = matrix(NA, nrow = nrep, ncol = sum(n_levels))
          colnames(within.imb) = var_levels
        } else {
          Fn0 = as.matrix(cbind(model.matrix(~ ., data[,f.index]), data[,!f.index]))
          Fn1 = as.matrix(model.matrix(~ ., data[,f.index]))
          Fn2 = as.matrix(data[,!f.index])

          for (nc in 1:n_cov_cat) {
            lev = unique(data[, f.index][,nc])
            var_levels[[nc]] = c(colnames(data[, f.index])[nc], as.character(lev))
            n_levels[nc] = length(lev)
          }

          var.l.names = paste(var_levels[[1]][1], var_levels[[1]][-1], sep = "_")
          for (nc in 2:n_cov_cat) {
            var.l.names = c(var.l.names, paste(var_levels[[nc]][1], var_levels[[nc]][-1], sep = "_"))
          }
          within.imb = matrix(NA, nrow = nrep, ncol = sum(n_levels))
          colnames(within.imb) = var.l.names
        }

        Imb.measures = matrix(NA, nrow = nrep, ncol = 3)
        dif.cont = matrix(NA, nrow = nrep, ncol = n_cov_quant)
        colnames(dif.cont) = colnames(data[,!f.index])

        for (r in 1:nrep) {

          delta = numeric(n)
          ind.temp = sample(1:n, n, T)
          Fn = Fn0[ind.temp,]
          data.r = data[ind.temp,]

          delta[1:10][sample(1:10, 5, F)] = 1

          for (i in 11:n) {
            imb = numeric()
            for(uu in 2:ncol(Fn1)){
              seld = Fn[1:(i-1),uu] == Fn[i,uu]
              imb[uu-1] = (sum(delta[1:(i-1)][seld] == 1) - sum(delta[1:(i-1)][seld] == 0) ) / (i - 1)
            }

            sel = delta[1:(i-1)] == 1
            for(u in ncol(Fn1)+1:ncol(Fn2)){

              s1 = sd(Fn[1:(i-1),u])
              tA = tB = (Fn[i,u] ) / s1
              g1 = (Fn[1:(i-1),u][sel]) / s1
              g2 = (Fn[1:(i-1),u][!sel]) / s1
              imb[u-1] = myKDE(tA, tB, g1, g2, s1) / (i - 1)

            }

            delta[i] = alloc.fun(mean(imb), p)

          }

          start = stop = 0

          for (nc in 1:n_cov_cat) {
            start =  stop + 1
            stop = start + n_levels[nc] - 1
            imb.temp = numeric()

            if(n_cov_cat==1){

              for (j in 1:n_levels) {
                N.tot = sum(data.r[, f.index] == var_levels[j])
                A.tot = sum(delta[data.r[, f.index] == var_levels[j]])
                imb.temp[j] = 2*A.tot - N.tot
              }

              within.imb[r, ] = imb.temp
              #names(within.imb) = paste(var_levels[-1])
            } else{
              for (j in 1:n_levels[nc]) {
                N.tot = sum(data.r[, f.index][,nc] == var_levels[[nc]][-1][j])
                A.tot = sum(delta[data.r[, f.index][,nc] == var_levels[[nc]][-1][j]])
                imb.temp[j] = 2*A.tot - N.tot
              }
              within.imb[r, start:stop] = imb.temp
            }
          }
          ###

          rownames(data.r) = NULL
          colnames(data.r) = var_names
          Imb.measures[r, ] = Imb.m.sim(delta, Fn) # names FN
          res.data[[r]] = list(data = data.r,
                               Assignment = factor(delta,
                                                   labels = c("B", "A"),
                                                   levels = 0:1))

          dif.temp = numeric()
          for (vc in (1:ncol(data))[!f.index]) {
            dif.temp = append(dif.temp,
                              diff(tapply(data.r[, vc],
                                          list(delta),
                                          mean, na.rm = T)))
          }
          dif.cont[r,] = dif.temp

        }
        colnames(Imb.measures) <- c("Loss", "Mahal", "overall.imb")

      }else{
        case = 2
        data = as.data.frame(data)
        n = nrow(data)
        Fn0 = as.matrix(cbind(1, data))
        n_cov_quant = ncol(data)
        Imb.measures = matrix(NA, nrow = nrep, ncol = 3)
        dif.cont = matrix(NA, nrow = nrep, ncol = n_cov_quant)
        colnames(dif.cont) = colnames(data)
        var_names = colnames(data)

        for (r in 1:nrep) {

          delta = numeric(n)
          ind.temp = sample(1:n, n, T)
          Fn = Fn0[ind.temp,]
          data.r = as.data.frame(data[ind.temp,])

          delta[1:10][sample(1:10, 5, F)] = 1

          for (i in 11:n) {
            imb = numeric()

            sel = delta[1:(i-1)] == 1
            for(u in 2:ncol(Fn)){

              s1 = sd(Fn[1:(i-1),u])
              tA = tB = (Fn[i,u] ) / s1
              g1 = (Fn[1:(i-1),u][sel]) / s1
              g2 = (Fn[1:(i-1),u][!sel]) / s1
              imb[u-1] = myKDE(tA, tB, g1, g2, s1) / (i - 1)

            }

            delta[i] = alloc.fun(mean(imb), p)

          }

          rownames(data.r) = NULL
          colnames(data.r) = var_names
          Imb.measures[r, ] = Imb.m.sim(delta, Fn) # names FN
          res.data[[r]] = list(data = data.r,
                               Assignment = factor(delta,
                                                   labels = c("B", "A"),
                                                   levels = 0:1))

          dif.temp = numeric()
          for (vc in 1:ncol(data)) {
            dif.temp = append(dif.temp,
                              diff(tapply(data.r[, vc],
                                          list(delta),
                                          mean, na.rm = T)))
          }
          dif.cont[r,] = dif.temp

        }# end MC
        colnames(Imb.measures) <- c("Loss", "Mahal", "overall.imb")

      }



      if(case == 1){
        res$summary.info = list(
          Design = "Kernel-Minimization",
          Sample_size = n,
          n_cov = ncol(data),
          var_names = paste(var_names, collapse = " "),
          n_quantitative_variables = n_cov_quant,
          n_categorical_variables = n_cov_cat,
          n_levels = paste(n_levels, collapse = " "),
          n.rep = nrep
        )
        res$Imbalances = list(
          Imb.measures = Imb.measures,
          within.imb = within.imb,
          diff_mean = dif.cont
        )
        res$out = res.data

      }else{
        res$summary.info = list(
          Design = "Kernel-Minimization",
          Sample_size = n,
          n_cov = ncol(data),
          var_names = paste(var_names, collapse = " "),
          n.rep = nrep
        )
        res$Imbalances = list(
          Imb.measures = Imb.measures,
          diff_mean = dif.cont
        )
        res$out = res.data

      }
      class(res) = "covadapsim"
      if (print.results) {
        print_mixed.sim(res, case)
      }

      invisible(res)

    }
  }

}
