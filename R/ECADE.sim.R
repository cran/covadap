#Efficient Covariate-Adaptive Design
ECADE.sim <- function(data = NULL,
                      covar = NULL,
                      n = NULL,
                      all.cat,
                      nrep = 1000,
                      rho = 0.85,
                      alloc.function = "Efron",
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
  if(is.null(rho) && alloc.function=="Efron"){
    stop("provide a value for rho")
  }

  if((!is.numeric(rho) | (rho <= 0.5 | rho >= 1)) && alloc.function=="Efron"){
    stop("rho must be a postive number in (0.5; 1)")
  }

  type.alloc = pmatch(tolower(alloc.function), c("efron", "norm1", "norm2"))
  if(is.na(type.alloc)){
    stop("alloc function must be one of efron, norm1, norm2")
  }

  if(type.alloc == 1) {
    alloc.fun = function(imb, rho) {
      imb = round(imb, 10)
      if (imb < 0) {
        out =  rbinom(1, 1, rho)

      } else if (imb > 0) {
        out = rbinom(1, 1, 1 - rho)

      } else{
        out = rbinom(1, 1, .5)

      }

      return(out)
    }

  } else if (type.alloc == 2) {
    alloc.fun = function(imb, rho) {
      out = rbinom(1, 1, prob =  0.1 + 0.8 * pnorm(imb, lower.tail = F))

      return(out)
    }

  } else if (type.alloc == 3) {
    alloc.fun = function(imb, rho) {
      out = rbinom(1, 1, prob =  0.2 + 0.6 * pnorm(imb, lower.tail = F))

      return(out)
    }

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
      obs.strata = dr$obs.strata # strata with original labels
      num.level = matrix(sapply(obs.strata, unclass), nrow = nrow(obs.strata)) # strata with numeric labels
      num.strata = model.matrix(~ ., data = obs.strata)
      var_names = colnames(obs.strata)
      prob.s = rep(1 / n.strata, n.strata)
    }

    res = res.data = list()
    within.imb = matrix(NA, nrow = nrep, ncol = sum(n_levels))
    var.l.names = paste(var_levels[[1]][1], var_levels[[1]][-1], sep = "_")
    if(n_cov>1){
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
      num_strata = as.data.frame(num.level[strata,])

      delta[1] = rbinom(1,1,.5)
      sel = strata[1]
      N.strata[sel] = N.strata[sel] + 1
      A.strata[sel] = A.strata[sel] + delta[1]
      D.strata = 2 * A.strata - N.strata
      bt = t(2*delta[1] - 1)%*%Fn[1, ]
      c.temp = (Fn[1, ])%*%t(Fn[1, ])
      for (i in 2:n) {

        sel = strata[i]
        c.temp = c.temp + Fn[i, ]%*%t(Fn[i, ])
        Fc = p.inverse(c.temp)
        imb = (4*Fn[i, ]%*%(i*Fc)%*%t(bt))[1,1]
        delta[i] = alloc.fun(imb, rho)
        bt = bt + (2 * delta[i] - 1)*Fn[i, ]
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
      Design = "Efficient Covariate-Adaptive Design",
      Sample_size = n,
      n_cov = n_cov,
      n_levels = paste(n_levels, collapse = " "),
      var_names = paste(var_names, collapse = " "),
      n.rep = nrep,
      allocation_function = c("Efron", "Norm_0.1", "Norm_0.2")[type.alloc]
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
          num.strata= model.matrix( ~ obs.strata[, n_levels > 1])
        } else {num.strata = model.matrix( ~ ., data = obs.strata[, n_levels > 1])
        }

        n.strata = nrow(obs.strata)
        n_cov = ncol(obs.strata)

        var_names = names(c(covar[[1]], covar[[2]]))
        type.quant = covar[[which(type.cov == 2)]]
        quant.par = matrix(unlist(type.quant), ncol = 2, byrow = TRUE)
        n_cov_tot = length(var_names)
        n_cov_quant = length(type.quant)
        n_cov_cat = n_cov_tot-n_cov_quant

        res = res.data = list()
        within.imb = matrix(NA, nrow = nrep, ncol = sum(n_levels))
        var.l.names = paste(var_levels[[1]][1], var_levels[[1]][-1], sep = "_")

        if(length(n_levels) ==1){
          var.l.names = var.l.names
        } else {
          for (nc in 2:n_cov) {
            var.l.names = c(var.l.names, paste(var_levels[[nc]][1], var_levels[[nc]][-1], sep = "_"))
          }}

        colnames(within.imb) = var.l.names
        Imb.measures = matrix(NA, nrow = nrep, ncol = 3)
        dif.cont = matrix(NA, nrow = nrep, ncol = n_cov_quant)
        colnames(dif.cont) = names(type.quant)

        for (r in 1:nrep) {

          delta = numeric(n)
          strata = sample(1:n.strata, n, T)
          Fn = cbind(as.matrix(num.strata[strata,]), matrix(rnorm(n_cov_quant*n, mean = quant.par[,1], sd = quant.par[,2]), nrow = n, byrow = T ))
          num_strata = num.level[strata,]

          delta[1] = rbinom(1,1,.5)
          bt = t(2*delta[1] - 1)%*%Fn[1, ]
          c.temp = (Fn[1, ])%*%t(Fn[1, ])

          for (i in 2:n) {
            sel = strata[i]
            c.temp = c.temp + Fn[i, ]%*%t(Fn[i, ])
            Fc = p.inverse(c.temp)
            imb = (4*Fn[i, ]%*%(i*Fc)%*%t(bt))[1,1]
            delta[i] = alloc.fun(imb, rho)
            bt = bt + (2 * delta[i] - 1)*Fn[i, ]
          }

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
          Design = "Efficient Covariate-Adaptive Design",
          Sample_size = n,
          n_cov = n_cov_tot,
          var_names = paste(var_names, collapse = " "),
          n_quantitative_variables = n_cov_quant,
          n_categorical_variables = n_cov,
          n_levels = paste(n_levels, collapse = " "),
          n.rep = nrep,
          allocation_function = c("Efron", "Norm_0.1", "Norm_0.2")[type.alloc]
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
        quant.par = matrix(unlist(covar), ncol = 2, byrow = TRUE)
        type.quant = covar[[which(type.cov == 2)]]

        n_cov_tot = length(var_names)
        n_cov_quant = length(type.quant)
        n_cov_tot = length(var_names)

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
          delta[1] = rbinom(1, 1, .5)
          bt = t(2 * delta[1] - 1) %*% Fn[1,]
          c.temp = (Fn[1, ])%*%t(Fn[1, ])

          for (i in 2:n) {
            c.temp = c.temp + Fn[i, ]%*%t(Fn[i, ])
            Fc = p.inverse(c.temp)
            imb = (4*Fn[i, ]%*%(i*Fc)%*%t(bt))[1,1]
            delta[i] = alloc.fun(imb, rho)
            bt = bt + (2 * delta[i] - 1)*Fn[i, ]
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
          Design = "Efficient Covariate-Adaptive Design",
          Sample_size = n,
          n_cov = n_cov_tot,
          var_names = paste(var_names, collapse = " "),
          n.rep = nrep,
          allocation_function = c("Efron", "Norm_0.1", "Norm_0.2")[type.alloc]
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

      if(sum(f.index) > 0){

        case = 1
        n_cov_cat = sum(f.index)
        n_cov_quant = ncol(data) - n_cov_cat
        var_levels = list()
        n_levels = numeric()

        if(n_cov_cat==1){
          Fn0 = as.matrix(cbind(model.matrix(~ data[,f.index]), data[,!f.index]))
          lev = unique(data[, f.index])
          var_levels = c(colnames(data[, f.index]), as.character(lev))
          n_levels = length(lev)
          within.imb = matrix(NA, nrow = nrep, ncol = sum(n_levels))
          colnames(within.imb) = var_levels
        } else {

          Fn0 = as.matrix(cbind(model.matrix(~ ., data[,f.index]), data[,!f.index]))

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

      }else{
        case = 2
        data = as.data.frame(data)
        Fn0 = as.matrix(cbind(1, data))
        n_cov_quant = ncol(data)
        Imb.measures = matrix(NA, nrow = nrep, ncol = 3)
        dif.cont = matrix(NA, nrow = nrep, ncol = n_cov_quant)
        colnames(dif.cont) = colnames(data)

      }

      n = nrow(data)
      var_names = colnames(data)

      for (r in 1:nrep) {

        delta = numeric(n)
        ind.temp = sample(1:n, n, T)
        Fn = Fn0[ind.temp,]
        data.r = as.data.frame(data[ind.temp,])
        delta[1] = rbinom(1,1,.5)
        bt = t(2*delta[1] - 1)%*%Fn[1, ]
        c.temp = (Fn[1, ])%*%t(Fn[1, ])

        for (i in 2:n) {
          c.temp = c.temp + Fn[i, ]%*%t(Fn[i, ])
          Fc = p.inverse(c.temp)
          imb = (4*Fn[i, ]%*%(i*Fc)%*%t(bt))[1,1]
          delta[i] = alloc.fun(imb, rho)
          bt = bt + (2 * delta[i] - 1)*Fn[i, ]
        }

        if(case == 1){

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
        }else{

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

        }


      }
      colnames(Imb.measures) <- c("Loss", "Mahal", "overall.imb")

      if(case == 1){
        res$summary.info = list(
          Design = "Efficient Covariate-Adaptive Design",
          Sample_size = n,
          n_cov = ncol(data),
          var_names = paste(var_names, collapse = " "),
          n_quantitative_variables = n_cov_quant,
          n_categorical_variables = n_cov_cat,
          n_levels = paste(n_levels, collapse = " "),
          n.rep = nrep,
          allocation_function = c("Efron", "Norm_0.1", "Norm_0.2")[type.alloc]
        )
        res$Imbalances = list(
          Imb.measures = Imb.measures,
          within.imb = within.imb,
          diff_mean = dif.cont
        )
        res$out = res.data

      }else{
        res$summary.info = list(
          Design = "Efficient Covariate-Adaptive Design",
          Sample_size = n,
          n_cov = ncol(data),
          var_names = paste(var_names, collapse = " "),
          n.rep = nrep,
          allocation_function = c("Efron", "Norm_0.1", "Norm_0.2")[type.alloc]
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
