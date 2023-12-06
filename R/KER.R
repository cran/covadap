# Covariate-Adaptive randomization by Ma and Hu
KER <- function(data,
                all.cat,
                p = 0.8,
                print.results = TRUE){

  if(missing(all.cat)){
    stop("Specify whether all covariates are categorical")
  }

  if(missing(data)){

    stop("data must be provided")

  }

  if(!all.cat && !is.data.frame(data)){

    stop("If all.cat = FALSE then data must be provided as a data.frame")

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

    delta[1] = rbinom(1,1,.5)
    sel = strata[1]
    N.strata[sel] = N.strata[sel] + 1
    A.strata[sel] = A.strata[sel] + delta[1]
    D.strata = 2 * A.strata - N.strata
    Fn = dr$Fn

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
        Design = "Kernel-Minimization",
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
    if (print.results) {
      print_covadap(res)
    }

    invisible(res)

  }else{

    f.index = sapply(data, is.factor)

    if(sum(f.index) == ncol(data)){
      stop("for categorical variables only all.cat must be TRUE")
    }
    n = nrow(data)
    if(n<=10){
      stop("A larger n is required")
    }
    delta = numeric(n)
    n_cov = ncol(data)

    if(sum(f.index) > 0){
      n_cov_cat = sum(f.index)
      var_levels = list()
      n_levels = numeric()

      if(sum(f.index)==1){
        Fn1 = as.matrix(model.matrix(~ data[,f.index]))

        lev = unique(data[, f.index])
        nn <- colnames(data)
        var_levels = c(nn[f.index], as.character(lev))
        n_levels = length(lev)

      } else {
        Fn1 = as.matrix(model.matrix(~ ., data[,f.index]))
        for (nc in 1:n_cov_cat) {

          lev = unique(data[, f.index][,nc])
          var_levels[[nc]] = c(colnames(data[, f.index])[nc], as.character(lev))
          n_levels[nc] = length(lev)

        }
      }

      Fn2 = as.matrix(data[,!f.index])
      Fn = as.matrix(cbind(Fn1, Fn2))

      case = 1

      delta[1:10][sample(1:10, 5, F)] = 1

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

      whit.cov = list()

      for (nc in 1:n_cov_cat) {
        imb.temp = numeric()
        if(n_cov_cat==1){

          for (j in 1:n_levels[nc]) {
            N.tot = sum(data[, f.index] == var_levels[-1][j])
            A.tot = sum(delta[data[, f.index] == var_levels[-1][j]])
            imb.temp[j] = 2*A.tot - N.tot
          }
          whit.cov = imb.temp
          names(whit.cov) = var_levels[-1]

        } else {
          for (j in 1:n_levels[nc]) {
            N.tot = sum(data[, f.index][,nc] == var_levels[[nc]][-1][j])
            A.tot = sum(delta[data[, f.index][,nc] == var_levels[[nc]][-1][j]])
            imb.temp[j] = 2*A.tot - N.tot
          }
          whit.cov[[nc]] = imb.temp
          names(whit.cov[[nc]]) = paste(var_levels[[nc]][-1])
        }
      }

      Imb.measures = Imb.m(delta, Fn)
      dif.cont = numeric()

      for(vc in which(!f.index == T)){
        dif.cont = append(dif.cont, diff(tapply(data[, vc], list(delta), mean, na.rm = T)))
      }

      names(dif.cont) = colnames(data)[!f.index]

      res = list(
        summary.info = list(
          Design = "Kernel-Minimization",
          Sample_size = n,
          n_cov = n_cov,
          n_categorical_variables = n_cov_cat,
          n_levels = paste(n_levels, collapse = " "),
          n_quantitative_variables = sum(!f.index),
          var_names = paste(colnames(data), collapse = " "),
          cov_levels_names = var_levels
        ),

        Assignments = factor(delta, labels = c("B", "A"), levels = 0:1),
        Imbalances.summary = Imb.measures,
        Imbalances = list(
          Imb.measures = Imb.measures[-3],
          Overall.imb = Imb.measures[[3]],
          Within.cov = whit.cov
        ),

        diff_mean = dif.cont,
        data = data
      )


    }else{

      case = 2
      Fn = as.matrix(cbind(1, data))
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


      Imb.measures = Imb.m(delta, Fn)

      dif.cont = numeric()

      for (vc in which(!f.index == T)) {
        dif.cont = append(dif.cont, diff(tapply(data[, vc], list(delta), mean, na.rm = T)))

      }

      names(dif.cont) = colnames(data)[!f.index]

      res = list(
        summary.info = list(
          Design = "Kernel-Minimization",
          Sample_size = n,
          n_cov = n_cov,
          n_quantitative_variables = sum(!f.index),
          var_names = paste(colnames(data), collapse = " ")

        ),

        Assignments = factor(delta, labels = c("B", "A"), levels = 0:1),
        Imbalances.summary = Imb.measures,
        Imbalances = list(Imb.measures = Imb.measures[-3],
                          Overall.imb = Imb.measures[[3]]),
        diff_mean = dif.cont,
        data = data

      )

    }
    class(res) = "covadap"
    if (print.results) {
      print_mixed.covadap(res, case)
    }

    invisible(res)

  }

}

