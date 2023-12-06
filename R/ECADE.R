# Efficient Covariate-Adaptive Design
ECADE <- function(data,
                  all.cat,
                  rho = 0.85,
                  alloc.function = "Efron",
                  print.results = TRUE){

  if(missing(data)){

    stop("data must be provided")

  }

  if(missing(all.cat)){

    stop("Specify whether all covariates are categorical")

  }

  if(!all.cat && !is.data.frame(data)){

    stop("If all.cat = FALSE then data must be provided as a data.frame")

  }

  if(is.null(rho) && alloc.function=="Efron"){

    stop("Provide a value for rho")

  }

  if((!is.numeric(rho) | (rho <= 0.5 | rho >= 1)) && alloc.function=="Efron"){

    stop("rho must be a postive number in (0.5; 1)")

  }

  type.alloc = pmatch(tolower(alloc.function), c("efron", "norm1", "norm2"))

  if(is.na(type.alloc)){
    stop("alloc function must be one of efron, norm1, norm2")
  }

  if(type.alloc == 1){

    alloc.fun =function(imb, rho){
      imb = round(imb, 10)
      if(imb < 0){

        out =  rbinom(1,1, rho)

      }else if(imb > 0){

        out = rbinom(1,1, 1 - rho)

      }else{

        out = rbinom(1,1, .5)

      }

      return(out)
    }

  }else if(type.alloc == 2){

    alloc.fun = function(imb, rho){

      out = rbinom(1,1, prob =  0.1 + 0.8 * pnorm(imb, lower.tail = F) )

      return(out)
    }

  }else if(type.alloc == 3){

    alloc.fun = function(imb, rho){

      out = rbinom(1,1, prob =  0.2 + 0.6 * pnorm(imb, lower.tail = F) )

      return(out)
    }

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
    bt = t(2*delta[1] - 1)%*%Fn[1, ]

    for (i in 2:n) {

      sel = strata[i]
      Fc = p.inverse(t(Fn[1:i, ])%*%Fn[1:i, ])
      imb = (4*Fn[i, ]%*%(i*Fc)%*%t(bt))[1,1]
      delta[i] = alloc.fun(imb, rho)
      bt = bt + (2 * delta[i] - 1)*Fn[i, ]
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
        imb.temp[j] = 2 * A.tot - N.tot
      }

      whit.cov[[nc]] = imb.temp
      names(whit.cov[[nc]]) = paste(dr$var_levels[[nc]][-1])

    }

    names(whit.cov) = colnames(dr$dm)
    Imb.measures = Imb.m(delta, dr$Fn)

    res = list(

      summary.info = list(
        Design = "Efficient Covariate-Adaptive Design",
        Sample_size = n,
        n_cov = n_cov,
        n_levels = paste(n_levels, collapse = " "),
        var_names = paste(colnames(obs.strata), collapse = " "),
        cov_levels_names = dr$var_levels,
        allocation_function = c("Efron", "Norm_0.1", "Norm_0.2")[type.alloc]
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
    if(sum(f.index) ==1 ){
      Fn = as.matrix(cbind(model.matrix(~ data[,f.index]), data[,!f.index]))
    } else if(sum(f.index) > 1){
      Fn = as.matrix(cbind(model.matrix(~. , data[,f.index]), data[,!f.index]))
    }else{
      Fn = as.matrix(cbind(1, data))
    }

    n = nrow(data)
    if(n < 3){
      stop("at least 3 observations are needed")
    }
    delta = numeric(n)
    n_cov = ncol(data)

    delta[1] = rbinom(1,1,.5)
    bt = t(2*delta[1] - 1)%*%Fn[1, ]

    for (i in 2:n) {

      Fc = p.inverse(t(Fn[1:i, ])%*%Fn[1:i, ])
      imb = (4*Fn[i, ]%*%(i*Fc)%*%t(bt))[1,1]
      delta[i] = alloc.fun(imb, rho)
      bt = bt + (2 * delta[i] - 1)*Fn[i, ]

    }

    if(sum(f.index) == length(f.index)){

      stop("With only categorical covariates all.cat must be set to TRUE")

    }else if(sum(f.index) > 0){

      n_cov_cat = sum(f.index)

      if(n_cov_cat==1){
        lev = unique(data[, f.index])
        nn = colnames(data)
        var_levels = c(nn[f.index], as.character(lev))
        n_levels = length(lev)

      } else {
        n_levels = numeric()
        var_levels = list()
        for (nc in 1:n_cov_cat) {

          lev = unique(data[, f.index][,nc])
          var_levels[[nc]] = c(colnames(data[, f.index])[nc], as.character(lev))
          n_levels[nc] = length(lev)

        }
      }
      whit.cov = list()

      for (nc in 1:n_cov_cat) {

        imb.temp = numeric()
        if(n_cov_cat==1){
          for (j in 1:n_levels) {
            N.tot = sum(data[, f.index] == var_levels[-1][j])
            A.tot = sum(delta[data[, f.index] == var_levels[-1][j]])
            imb.temp[j] = 2*A.tot - N.tot
          }

          whit.cov = imb.temp
          names(whit.cov) = var_levels[-1]

        } else{
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
          Design = "Efficient Covariate-Adaptive Design",
          Sample_size = n,
          n_cov = n_cov,
          n_categorical_variables = n_cov_cat,
          n_levels = paste(n_levels, collapse = " "),
          n_quantitative_variables = sum(!f.index),
          var_names = paste(colnames(data), collapse = " "),
          cov_levels_names = var_levels,
          allocation_function = c("Efron", "Norm_0.1", "Norm_0.2")[type.alloc]
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
      case = 1

    }else{

      case = 2

      Imb.measures = Imb.m(delta, Fn)

      dif.cont = numeric()

      for (vc in which(!f.index == T)) {
        dif.cont = append(dif.cont, diff(tapply(data[, vc], list(delta), mean, na.rm = T)))

      }

      names(dif.cont) = colnames(data)[!f.index]

      res = list(
        summary.info = list(
          Design = "Efficient Covariate-Adaptive Design",
          Sample_size = n,
          n_cov = n_cov,
          n_quantitative_variables = sum(!f.index),
          var_names = paste(colnames(data), collapse = " "),
          allocation_function = c("Efron", "Norm_0.1", "Norm_0.2")[type.alloc]

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

