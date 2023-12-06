# D_A-optimum biased coin design
DABCD <- function(data,
                  all.cat,
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

  #### all cat
  if(all.cat){
    dr = data_preproc(data)
    n = dr$n
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
    Fc = p.inverse((Fn[1, ])%*%t(Fn[1, ]))
    bt = t(2*delta[1] - 1)%*%Fn[1, ]

    for (i in 2:n) {

      sel = strata[i]
      cn1 = (t(Fn[i,])%*%Fc%*%t(bt))[1]
      P.atk = 0.5 - cn1/(1+cn1^2)
      delta[i] = rbinom(1, 1, P.atk)
      Fc = p.inverse(t(Fn[1:i, ])%*%Fn[1:i, ])
      bt = bt + (2 * delta[i] - 1)*Fn[i, ]
      #bt <- t(2*delta[1:i] - 1)%*%Fn[1:i, ]
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
        Design = "D_A-optimum BCD",
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

  }else{ # mixed covariates

    f.index = sapply(data, is.factor)
    if(sum(f.index) ==1 ){
      Fn = as.matrix(cbind(model.matrix(~ data[,f.index]), data[,!f.index]))
    } else if(sum(f.index) > 1){
      Fn = as.matrix(cbind(model.matrix(~. , data[,f.index]), data[,!f.index]))
    }else{
      Fn = as.matrix(cbind(1, data))
    }

    n = nrow(data)
    delta = numeric(n)
    n_cov = ncol(data)

    delta[1] = rbinom(1,1,.5)
    Fc = p.inverse((Fn[1, ])%*%t(Fn[1, ]))
    bt = t(2*delta[1] - 1)%*%Fn[1, ]

    for (i in 2:n) {

      cn1 = (t(Fn[i,])%*%Fc%*%t(bt))[1]
      P.atk = 0.5 - cn1/(1+cn1^2)
      delta[i] = rbinom(1, 1, P.atk)
      Fc = p.inverse(t(Fn[1:i, ])%*%Fn[1:i, ])
      bt = bt + (2 * delta[i] - 1)*Fn[i, ]

    }

    if(sum(f.index) > 0){
      case = 1

      if( sum(f.index) ==1 ){
        n_cov_cat = sum(f.index)
        lev = unique(data[, f.index])
        nn <- colnames(data)
        var_levels = c(nn[f.index], as.character(lev))
        n_levels = length(lev)

        whit.cov = list()
        imb.temp = numeric()

        for (j in 1:n_levels) {

          N.tot = sum(data[, f.index] == var_levels[-1][j])
          A.tot = sum(delta[data[, f.index] == var_levels[-1][j]])
          imb.temp[j] = 2*A.tot - N.tot
        }

        whit.cov = imb.temp
        names(whit.cov) = paste(var_levels[-1])

      } else {
        n_cov_cat = sum(f.index)
        var_levels = list()
        n_levels = numeric()

        for (nc in 1:n_cov_cat) {

          lev = unique(data[, f.index][,nc])
          var_levels[[nc]] = c(colnames(data[, f.index])[nc], as.character(lev))
          n_levels[nc] = length(lev)

        }

        whit.cov = list()

        for (nc in 1:n_cov_cat) {

          imb.temp = numeric()

          for (j in 1:n_levels[nc]) {

            N.tot = sum(data[, f.index][,nc] == var_levels[[nc]][-1][j])
            A.tot = sum(delta[data[, f.index][,nc] == var_levels[[nc]][-1][j]])
            imb.temp[j] = 2*A.tot - N.tot

          }

          whit.cov[[nc]] = imb.temp
          names(whit.cov[[nc]]) = paste(var_levels[[nc]][-1])

        }
      }

      names(whit.cov) = colnames(data[, f.index])
      Imb.measures = Imb.m(delta, Fn)
      dif.cont = numeric()

      for(vc in which(!f.index == T)){
        dif.cont = append(dif.cont, diff(tapply(data[, vc], list(delta), mean, na.rm = T)))
      }

      names(dif.cont) = colnames(data)[!f.index]

      res = list(
        summary.info = list(
          Design = "D_A-optimum BCD",
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

      Imb.measures = Imb.m(delta, Fn)

      dif.cont = numeric()

      for (vc in which(!f.index == T)) {
        dif.cont = append(dif.cont, diff(tapply(data[, vc], list(delta), mean, na.rm = T)))

      }

      names(dif.cont) = colnames(data)[!f.index]

      res = list(
        summary.info = list(
          Design = "D_A-optimum BCD",
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

