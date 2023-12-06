######################### COVADAP ##################
# Internal functions
# data processing
data_preproc  <-  function(x) {
  if (length(x[is.na(x)]) > 0) {
    x[is.na(x)] = "NA"
  }
  x = as.data.frame(x)
  xp = data.frame(lapply(x, as.factor))
  sx = unique(xp)
  index_f = sapply(sx, function(x)
    length(levels(x)) > 1)
  row.names(sx) = 1:nrow(sx)
  matched = match(interaction(xp), interaction(sx))
  if (sum(index_f) == 0) {
    sm = matrix(1, nrow(sx), 1)
  } else{
    sm = model.matrix(~ ., data = as.data.frame(xp[, index_f]))
  }
  n = nrow(xp)
  n_cov = ncol(xp)
  var_levels = list()
  n_levels = numeric()
  data_1 = matrix(NA, nrow = n, ncol = n_cov)
  for (p in 1:n_cov) {
    lev = unique(xp[, p])
    var_levels[[p]] = c(colnames(xp)[p], as.character(lev))
    n_levels[p] = length(lev)
    data_1[, p] = match(xp[, p], lev)
  }

  return(
    list(
      dm = xp,
      Fn = sm,
      n = n,
      n_cov = n_cov,
      n.strata = nrow(sx),
      strata = matched,
      obs.strata = sx,
      num_strata = data_1,
      var_levels = var_levels,
      n_levels = n_levels
    )
  )
}

# pseudo-inverse
p.inverse = function (m) {
  s = svd(m)
  tol = max(dim(m)) * max(s$d) * .Machine$double.eps
  Positive = s$d > tol
  d = s$d[Positive]
  u = s$u[, Positive, drop = FALSE]
  v = s$v[, Positive, drop = FALSE]

  if (length(d) == 0) {
    return(array(0, dim(m)[2:1]))
  } else{
    return(v %*% (1 / d * t(u)))
  }
}

# Loss and Mahalanobis
Imb.m <- function(delta, Fn) {
  n = nrow(Fn)
  if (length(unique(delta)) == 1) {
    mahal = "at least one observation \n              per group is needed"
    Fc = p.inverse(t(Fn) %*% Fn)
    bt.CA = t(2 * delta - 1) %*% Fn
    loss = (bt.CA %*% Fc %*% t(bt.CA))[1]
  } else if (ncol(Fn) == 1) {
    mahal = "at least one non-unique \n               variable is needed"
    Fc = p.inverse(t(Fn) %*% Fn)
    bt.CA = t(2 * delta - 1) %*% Fn
    loss = (bt.CA %*% Fc %*% t(bt.CA))[1]
  } else{
    Fc = p.inverse(t(Fn) %*% Fn)
    bt.CA = t(2 * delta - 1) %*% Fn
    loss = (bt.CA %*% Fc %*% t(bt.CA))[1]
    Z = as.matrix(Fn[, -1])
    AB = apply(as.matrix(Z[delta == 1, ]), 2, mean) - apply(as.matrix(Z[delta == 0, ]), 2, mean) #apply means!
    ccl = p.inverse(cov(Z))
    tot.a = sum(delta == 1)
    mahal = ((t(AB) %*% ccl %*% AB)[1]) * (n - tot.a) * tot.a / n
  }
  return(list(
    Loss = loss,
    Mahal = mahal,
    overall.imb = 2 * sum(delta) - n
  ))
}

# Loss and Mahalanobis for MC replicates
Imb.m.sim <- function(delta, Fn) {
  n = nrow(Fn)
  if (n < 2) {
    loss = NA
    mahal = NA
  } else if (length(unique(delta)) == 1) {
    mahal = NA
    Fc = p.inverse(t(Fn) %*% Fn)
    bt.CA = t(2 * delta - 1) %*% Fn
    loss = (bt.CA %*% Fc %*% t(bt.CA))[1]
  } else if (ncol(Fn) == 1) {
    mahal = NA
    Fc = p.inverse(t(Fn) %*% Fn)
    bt.CA = t(2 * delta - 1) %*% Fn
    loss = (bt.CA %*% Fc %*% t(bt.CA))[1]
  } else{
    Fc = p.inverse(t(Fn) %*% Fn)
    bt.CA = t(2 * delta - 1) %*% Fn
    loss = (bt.CA %*% Fc %*% t(bt.CA))[1]
    Z = as.matrix(Fn[, -1])
    # if(ncol(as.matrix(Z))==1){
    #   mahal = NA
    # } else {
      AB = apply(as.matrix(Z[delta == 1, ]), 2, mean, na.rm = TRUE) - apply(as.matrix(Z[delta == 0, ]), 2, mean, na.rm = TRUE)
      ccl = p.inverse(cov(Z))
      tot.a = sum(delta == 1)
      mahal = ((t(AB) %*% ccl %*% AB)[1]) * (n - tot.a) * tot.a / n
    #}
  }
  return(c(
    Loss = loss,
    Mahal = mahal,
    overall.imb = 2 * sum(delta) - n
  ))
}


# print.covadap
print_covadap <- function(res) {
  tt1 = res$Strata.measures
  if(nrow(tt1) == 1){
    n.d = sum(pmax(max(nchar(
      sapply(tt1, as.character)
    )),
    nchar(colnames(tt1)))) + nchar(nrow(tt1)) + ncol(tt1) + 1
  }else{
    n.d = sum(pmax(apply(nchar(
      sapply(tt1, as.character)
    ), 2, max),
    nchar(colnames(tt1)))) + nchar(nrow(tt1)) + ncol(tt1) + 1
  }
  end.line = paste(rep("_", n.d),  collapse = "")
  cat(end.line, sep = "\n")

  cat(
    unlist(res$summary.info[-6]),
    sep = "\n",
    labels = paste0(names(res$summary.info[-6]), ":"),
    fill = length(res$summary.info[-6])
  )
  cat(end.line, sep = "\n")
  cat("IMBALANCE MEASURES", sep = "\n")
  cat(" ", sep = "\n")

  r.imb = res$Imbalances
  #index = sapply(r.imb, is.numeric)
  #r.imb[index] = round(unlist(r.imb[index]), 2)
  indx.1 = sapply(r.imb$Imb.measures, is.numeric)
  r.imb$Imb.measures[indx.1] = sapply(r.imb$Imb.measures[indx.1], round, digits = 2)
  # a = list(lapply(r.imb$Imb.measures, round, digits=2), r.imb$Overall.imb)
  a = list(r.imb$Imb.measures, r.imb$Overall.imb)
  cat(paste0(c("Loss", "Mahalanobis", "Overall Imb"), ": ", unlist(a)), sep = "\n", fill = 1)
  cat(end.line, sep = "\n")
  cat("STRATA IMBALANCES", sep = "\n")
  cat(" ", sep = "\n")
  print(tt1)
  cat(end.line, sep = "\n")
  cat("WITHIN COVARIATE IMBALANCES", sep = "\n")
  cat(" ", sep = "\n")
  cat(
    unlist(res$Imbalances$Within.cov),
    labels = paste0(names(unlist(res$Imbalances$Within.cov)), ":"),
    fill = length(res$Imbalances$Within.cov),
    sep = "\n"
  )
}


# print for simulation results
print_sim <- function(res) {
  tt1 = cbind(
    res$Imbalances$obs.strata,
    N.strata = round(apply(res$Imbalances$strata.N, 2, mean, na.rm = TRUE),2),
    A.strata = round(apply(res$Imbalances$strata.A, 2, mean, na.rm = TRUE),2),
    D.strata = round(apply(abs(res$Imbalances$strata.imb), 2, mean, na.rm = TRUE),2)
  )
  if(nrow(tt1)== 1){
    n.d = sum(pmax(max(nchar(
      sapply(tt1, as.character)
    )),
    nchar(colnames(tt1)))) + nchar(nrow(tt1)) + ncol(tt1) + 1
  }else{
  n.d = sum(pmax(apply(nchar(
    sapply(tt1, as.character)
  ), 2, max),
  nchar(colnames(tt1)))) + nchar(nrow(tt1)) + ncol(tt1) + 1
  }
  end.line = paste(rep("_", n.d),  collapse = "")
  cat(end.line, sep = "\n")
  cat(
    unlist(res$summary.info),
    sep = "\n",
    labels = paste0(names(res$summary.info), ":"),
    fill = length(res$summary.info)
  )

  cat(end.line, sep = "\n")
  cat("MEAN ABSOLUTE IMBALANCE MEASURES",
      sep = "\n",
      fill = 1)
  cat(" ", sep = "\n")
  cat(paste0(c("Loss", "Mahalanobis", "Overall Imb"), ": ",
             round(apply(
               abs(res$Imbalances$Imb.measures), 2, mean, na.rm = TRUE
             ), 3)),
      sep = "\n", fill = 1)
  cat(end.line, sep = "\n")
  cat("MEAN ABSOLUTE STRATA IMBALANCES", sep = "\n")
  cat(" ", sep = "\n")
  print(tt1)
  cat(end.line, sep = "\n")
  cat("MEAN ABSOLUTE WITHIN COVARIATE IMBALANCES", sep = "\n")
  cat(" ", sep = "\n")
  cat(
    round(apply(abs(res$Imbalances$within.imb), 2, mean, na.rm = TRUE),2),
    labels = paste0(colnames(res$Imbalances$within.imb), ":"),
    fill = ncol(res$Imbalances$within.imb),
    sep = "\n"
  )
}

# Print for mixed profile and simulation
print_mixed.covadap <- function(res, case) {

  if(case == 1){

    n.c = max( nchar( unlist(res$summary.info[-8]) ) + nchar(names(res$summary.info[-8])),
               nchar("WITHIN COVARIATE IMBALANCES"), nchar( colnames(res$Imbalances$Within.cov)) + 6  )

    end.line = paste(rep("_", n.c + 3),  collapse = "")
    cat(end.line, sep = "\n")

    cat(
      unlist(res$summary.info[-8]),
      sep = "\n",
      labels = paste0(names(res$summary.info[-8]), ":"),
      fill = length(res$summary.info[-8])
    )

    cat(end.line, sep = "\n")

    cat("IMBALANCE MEASURES", sep = "\n")

    cat(" ", sep = "\n")

    r.imb = res$Imbalances.summary

    index = sapply(r.imb, is.numeric)

    r.imb[index] = round(unlist(r.imb[index]), 2)

    a <- list(lapply(r.imb$Loss, round, digits=2),
              lapply(r.imb$Mahal, round, digits=2),
              r.imb$overall.imb)


    cat(paste0(c("Loss", "Mahalanobis", "Overall Imb"), ": ",unlist(a)), sep = "\n", fill = 1)

    cat(end.line, sep = "\n")

    cat("WITHIN COVARIATE IMBALANCES", sep = "\n")

    cat(" ", sep = "\n")

    if(!is.null(res$summary.info$n_categorical_variables) && res$summary.info$n_categorical_variables==1){
      cat(
        unlist(res$Imbalances$Within.cov),
        labels = paste0(res$summary.info$cov_levels_names[-1], ":"),
        fill = length(res$Imbalances$Within.cov),
        sep = "\n"
      )
    } else {
      cat(" ", sep = "\n")
      cat(
        unlist(res$Imbalances$Within.cov),
        labels = paste0(names(unlist(res$Imbalances$Within.cov)), ":"),
        fill = length(res$Imbalances$Within.cov),
        sep = "\n"
      )
    }

    cat(end.line, sep = "\n")
    cat(" ", sep = "\n")

    cat("DIFFERENCE IN MEANS", sep = "\n")

    cat(" ", sep = "\n")

    print(round(res$diff_mean, 3))

  }else{

    n.c = max( nchar( unlist(res$summary.info[-8]) ) + nchar(names(res$summary.info[-8])),
               nchar("DIFFERENCE IN MEANS")  )

    end.line = paste(rep("_", n.c + 3),  collapse = "")

    cat(end.line, sep = "\n")

    cat(
      unlist(res$summary.info),
      sep = "\n",
      labels = paste0(names(res$summary.info), ":"),
      fill = length(res$summary.info)
    )

    cat(end.line, sep = "\n")
    cat(" ", sep = "\n")

    cat("IMBALANCE MEASURES", sep = "\n")

    cat(" ", sep = "\n")

    r.imb = res$Imbalances
    index = sapply(r.imb, is.numeric)
    r.imb[index] = round(unlist(r.imb[index]), 2)
    a <- list(lapply(r.imb$Imb.measures, round, digits=2), r.imb$Overall.imb)

    cat(paste0(c("Loss", "Mahalanobis", "Overall Imb"), ": ",unlist(a)), sep = "\n", fill = 1)

    cat(end.line, sep = "\n")

    cat(" ", sep = "\n")

    cat("DIFFERENCE IN MEANS", sep = "\n")

    cat(" ", sep = "\n")

    print(round(res$diff_mean, 3))

  }

}


print_mixed.sim <- function(res, case) {

  if(case == 1){


    n.c = max( nchar( unlist(res$summary.info[-8]) ) + nchar(names(res$summary.info[-8])),

               nchar("WITHIN COVARIATE IMBALANCES"), nchar( colnames(res$Imbalances$within.imb)) + 6  )

    end.line = paste(rep("_", n.c + 3),  collapse = "")

    cat(paste(rep("_", n.c + 3),  collapse = ""))

    cat(" ", sep = "\n")

    cat(

      unlist(res$summary.info[-8]),

      sep = "\n",

      labels = paste0(names(res$summary.info[-8]), ":"),

      fill = length(res$summary.info[-8])

    )

    cat(end.line, sep = "\n")

    cat("IMBALANCE MEASURES", sep = "\n")

    cat(" ", sep = "\n")

    r.imb = round(apply(abs(res$Imbalances$Imb.measures), 2, mean, na.rm = TRUE),2)

    cat(paste0(c("Loss", "Mahalanobis", "Overall Imb"), ": ",unlist(r.imb)), sep = "\n", fill = 1)

    cat(end.line, sep = "\n")
    cat("WITHIN COVARIATE IMBALANCES", sep = "\n")

    cat(" ", sep = "\n")
    cat(
      round(apply(abs(res$Imbalances$within.imb), 2, mean, na.rm = TRUE),2),
      labels = paste0(colnames(res$Imbalances$within.imb), ":"),
      fill = ncol(res$Imbalances$within.imb),
      sep = "\n"
    )

    cat(end.line, sep = "\n")
    cat("DIFFERENCE IN MEANS", sep = "\n")

    cat(" ", sep = "\n")

    print(round(apply(res$Imbalances$diff_mean, 2, mean, na.rm = TRUE), 3))

  }else{

    n.c = max( nchar( unlist(res$summary.info[-8]) ) + nchar(names(res$summary.info[-8])),

               nchar(" DIFFERENCE IN MEANS")  )


    end.line = paste(rep("_", n.c + 3),  collapse = "")
    cat(end.line, sep = "\n")

    cat(

      unlist(res$summary.info),

      sep = "\n",

      labels = paste0(names(res$summary.info), ":"),

      fill = length(res$summary.info)

    )

    cat(end.line, sep = "\n")

    cat("IMBALANCE MEASURES", sep = "\n")

    cat(" ", sep = "\n")

    r.imb = round(apply(abs(res$Imbalances$Imb.measures), 2, mean, na.rm = TRUE),2)

    cat(paste0(c("Loss", "Mahalanobis", "Overall Imb"), ": " ,unlist(r.imb)), sep = "\n", fill = 1)

    cat(end.line, sep = "\n")

    cat("DIFFERENCE IN MEANS", sep = "\n")

    cat(" ", sep = "\n")

    print(round(apply(res$Imbalances$diff_mean, 2, mean, na.rm = TRUE), 3))

  }

}
