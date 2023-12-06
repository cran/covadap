# Summary of Covariate-Adaptive Designs
summary_covadap <- function(res){
  if(inherits(res, what = "covadap")){ #is.null(res$summary.info$n.rep)
    if(is.null(res$diff_mean)){
      print_covadap(res)
    }else{
      if(is.null(res$summary.info$n_categorical_variables)){
        print_mixed.covadap(res, case = 2)
      }else{
        print_mixed.covadap(res, case = 1)
      }
    }
  }else{# sim
    if(is.null(res$Imbalances$dif.cont)){
      print_sim(res)
    }else{
      if(is.null(res$Imbalances$within.imb)){
        print_mixed.sim(res, case = 2)
      }else{
        print_mixed.sim(res, case = 1)
      }
    }
  }
}
