function(train, test, model_meth){

  # model selection
  if(model_meth == "randomForest"){
    mdl <- randomForest(train, 1, RCF_dev, baseline_fit, as.numeric(RCF_rad))
  } else if(model_meth == "svm"){
    bl_rmv <- t(apply(spectra, 1, savgol, fl = as.numeric(filter_length),
                      forder = as.numeric(poly_order), dorder = as.numeric(deriv_order)))
  } else if(model_meth == "adaboost"){
    temp <- spc.rubberband(spectra)
    bl_rmv <- spectra - temp
  } else if(model_meth == ""){
    bl_rmv <- poly_baseline_removal(spectra, poly.order = poly_order)
  } else if(model_meth == "randomForest"){
    bl_rmv <- spectra
  }

  return(as.data.frame(mdl))
}
