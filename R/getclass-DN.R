
getclass.DN <- function(model){
  mfamily <- NA
  mclass <- attr(model, "class")[1]
  if (mclass == "coxph.null")
    stop("Error in model syntax: the model is null")

  if (!mclass %in% c("lm", "glm", "coxph", "ols", "lrm", "Glm", "cph", "gam", "Gam", "glmnet"))
    stop("Unrecognized model object type.")

  if (mclass %in% c("elnet", "lognet", "multnet", "fishnet", "coxnet", "mrelnet")){
    mclass <- "glmnet"
    mfamily <- attr(model, "class")[1]
  }

  if (mclass %in% c("glm", "Glm"))
    mfamily <- model$family$family
  if (mclass == "lrm")
    mfamily <- mclass

  list(model.class = mclass, model.family = mfamily)
}
