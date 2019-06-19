
DynNom <- function(model, data = NULL, clevel = 0.95, m.summary = c("raw", "formatted"),
                   covariate = c("slider", "numeric"), ptype = c("st", "1-st"),
                   DNtitle = NULL, DNxlab = NULL, DNylab = NULL, DNlimits = NULL, KMtitle = NULL, KMxlab = NULL, KMylab = NULL) {

  mclass <- getclass.DN(model)$model.class
  mfamily <- getclass.DN(model)$model.family

  if (mclass %in% c("coxph", "cph")){
    Surv.in <- length(model$terms[[2]]) != 1
  }

  if (mclass %in% c("ols", "Glm", "lrm", "cph")){
    model <- update(model, x=T, y=T)
  }

  if (!is.data.frame(data)){
    if (any(class(try(getdata.DN(model), silent = TRUE)) == "try-error")){
      stop("Dataset needs to be provided in a data.frame format")
    } else{
      data <- getdata.DN(model)
    }
  }
  covariate <- match.arg(covariate)
  m.summary <- match.arg(m.summary)
  ptype <- match.arg(ptype)

  if (mclass %in% c("lm", "glm", "ols", "Glm", "lrm", "gam", "Gam")){
    Terms.T <- all(all.vars(model$terms) %in% names(data))
  }
  if (mclass %in% c("coxph")){
    if (Surv.in){
      Terms.T <- all(all.vars(model$terms)[-c(1:2)] %in% names(data))
    } else{
      Terms.T <- all(all.vars(model$terms)[-1] %in% names(data))
    }
  }
  if (mclass %in% c("cph")){
    Terms.T <- all(names(model$Design$units) %in% names(data))
  }
  if (!Terms.T)
    stop("Error in model syntax: some of model's terms do not match to variables' name in dataset")

  if (!is.null(DNlimits) & !length(DNlimits)==2)
    stop("A vector of 2 is required as 'DNlimits'")

  if (is.null(DNtitle))
    DNtitle <- "Dynamic Nomogram"

  if (is.null(DNxlab)){
    DNxlab <- ifelse((mclass %in% c("glm") & mfamily %in% c("binomial", "quasibinomial")) | mclass == "lrm" | mclass %in% c("coxph", "cph"),
                     "Probability", "Response variable")
  }

  if (mclass %in% c("coxph", "cph")){
    if (is.null(KMtitle)){
      if (ptype == "st") {
        KMtitle <- "Estimated Survival Probability"
      } else{
        KMtitle <- "Estimated Probability"
      }
    }
    if (is.null(KMxlab)){
      KMxlab <- "Follow Up Time"
    }
    if (is.null(KMylab)){
      if (ptype == "st") {
        KMylab <- "S(t)"
      } else{
        KMylab <- "F(t)"
      }
    }
  }

  if (mclass %in% c("lm", "glm", "ols", "Glm", "lrm", "gam", "Gam")) {
    DynNom.core(model, data, clevel, m.summary, covariate, DNtitle, DNxlab, DNylab, DNlimits)
  }
  if (mclass %in% c("coxph", "cph")) {
    DynNom.surv(model, data, clevel, m.summary, covariate, ptype, DNtitle, DNxlab, DNylab, KMtitle, KMxlab, KMylab)
  }
}
