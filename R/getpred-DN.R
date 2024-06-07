
getpred.DN <- function(model, newd, set.rms=F){
  inrange <- T
  mclass <- getclass.DN(model)$model.class

  if (!mclass %in% c("lm", "glm", "coxph", "ols", "lrm", "Glm", "cph", "gam", "glmnet"))
    stop("Unrecognized model object type.")

  if (mclass %in% c("ols", "lrm", "Glm", "cph")){
    if (set.rms==T){
      model <- update(model, x=T, y=T, data=data)
    } else{
      model <- update(model, x=T, y=T)
    }
  }

  if (mclass %in% c("ols", "lrm", "Glm")){
    m.pred <- predict(model, newdata = newd, se.fit = TRUE)
    mpred <- m.pred$linear.predictors
    se.pred <- m.pred$se.fit[[1]]
  }

  if (mclass %in% c("lm", "glm", "gam")){
    mpred <- broom::augment(x=model, newdata = newd, se_fit=T)$.fitted
    se.pred <- broom::augment(x=model, newdata = newd, se_fit=T)$.se.fit
  }

  if (mclass %in% c("coxph")){
    mpred <- broom::augment(x=model, newdata = newd, se_fit=T, type.predict = "expected")$.fitted
    se.pred <- broom::augment(x=model, newdata = newd, se_fit=T, type.predict = "expected")$.se.fit
  }

  if (mclass %in% c("cph")){
    strata.l <- levels(model$strata)
    if (length(model$strata) != length(levels(attr(predict(model, newd, type='x', expand.na=FALSE), 'strata')))){
      levels(model$strata) <- levels(attr(predict(model, newd, type='x', expand.na=FALSE), 'strata'))
    }
    m.pred <- suppressWarnings({ survest(model, newdata = newd, times = newd[,all.vars(model$terms)[1]]) })
    mpred <- -log(m.pred$surv)
    se.pred <- m.pred$std.err

    if (mpred == 0){
      inrange = F
    }
  }
  list(pred = mpred, SEpred = se.pred, InRange = inrange)
}
