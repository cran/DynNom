\name{getpred.DN}
\alias{getpred.DN}
\title{
  Extract predictions from a Model Object
}
\description{
  \code{getpred.DN} extracts class, family and inverse of link function from a model object (supported in \code{DynNom}).
}
\usage{
getpred.DN(model, newd, set.rms=F)
}
\arguments{
  \item{model}{
  an \code{lm}, \code{glm}, \code{coxph}, \code{ols}, \code{Glm}, \code{lrm}, \code{cph}, \code{mgcv::gam} or \code{gam::gam} model objects.
}
\item{newd}{
  a data frame of predictors for prediction
}
\item{set.rms}{
  a logical value indicating if data should be updated in the model object (required for \code{rms} model objects in \code{DNbuilder}).
}
}
\value{
  A list including the prediction (\code{pred}) and the standard error of prediction (\code{SEpred}).
}
\seealso{
\code{\link{DynNom}}, \code{\link{DNbuilder}}
}
\examples{
fit1 <- glm(Survived ~ Age + Class + Sex, data = as.data.frame(Titanic),
  weights = Freq, family = binomial("probit"))
getpred.DN(fit1, newd = data.frame(Class="1st", Sex="Male", Age="Child"))
}
