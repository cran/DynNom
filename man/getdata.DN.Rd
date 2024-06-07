\name{getdata.DN}
\alias{getdata.DN}
\title{
  Extract dataset from a model object
}
\description{
  \code{getdata.DN} extracts dataset that was used to produce the model object (supported in \code{DynNom}).
}
\usage{
getdata.DN(model)
}
\arguments{
  \item{model}{
  an \code{lm}, \code{glm}, \code{coxph}, \code{ols}, \code{Glm}, \code{lrm}, \code{cph} or \code{mgcv::gam} model objects.
}
}
\value{
  A data.frame containing the dataset used in the fitted model object.
}
\seealso{
\code{\link{DynNom}}, \code{\link{DNbuilder}}
}
\examples{
fit1 <- glm(Survived ~ Age + Class + Sex, data = as.data.frame(Titanic),
  weights = Freq, family = binomial("probit"))
getdata.DN(fit1)

library(survival)
fit2 <- coxph(Surv(time, status) ~ age + strata(sex) + ph.ecog, data = lung)
getdata.DN(fit2)
}
