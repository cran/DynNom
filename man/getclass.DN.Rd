\name{getclass.DN}
\alias{getclass.DN}
\title{
  Extract class and family of a model object
}
\description{
  \code{getclass.DN} extracts class and family of a model object (supported in \code{DynNom}).
}
\usage{
getclass.DN(model)
}
\arguments{
  \item{model}{
  an \code{lm}, \code{glm}, \code{coxph}, \code{ols}, \code{Glm}, \code{lrm}, \code{cph} or \code{mgcv::gam} model objects.
}
}
\value{
  A list including the model class and the family name of the model (if relevant).
}
\seealso{
\code{\link{DynNom}}, \code{\link{DNbuilder}}
}
\examples{
fit1 <- glm(Survived ~ Age + Class + Sex, data = as.data.frame(Titanic),
  weights = Freq, family = binomial("probit"))
getclass.DN(fit1)

library(survival)
fit2 <- coxph(Surv(time, status) ~ age + strata(sex) + ph.ecog, data = lung)
getclass.DN(fit2)
}
