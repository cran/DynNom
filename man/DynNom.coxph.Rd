\name{DynNom.coxph}
\alias{DynNom.coxph}
\title{
Dynamic Nomograms for Proportional Hazard Models
}
\description{
DynNom.coxph() displays the results of a coxph() model object as a dynamic nomogram in an 'RStudio' panel or web browser.
}
\usage{
DynNom.coxph(model, data,
             clevel = 0.95, covariate = c("slider", "numeric"),
             ptype = c("st", "1-st"))
}
\arguments{
  \item{model}{
  a coxph model object
}
  \item{data}{
  dataframe containing the accompanying data
}
  \item{clevel}{
  confidence level required
}
  \item{covariate}{
  The option to choose the type of covariate(s) input control widget for numeric values. If "slider" (the default) is chosen a shiny application with slider control widget is used while if "numeric" is chosen a numeric values input control will be displayed.
}
  \item{ptype}{
  If "st" (the default) is chosen, a plot of the estimated survivor function, S(t),  is displayed.  If "1-st" is chosen a plot of 1- S(t) is displayed.
}
}
\value{
A dynamic nomogram in a shiny application
}
\seealso{
\code{\link{coxph}}, \code{\link{predict.coxph}}
}
\examples{
\dontrun{
data.kidney <- kidney
# always make sure that the categorical variables are in a factor class
data.kidney$sex <- as.factor(data.kidney$sex)
levels(data.kidney$sex) <- c("male", "female")

model1 <- coxph(Surv(time, status)~ age + sex + disease, data.kidney)
DynNom(model1, data.kidney)
DynNom(model1, data.kidney, ptype = "1-st")

# a cox model including a strata term
data(lung)
model1 <- coxph(Surv(time, status) ~ age + strata(sex) + ph.ecog , data=lung)
DynNom(model1, lung)
}

if (interactive()) {
data.ovary <- ovarian
data.ovary$resid.ds <- as.factor(data.ovary$resid.ds)
levels(data.ovary$resid.ds) <- c("no", "yes")
data.ovary$rx <- as.factor(data.ovary$rx)
data.ovary$ecog.ps <- as.factor(data.ovary$ecog.ps)

fit1 <- coxph(Surv(futime, fustat) ~ age + resid.ds *rx + ecog.ps, data = data.ovary)
DynNom(fit1, data.ovary)
}
}
\keyword{ dynamic nomograms }
\keyword{ proportional hazard models }
\keyword{ prediction }
\keyword{ shiny }
