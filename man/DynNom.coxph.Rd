\name{DynNom.coxph}
\alias{DynNom.coxph}
\title{
Dynamic Nomograms for Proportional Hazards Models
}
\description{
\code{DynNom.coxph} displays the results of a \code{coxph} model object as a dynamic nomogram in an 'RStudio' panel or web browser.
}
\usage{
DynNom.coxph(model, data, clevel = 0.95, m.summary = c("raw", "formatted"),
                   covariate = c("slider", "numeric"), ptype = c("st", "1-st"))
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
  \item{m.summary}{
  The option to choose the format of the model output in the 'Summary Model' tab. If "raw" (the default) is chosen the result of summary(model) will be display while if "formatted" is chosen the model summary using the stargazer package will be displayed.
}
  \item{covariate}{
  The option to choose the type of covariate(s) input control widget for numeric values. If "slider" (the default) is chosen a shiny application with slider control widgets are used while if "numeric" is chosen numeric values input controls will be displayed.
}
  \item{ptype}{
  If "st" (the default) is chosen, a plot of the estimated survivor function, S(t),  is displayed.  If "1-st" is chosen a plot of 1- S(t) is displayed.
}
}
\value{
A dynamic nomogram in a shiny application which recognises all the predictors in the model and uses them to build a sidebar panel. It sets up drop down menus for factors and sliders set at the mean and bounded by the range for covariates.

The individual predictions with a relative confidence interval are calculated using the \code{predict} function, displaying graphically as either the Kaplan-Meier in the \code{Estimated S(t)} tab or the \code{Predicted Survival} tab. Table of individual predictions and model output are available in the \code{Numerical Summary} and \code{Model Summary} tabs, respectively.
}
\author{
Amirhossein Jalali, Davood Roshan, Alberto Alvarez-Iglesias, John Newell
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

model1 <- coxph(Surv(time, status) ~ age + sex + disease, data.kidney)
DynNom(model1, data.kidney)
DynNom(model1, data.kidney, ptype = "1-st")

# a cox model including a strata term
data(lung)
model2 <- coxph(Surv(time, status) ~ age + strata(sex) + ph.ecog , data = lung)
DynNom(model2, lung)
}

if (interactive()) {
data.ovary <- ovarian
data.ovary$resid.ds <- as.factor(data.ovary$resid.ds)
levels(data.ovary$resid.ds) <- c("no", "yes")
data.ovary$rx <- as.factor(data.ovary$rx)
data.ovary$ecog.ps <- as.factor(data.ovary$ecog.ps)

model3 <- coxph(Surv(futime, fustat) ~ age + resid.ds * rx + ecog.ps, data = data.ovary)
DynNom(model3, data.ovary)
}
}
\keyword{ dynamic nomograms }
\keyword{ Cox Proportional Hazards Model }
\keyword{ individual prediction }
\keyword{ shiny }
