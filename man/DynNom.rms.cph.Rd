\name{DynNom.cph}
\alias{DynNom.cph}
\title{
Dynamic Nomograms for Cox Proportional Hazards Models and Extensions
}
\description{
\code{DynNom.cph} displays the results of a \code{cph} model object from \code{rms} package as a dynamic nomogram in an 'RStudio' panel or web browser.
}
\usage{
DynNom.cph(model, data,
             clevel = 0.95, covariate = c("slider", "numeric"),
             ptype = c("st", "1-st"))
}
\arguments{
  \item{model}{
  a \code{cph} model object which accepts a variety of transformation functions such as \code{asis}, \code{pol}, \code{lsp}, \code{rcs}, \code{catg}, \code{scored}, \code{strat} and \code{matrx} as defined in the \code{rms} package.
}
  \item{data}{
  dataframe containing the accompanying data
}
  \item{clevel}{
  confidence level required
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
Davood Roshan, Amirhossein Jalali, Alberto Alvarez-Iglesias, John Newell
}
\seealso{
\code{\link{cph}}, \code{\link{predict.cph}}, \code{\link{rms}}
}
\examples{
\dontrun{
# example 1
data = veteran
model1 <- cph((Surv(log(time), status)) ~ rcs(age, 4) * strat(trt) +
              diagtime * strat(prior) + lsp(karno, 60), data = veteran)
model1 <- update(model, x = T, y = T, surv = T)
DynNom.cph(model1, data)
DynNom(model1, data, ptype = "1-st")

# example 2
data(lung)
Surv = Surv(time, status)
model2 <- cph(Surv ~ age + strat(sex) + ph.ecog , data = lung)
DynNom.cph(model2, lung)
}

if (interactive()) {
data.ovary <- ovarian
data.ovary$resid.ds <- as.factor(data.ovary$resid.ds)
levels(data.ovary$resid.ds) <- c("no", "yes")
data.ovary$rx <- as.factor(data.ovary$rx)
data.ovary$ecog.ps <- as.factor(data.ovary$ecog.ps)

model3 <- cph(Surv(futime, fustat) ~ age + resid.ds * rx + ecog.ps, data = data.ovary)
DynNom(model3, data.ovary)
}
}
\keyword{ dynamic nomograms }
\keyword{ cox proportional hazards model }
\keyword{ individual prediction }
\keyword{ shiny }
