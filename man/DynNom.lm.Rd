\name{DynNom.lm}
\alias{DynNom.lm}
\title{
Dynamic Nomograms for Linear Models
}
\description{
\code{DynNom.lm} displays the results of an \code{lm} model object as a dynamic nomogram in an 'RStudio' panel or web browser.
}
\usage{
DynNom.lm(model, data, clevel = 0.95, m.summary = c("raw", "formatted"),
                      covariate = c("slider", "numeric"))
}
\arguments{
  \item{model}{
  an \code{lm} model object
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
}
\value{
A dynamic nomogram in a shiny application which recognises all the predictors in the model and uses them to build a sidebar panel. It sets up drop down menus for factors and sliders set at the mean and bounded by the range for covariates.

The individual predictions with a relative confidence interval are calculated using the \code{predict} function, displaying either graphically as an interactive plot in the \code{Graphical Summary} tab or a table in the \code{Numerical Summary} tab. A table of model output is also available in the \code{Model Summary} tab.
}
\author{
Amirhossein Jalali, Alberto Alvarez-Iglesias, John Newell
}
\seealso{
\code{\link{lm}}, \code{\link{predict.lm}}
}
\examples{
\dontrun{
# a linear regression model
model1 <- lm(Fertility ~ Agriculture + Education + Catholic , data = swiss)
DynNom(model1, swiss)

model2 <- lm(uptake ~ Plant + conc + Plant * conc, data = CO2)
DynNom(model2, CO2)
}

if (interactive()) {
data1 <- data.frame(state.x77)
fit1 <- lm(formula = Life.Exp ~ ., data = data1)
DynNom(fit1, data1)
}
}
\keyword{ dynamic nomograms }
\keyword{ linear models }
\keyword{ individual prediction }
\keyword{ shiny }
