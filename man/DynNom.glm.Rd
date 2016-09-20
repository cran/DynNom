\name{DynNom.glm}
\alias{DynNom.glm}
\title{
Dynamic Nomograms for Generalized Linear Models
}
\description{
\code{DynNom.glm} displays the results of a \code{glm} model object as a dynamic nomogram in an 'RStudio' panel or web browser.
}
\usage{
DynNom.glm(model, data,
           clevel = 0.95, covariate = c("slider", "numeric"))
}
\arguments{
  \item{model}{
  a \code{glm} model object
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
}
\value{
A dynamic nomogram in a shiny application which recognises all the predictors in the model and uses them to build a sidebar panel. It sets up drop down menus for factors and sliders set at the mean and bounded by the range for covariates.

The individual predictions with a relative confidence interval are calculated using the \code{predict} function, displaying either graphically as an interactive plot in the \code{Graphical Summary} tab or a table in the \code{Numerical Summary} tab. A table of model output is also available in the \code{Model Summary} tab.
}
\author{
Amirhossein Jalali, Alberto Alvarez-Iglesias, John Newell
}
\seealso{
\code{\link{glm}}, \code{\link{predict.glm}}
}
\examples{
\dontrun{
# a generilized linear model
model1 <- glm(Fertility ~ Agriculture + Education + Catholic, data = swiss)
DynNom(model1, swiss, clevel = 0.9)

# a logistic regression model
data1 =as.data.frame(Titanic)
model2 <- glm(Survived ~ Age + Class + Sex, data = data1, weights = Freq,
              family = binomial("probit"))
DynNom(model2, as.data.frame(Titanic), clevel = 0.9)
}

if (interactive()) {
# a poisson regression model
model3 <- glm(event ~ mag + dist + accel, data = attenu, family = poisson)
DynNom(model3, attenu, covariate = "numeric")
}
}
\keyword{ dynamic nomograms }
\keyword{ generalized linear models }
\keyword{ individual prediction }
\keyword{ shiny }
