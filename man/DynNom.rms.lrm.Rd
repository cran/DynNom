\name{DynNom.lrm}
\alias{DynNom.lrm}
\title{
Dynamic Nomograms for Logistic Regression Models from the \code{rms} package
}
\description{
\code{DynNOm.lrm} displays the results of a \code{lrm} model object from the \code{rms} package as a dynamic nomogram in an 'RStudio' panel or web browser.
}
\usage{
DynNom.lrm(model, data, clevel = 0.95, m.summary = c("raw", "formatted"),
                       covariate = c("slider", "numeric"))
}
\arguments{
  \item{model}{
  a \code{lrm} model object which accepts a variety of transformation functions such as \code{asis}, \code{pol}, \code{lsp}, \code{rcs}, \code{catg}, \code{scored}, \code{strat} and \code{matrx} defined in \code{rms} package.
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
Davood Roshan, Amirhossein Jalali, Alberto Alvarez-Iglesias, John Newell
}
\seealso{
\code{\link{lrm}}, \code{\link{predict.lrm}}, \code{\link{rms}}
}
\examples{
\dontrun{
# example
n <- 1000
set.seed(17)
age            <- rnorm(n, 50, 10)
blood.pressure <- rnorm(n, 120, 15)
cholesterol    <- rnorm(n, 200, 25)
sex            <- factor(sample(c('female', 'male'), n, TRUE))
label(age)            <- 'Age'      # label is in Hmisc
label(cholesterol)    <- 'Total Cholesterol'
label(blood.pressure) <- 'Systolic Blood Pressure'
label(sex)            <- 'Sex'
units(cholesterol)    <- 'mg/dl'
units(blood.pressure) <- 'mmHg'

ch <- cut2(cholesterol, g = 40, levels.mean = TRUE)

d <- data.frame(age = seq(0, 90, by = 10))

L <- .4 * (sex == 'male') + .045 * (age - 50) +
  (log(cholesterol - 10) - 5.2) * ( -2 * (sex == 'female') + 2 * (sex == 'male'))
y <- ifelse(runif(n) < plogis(L), 1, 0)
cholesterol[1:3] <- NA

ddist <- datadist(age, blood.pressure, cholesterol, sex)
options(datadist = 'ddist')

data = data.frame(y = y, blood.pressure = blood.pressure, sex = sex, age = age,
cholesterol = cholesterol)
model <- lrm(y ~ blood.pressure + sex * (age + rcs(cholesterol, 4)),
             x = TRUE, y = TRUE, m.summary = "formatted")

DynNom.lrm(model, data, m.summary = "formatted")
}

if (interactive()) {
fit <- lrm(formula = vs ~ wt + disp, data = mtcars)
DynNom.lrm(fit, mtcars, clevel = 0.9)
}
}
\keyword{ dynamic nomograms }
\keyword{ logistic regression models }
\keyword{ individual prediction }
\keyword{ shiny }
