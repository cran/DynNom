\name{DynNom.ols}
\alias{DynNom.ols}
\title{
Dynamic Nomograms for Linear Models from the \code{rms} package
}
\description{
\code{DynNom.ols} displays the results of an \code{ols} model object from the \code{rms} package as a dynamic nomogram in an 'RStudio' panel or web browser.
}
\usage{
DynNom.ols(model, data, clevel = 0.95, m.summary = c("raw", "formatted"),
                       covariate = c("slider", "numeric"))
}
\arguments{
  \item{model}{
  an \code{ols} model object which accepts a variety of transformation functions such as \code{asis}, \code{pol}, \code{lsp}, \code{rcs}, \code{catg}, \code{scored}, \code{strat} and \code{matrx} defined in \code{rms} package.
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
\code{\link{ols}}, \code{\link{predict.ols}}, \code{\link{rms}}
}
\examples{
\dontrun{
# example 1
x1 <- runif(200)
x2 <- runif(200)
x3 <- runif(200)
x4 <- runif(200)
y <- x1 + x2 + rnorm(200)
data = data.frame(x1, x2, x3, x4, y)
f    <- ols(y ~ rcs(x1, 4) + x2 + x3 + x4)
DynNom.ols(f, data)

# example 2
data1 = as.data.frame(Titanic)
year <- sample(c(1:5), 32, replace = TRUE)
data <- data.frame(year, data1)
model <- ols(year ~ Age + Class + Sex, data = data, weights = Freq)
DynNom.ols(model, data)
}

if (interactive()) {
data1 <- data.frame(state.x77)
fit1 <- ols(Life.Exp ~ Population + Income + Murder + Frost , data = data1)
DynNom(fit1, data1)
}
}
\keyword{ dynamic nomograms }
\keyword{ linear models }
\keyword{ individual prediction }
\keyword{ shiny }
