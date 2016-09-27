\name{DynNom.Glm}
\alias{DynNom.Glm}
\title{
Dynamic Nomograms for Generalized Linear Models from the \code{rms} package
}
\description{
\code{DynNom.Glm} displays the results of a \code{Glm} model object from the \code{rms} package as a dynamic nomogram in an 'RStudio' panel or web browser.
}
\usage{
DynNom.Glm(model, data,
           clevel = 0.95, covariate = c("slider", "numeric"))
}
\arguments{
  \item{model}{
  a \code{Glm} model object which accepts a variety of transformation functions such as \code{asis}, \code{pol}, \code{lsp}, \code{rcs}, \code{catg}, \code{scored}, \code{strat} and \code{matrx} defined in \code{rms} package.
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
Davood Roshan, Amirhossein Jalali, Alberto Alvarez-Iglesias, John Newell
}
\seealso{
\code{\link{Glm}}, \code{\link{predict.Glm}}, \code{\link{rms}}
}
\examples{
\dontrun{
# example 1 - a generalized linear model
set.seed(1)
x1 <- runif(200)
x2 <- sample(0:3, 200, TRUE)
x3 <- sample(0:2, 200, TRUE)

distance <- (x1 + x2 / 3 + rnorm(200)) ^ 2
d <- datadist(x1, x2)
options(datadist = "d")
data1 = data.frame(distance, x1, x2, x3)
model1 <- Glm(distance ~ x3 + rcs(x1, 4) * scored(x2), data = data1)
DynNom.Glm(model1, data1)

# example 2 - a poisson regression model
counts <- c(18, 17, 15, 20, 10, 20, 25, 13, 12)
outcome <- gl(3, 1, 9)
treatment <- gl(3, 3)
data2 = data.frame(counts, outcome, treatment)
model2 <- Glm((2 * counts) ~ outcome + treatment, family = poisson(), data = data2)
DynNom.Glm(model2, data2)
}

if (interactive()) {
# a Gamma regression model
clotting <- data.frame(
    u = c(5, 10, 15, 20, 30, 40, 60, 80, 100),
    lot1 = c(118, 58, 42, 35, 27, 25, 21, 19, 18),
    lot2 = c(69, 35, 26, 21, 18, 16, 13, 12, 12),
    cat = c(rep("A",5), rep("B",4)))
model3 <- Glm(lot1 ~ log(u) + cat, data = clotting, family = Gamma)
DynNom.Glm(model3, clotting)
}
}
\keyword{ dynamic nomograms }
\keyword{ generalized linear models }
\keyword{ individual prediction }
\keyword{ shiny }
