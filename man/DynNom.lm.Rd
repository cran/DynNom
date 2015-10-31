\name{DynNom.lm}
\alias{DynNom.lm}
\title{
Dynamic Nomograms for Linear Models
}
\description{
DynNom.lm() displays the results of an lm() model object as a dynamic nomogram in an 'RStudio' panel or web browser.
}
\usage{
DynNom.lm(model, data,
          clevel = 0.95, covariate = c("slider", "numeric"))
}
\arguments{
  \item{model}{
  an lm model object
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
}
\value{
A dynamic nomogram in a shiny application
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
\keyword{ prediction }
\keyword{ shiny }
