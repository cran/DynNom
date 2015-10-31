\name{DynNom.glm}
\alias{DynNom.glm}
\title{
Dynamic Nomograms for Generalized Linear Models
}
\description{
DynNom.glm() displays the results of a glm() model object as a dynamic nomogram in an 'RStudio' panel or web browser.
}
\usage{
DynNom.glm(model, data,
           clevel = 0.95, covariate = c("slider", "numeric"))
}
\arguments{
  \item{model}{
  a glm model object
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
\keyword{ prediction }
\keyword{ shiny }
