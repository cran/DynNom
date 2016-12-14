\name{DNbuilder.glm}
\alias{DNbuilder.glm}
\title{
Shiny code for Dynamic Nomograms
}
\description{
\code{DNbuilder.glm} provides required scripts to deploy an \code{lm} model object as a dynamic nomogram on a server on the web such as the \url{http://shinyapps.io}.
}
\usage{
DNbuilder.glm(model, data,
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
A new folder in the current working directory called DynNomapp which contains all the required scripts to deploy this dynamic nomogram on a server on the web such as the \url{http://shinyapps.io}. This folder includes \code{ui.R}, \code{server.R} and \code{global.R} script files needed to build the application and \code{dataset.rds} which is the accompanying dataset and a user guide text file called \code{README.txt} which explains how to deploy the app using all these objects.
}
\author{
Amirhossein Jalali, Alberto Alvarez-Iglesias, John Newell
}
\seealso{
\link{glm}, \link{DynNom}, \link{DynNom.glm}
}
\examples{
\dontrun{
# a generilized linear model
model1 <- glm(Fertility ~ Agriculture + Education + Catholic, data = swiss)
DNbuilder.glm(model1, swiss, clevel = 0.9)

# a logistic regression model
data1 =as.data.frame(Titanic)
model2 <- glm(Survived ~ Age + Class + Sex, data = data1, weights = Freq,
              family = binomial("probit"))
DNbuilder(model2, as.data.frame(Titanic), clevel = 0.9)
}

if (interactive()) {
# a poisson regression model
model3 <- glm(event ~ mag + dist + accel, data = attenu, family = poisson)
DNbuilder(model3, attenu, covariate = "numeric")
}
}
\keyword{ dynamic nomograms }
\keyword{ generalized linear models }
\keyword{ individual prediction }
\keyword{ shiny }
