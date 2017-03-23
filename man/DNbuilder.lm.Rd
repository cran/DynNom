\name{DNbuilder.lm}
\alias{DNbuilder.lm}
\title{
Shiny code for Dynamic Nomograms
}
\description{
\code{DNbuilder.lm} provides required scripts to deploy an \code{lm} model object as a dynamic nomogram on a server on the web such as the \url{http://shinyapps.io}.
}
\usage{
DNbuilder.lm(model, data, clevel = 0.95, m.summary = c("raw", "formatted"),
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
A new folder in the current working directory called DynNomapp which contains all the required scripts to deploy this dynamic nomogram on a server on the web such as the \url{http://shinyapps.io}. This folder includes \code{ui.R}, \code{server.R} and \code{global.R} script files needed to build the application and \code{dataset.rds} which is the accompanying dataset and a user guide text file called \code{README.txt} which explains how to deploy the app using all these objects.
}
\author{
Amirhossein Jalali, Alberto Alvarez-Iglesias, John Newell
}
\seealso{
\link{lm}, \link{DynNom}, \link{DynNom.lm}
}
\examples{
\dontrun{
# a linear regression model
model1 <- lm(Fertility ~ Agriculture + Education + Catholic , data = swiss)
DNbuilder.lm(model1, swiss)

model2 <- lm(uptake ~ Plant + conc + Plant * conc, data = CO2)
DNbuilder(model2, CO2)
}

if (interactive()) {
data1 <- data.frame(state.x77)
fit1 <- lm(formula = Life.Exp ~ ., data = data1)
DNbuilder(fit1, data1)
}
}
\keyword{ dynamic nomograms }
\keyword{ linear models }
\keyword{ individual prediction }
\keyword{ shiny }
