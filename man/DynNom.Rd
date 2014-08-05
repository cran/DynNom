\name{DynNom}
\alias{DynNom}
\title{
A Dynamic Nomogram for Linear and Generalized Linear Models as Shiny Applications
}
\description{
The DynNom function makes it possible to display predicted values of an lm or glm model object as a dynamic nomogram that can be displayed in an R Studio panel or web browser.
}
\usage{
DynNom(model, data, clevel = 0.95, covariate = c("slider", "numeric"))
}
\arguments{
  \item{model}{
an lm or glm model object
}
  \item{data}{
dataframe containing the accompanying data
}
\item{clevel}{
confidence level required
}
\item{covariate}{
option to choose the value(s) of the covariate(s) using a slider or as a numerical input
}
}
\details{
This function takes as input a model object created using the lm or glm function (accompanied by the dataset) and automatically generates a dynamic nomogram as a shiny application \cr Predicted values, and accompanying interval estimates, are calculated using \code{\link{predict.lm}} or \code{\link{predict.glm}} functions on the linear scale and then transformed based on the link function used.
}
\value{
A shiny application
}
\author{
Amirhossein Jalali, Alberto Alvarez-Iglesias, John Newell
}
\references{
Levens, A. S. Nomography (2nd Edition). Fearon Publishers, Belmont, CA, USA, 1971. Banks, J. 2006. Nomograms. Encyclopedia of Statistical Sciences. 8. \cr
Easy web applications in R. \url{http://shiny.rstudio.com}
} 
\section{Please cite as:}{Jalali, A., Alvarez-Iglesias, A., Newell, J. (2014). A Dynamic Nomogram for Linear and Generalized Linear Models as Shiny Applications. R package version 1.0.
}
\seealso{
  \code{\link{glm}}, \code{\link{predict.lm}}, \code{\link{predict.glm}}
}
\examples{
\dontrun{
# a simple linear regression model
model1 <- lm(uptake ~ Plant + conc + Plant * conc, data = CO2)
DynNom(model1, CO2)

# a logistic regression model
swiss$ftl <- as.numeric(swiss$Fertility > 70)
model2 <- glm(ftl ~ Agriculture + Education + Catholic, data = swiss, family = binomial)
DynNom(model2, swiss, clevel = 0.9)
}

if (interactive()) {
# a poisson regression model
model3 <- glm(event ~ mag + station + dist + accel, data = attenu, family = poisson)
DynNom(model3, attenu, covariate = "numeric")
}
}