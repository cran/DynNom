\name{DynNom-package}
\alias{DynNom-package}
\docType{package}
\title{
A Dynamic Nomogram for Linear and Generalized Linear Models as Shiny Applications
}
\description{
The DynNom function makes it possible to display predicted values of an lm or glm model object as a dynamic nomogram that can be displayed in an R Studio panel or web browser.
}
\details{
\tabular{ll}{
Package: \tab DynNom\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-08-05\cr
License: \tab GPL-2\cr
}
}
\author{
Amirhossein Jalali, Alberto Alvarez-Iglesias, John Newell

Maintainer: Amirhossein Jalali <a.jalali2@nuigalway.ie>
}
\references{
Levens, A. S. Nomography (2nd Edition). Fearon Publishers, Belmont, CA, USA, 1971. Banks, J. 2006. Nomograms. Encyclopedia of Statistical Sciences. 8. \cr
Easy web applications in R. \url{http://shiny.rstudio.com}
}
\section{Please cite as:}{Jalali, A., Alvarez-Iglesias, A., Newell, J. (2014). A Dynamic Nomogram for Linear and Generalized Linear Models as Shiny Applications. R package version 1.0.
}
\examples{
# a simple linear regression model
model1 <- lm(uptake ~ Plant + conc + Plant * conc, data = CO2)
\donttest{DynNom(model1, CO2)}

# a logistic regression model
swiss$ftl <- as.numeric(swiss$Fertility > 70)
model2 <- glm(ftl ~ Agriculture + Education + Catholic, data = swiss, family = binomial)
\donttest{DynNom(model2, swiss, clevel = 0.9)}

# a poisson regression model
model3 <- glm(event ~ mag + station + dist + accel, data = attenu, family = poisson)
\donttest{DynNom(model3, attenu, covariate = "numeric")}
}