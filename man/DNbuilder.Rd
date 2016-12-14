\name{DNbuilder}
\alias{DNbuilder}
\title{
Shiny code for Dynamic Nomograms
}
\description{
  \code{DNbuilder} is a generic function to build dynamic nomograms and provide the required scripts for deploying them on a server on the web such as the \url{http://shinyapps.io}. \code{DNbuilder} supports \code{lm}, \code{glm}, \code{coxph} model objects.
}
\usage{
DNbuilder(model, data,
       clevel = 0.95, covariate = c("slider", "numeric"),
       ptype = c("st", "1-st"))
}
\arguments{
  \item{model}{
  an \code{lm}, \code{glm} or \code{coxph} model object
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
  \item{ptype}{
  This plot type option relates to coxph objects only. If "st" (the default) is chosen, a plot of the estimated survivor function, S(t),  is displayed.  If "1-st" is chosen a plot of 1- S(t) is displayed.
}
}
\value{
  A new folder in the current working directory called DynNomapp which contains all the required scripts to deploy this dynamic nomogram on a server on the web such as the \url{http://shinyapps.io}. This folder includes \code{ui.R}, \code{server.R} and \code{global.R} script files needed to build the application and \code{dataset.rds} which is the accompanying dataset and a user guide text file called \code{README.txt} which explains how to deploy the app using all these objects.
}
\references{
Banks, J. 2006. Nomograms. Encyclopedia of Statistical Sciences. 8. \cr
Easy web applications in R. \url{http://shiny.rstudio.com}
}
\author{
Amirhossein Jalali, Alberto Alvarez-Iglesias, John Newell
}
\seealso{
\link{lm}, \link{glm}, \link{coxph}, \link{DynNom}, \link{DynNom.lm}, \link{DynNom.glm}, \link{DynNom.coxph}
}
\examples{
\dontrun{
# simple linear regression models
model1 <- lm(uptake ~ Plant + conc + Plant * conc, data = CO2)
DNbuilder(model1, CO2)

# Generalized regression models
data1 =as.data.frame(Titanic)
model2 <- glm(Survived ~ Age + Class + Sex, data = data1, weights = Freq,
              family = binomial("probit"))
DNbuilder(model2, data1, clevel = 0.9)

# a proportional hazard model
data.kidney <- kidney
# always make sure that the categorical variables are in a factor class
data.kidney$sex <- as.factor(data.kidney$sex)
levels(data.kidney$sex) <- c("male", "female")

model3 <- coxph(Surv(time, status) ~ age + sex + disease, data.kidney)
DNbuilder(model3, data.kidney)
DNbuilder(model3, data.kidney, ptype = "1-st")
}

if (interactive()) {
# a poisson regression model
model4 <- glm(event ~ mag + station + dist + accel, data = attenu, family = poisson)
DynNom(model4, attenu, covariate = "numeric")
}
}
\keyword{ dynamic nomograms }
\keyword{ individual prediction }
\keyword{ shiny }
