\name{DynNom}
\alias{DynNom}
\title{
Dynamic Nomograms for Linear, Generalized Linear and proportional hazard Models
}
\description{
DynNom() is a generic function for displaying the results of an lm(), glm() or coxph() model object as a dynamic nomogram in an 'RStudio' panel or web browser.
}
\usage{
DynNom(model, data,
       clevel = 0.95, covariate = c("slider", "numeric"),
       ptype = c("st", "1-st"))
}
\arguments{
  \item{model}{
  an lm, glm or coxph model object
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
  \item{ptype}{
  This plot type option relates to coxph objects only. If "st" (the default) is chosen, a plot of the estimated survivor function, S(t),  is displayed.  If "1-st" is chosen a plot of 1- S(t) is displayed.
}
}
\value{
A dynamic nomogram in a shiny application
}
\references{
Banks, J. 2006. Nomograms. Encyclopedia of Statistical Sciences. 8. \cr
Easy web applications in R. \url{http://shiny.rstudio.com}
}
\author{
Amirhossein Jalali, Alberto Alvarez-Iglesias, Davood Roshan Sangachin, John Newell

Maintainer: Amirhossein Jalali <a.jalali2@nuigalway.ie>
}
\section{Please cite as:}{
Jalali, A., Alvarez-Iglesias, A., Roshan, D. Newell, J. (2015). Dynamic Nomograms for Linear, Generalized Linear and Proportional Hazard Models: The DynNom package. Version 2.0.
}
\seealso{
\code{\link{DynNom.lm}}, \code{\link{DynNom.glm}}, \code{\link{DynNom.coxph}}
}
\examples{
\dontrun{
# a simple linear regression model
model1 <- lm(uptake ~ Plant + conc + Plant * conc, data = CO2)
DynNom(model1, CO2)

# a logistic regression model
data1 =as.data.frame(Titanic)
model2 <- glm(Survived ~ Age + Class + Sex, data = data1, weights = Freq,
              family = binomial("probit"))
DynNom(model2, as.data.frame(Titanic), clevel = 0.9)

# a proportional hazard model
data.kidney <- kidney
# always make sure that the categorical variables are in a factor class
data.kidney$sex <- as.factor(data.kidney$sex)
levels(data.kidney$sex) <- c("male", "female")

model1 <- coxph(Surv(time, status)~ age + sex + disease, data.kidney)
DynNom(model1, data.kidney)
DynNom(model1, data.kidney, ptype = "1-st")
}

if (interactive()) {
# a poisson regression model
model3 <- glm(event ~ mag + station + dist + accel, data = attenu, family = poisson)
DynNom(model3, attenu, covariate = "numeric")
}
}
\keyword{ dynamic nomograms }
\keyword{ prediction }
\keyword{ shiny }
