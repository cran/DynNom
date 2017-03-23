\name{DynNom}
\alias{DynNom}
\title{
  Dynamic Nomograms for Linear, Generalized Linear and Proportional Hazards Models
}
\description{
  \code{DynNom} is a generic function for displaying the results of an statistical model object as a dynamic nomogram in an 'RStudio' panel or web browser. \code{DynNom} supports a variety of model objects; \code{lm}, \code{glm}, \code{coxph} and also \code{ols}, \code{Glm}, \code{lrm}, \code{cph} models in the \code{rms} package. It is a translational tool aiming to provide easy, informative individual predictions.
}
\usage{
DynNom(model, data, clevel = 0.95, m.summary = c("raw", "formatted"),
                   covariate = c("slider", "numeric"), ptype = c("st", "1-st"))
}
\arguments{
  \item{model}{
  an \code{lm}, \code{glm}, \code{coxph}, \code{ols}, \code{Glm}, \code{lrm} or \code{cph} model object
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
  \item{ptype}{
  This plot type option relates to coxph objects only. If "st" (the default) is chosen, a plot of the estimated survivor function, S(t),  is displayed.  If "1-st" is chosen a plot of 1- S(t) is displayed.
}
}
\value{
  A dynamic nomogram in a shiny application which recognises all the predictors in the model and uses them to build a sidebar panel. It sets up drop down menus for factors and sliders set at the mean and bounded by the range for covariates.

The individual predictions with a relative confidence interval are calculated using the \code{predict} function, displaying either graphically as an interactive plot in the \code{Graphical Summary} tab or a table in the \code{Numerical Summary} tab. A table of model output is also available in the \code{Model Summary} tab. In the case of the Cox proportional hazards model, estimated survivor/death function will be additionally plotted in an extra tab.
}
\references{
Banks, J. 2006. Nomograms. Encyclopedia of Statistical Sciences. 8. \cr
Easy web applications in R. \url{http://shiny.rstudio.com} \cr
Frank E Harrell Jr (2017). rms: Regression Modeling Strategies. R package version 4.5-0.
\url{https://CRAN.R-project.org/package=rms}
}
\author{
Amirhossein Jalali, Davood Roshan, Alberto Alvarez-Iglesias, John Newell

Maintainer: Amirhossein Jalali <a.jalali2@nuigalway.ie>
}
\section{Please cite as:}{
Jalali, A., Roshan, D., Alvarez-Iglesias, A., Newell, J. (2017). Dynamic Nomograms for Linear, Generalized Linear and Proportional Hazard Models. R package version 4.1.
}
\seealso{
\code{\link{DynNom.lm}}, \code{\link{DynNom.glm}}, \code{\link{DynNom.coxph}}, \code{\link{DynNom.ols}}, \code{\link{DynNom.lrm}}, \code{\link{DynNom.Glm}}, \code{\link{DynNom.cph}}
}
\examples{
\dontrun{
# simple linear regression models
model1 <- lm(uptake ~ Plant + conc + Plant * conc, data = CO2)
DynNom(model1, CO2)

data1 <- data.frame(state.x77)
model2 <- ols(Life.Exp ~ Population + Income + Illiteracy + Murder + HS.Grad +
Frost + Area,data=data1)
DynNom(model2, data1)

# Generalized regression models
data2 =as.data.frame(Titanic)
model3 <- glm(Survived ~ Age + Class + Sex, data = data2, weights = Freq,
              family = binomial("probit"))
DynNom(model3, data2, clevel = 0.9)

model4 <- lrm(formula= vs ~ wt + disp, data = mtcars)
DynNom(model4, mtcars, clevel = 0.9, m.summary = "formatted")

counts <- c(18, 17, 15, 20, 10, 20, 25, 13, 12)
outcome <- gl(3, 1, 9)
treatment <- gl(3, 3)
data2 = data.frame(counts, outcome, treatment)
model5 <- Glm((2 * counts) ~ outcome + treatment, family = poisson(), data = data2)
DynNom.Glm(model5, data2)

# a proportional hazard model
data.kidney <- kidney
# always make sure that the categorical variables are in a factor class
data.kidney$sex <- as.factor(data.kidney$sex)
levels(data.kidney$sex) <- c("male", "female")

model6 <- coxph(Surv(time, status) ~ age + sex + disease, data.kidney)
DynNom(model6, data.kidney)
DynNom(model6, data.kidney, ptype = "1-st")

model7 <-cph((Surv(log(time), status)) ~ rcs(age, 4) * strat(trt) +
              diagtime * strat(prior) + lsp(karno, 60), data = veteran)
DynNom(model7, veteran)
}

if (interactive()) {
# a poisson regression model
model8 <- glm(event ~ mag + station + dist + accel, data = attenu, family = poisson)
DynNom(model8, attenu, covariate = "numeric")
}
}
\keyword{ dynamic nomograms }
\keyword{ individual prediction }
\keyword{ shiny }
