\name{DNbuilder}
\alias{DNbuilder}
\alias{DNbuilder.core}
\alias{DNbuilder.surv}
\title{
  Publishing a dynamic nomogram
}
\description{
  \code{DNbuilder} is a generic function which builds required scripts to publish a dynamic nomogram on a web server such as the \url{https://www.shinyapps.io/}. This application can then access through a URL and be used independent of R. \code{DNbuilder} supports a large number of model objects from a variety of packages.
}
\usage{
DNbuilder(model, data = NULL, clevel = 0.95, m.summary = c("raw", "formatted"),
  covariate = c("slider", "numeric"), ptype = c("st", "1-st"),
  DNtitle = NULL, DNxlab = NULL, DNylab = NULL, DNlimits = NULL,
  KMtitle = NULL, KMxlab = NULL, KMylab = NULL)

DNbuilder.core(model, data, clevel, m.summary, covariate,
  DNtitle, DNxlab, DNylab, DNlimits)

DNbuilder.surv(model, data, clevel, m.summary, covariate,
  ptype, DNtitle, DNxlab, DNylab, KMtitle, KMxlab, KMylab)
}
\arguments{
  \item{model}{
  an \code{lm}, \code{glm}, \code{coxph}, \code{ols}, \code{Glm}, \code{lrm}, \code{cph}, \code{mgcv::gam} or \code{gam::gam} model objects.
}
  \item{data}{
  a dataframe of the accompanying dataset for the model (if required).
}
  \item{clevel}{
  a confidence level for constructing the confidence interval. If not specified, a 95\% level will be used.
}
  \item{m.summary}{
  an option to choose the type of the model output represented in the 'Summary Model' tab. "raw" (the default) returns an unformatted summary of the model; "formatted" returns a formatted table of the model summary using stargazer package.
}
  \item{covariate}{
  an option to choose the type of input control widgets used for numeric values. "slider" (the default) picks out \code{sliderInput}; "numeric" picks out \code{numericInput}.
}
  \item{ptype}{
  an option for \code{coxph} or \code{cph} model objects to choose the type of plot which displays in "Survival plot" tab. "st" (the default) returns plot of estimated survivor probability (S(t)). "1-st" returns plot of estimated failure probability (1-S(t)).
}
\item{DNtitle}{
  a character vector used as the app's title. If not specified, "Dynamic Nomogram" will be used.
}
\item{DNxlab}{
  a character vector used as the title for the x-axis in "Graphical Summary" tab. If not specified, "Probability" will be used for logistic model and Cox proportional model objects; or "Response variable" for other model objects.
}
\item{DNylab}{
  a character vector used as the title for the y-axis in "Graphical Summary" tab (default is NULL).
}
\item{DNlimits}{
  a vector of 2 numeric values used to set x-axis limits in "Graphical Summary" tab. Note: This also removes the 'Set x-axis ranges' widget in the sidebar panel.
}
\item{KMtitle}{
  a character vector used as KM plot's title in "Survival plot" tab. If not specified, "Estimated Survival Probability" for \code{ptype = "st"} and "Estimated Probability" for \code{ptype = "1-st"} will be used.
}
\item{KMxlab}{
  a character vector used as the title for the x-axis in "Survival plot" tab. If not specified, "Follow Up Time" will be used.
}
\item{KMylab}{
  a character vector used as the title for the y-axis in "Survival plot" tab. If not specified, "S(t)" for \code{ptype = "st"} and "F(t)" for \code{ptype = "1-st"} will be used.
}
}
\value{
  A new folder called 'DynNomapp' will be created in the current working directory which contains all the required scripts to deploy this dynamic nomogram on a host server such as the \url{https://www.shinyapps.io/}. This folder includes \code{ui.R}, \code{server.R}, \code{global.R} and \code{data.RData} which needs to publish the app. A user guide text file (\code{README.txt}) will be also added to explain how to deploy the app using these files.
}
\references{
Banks, J. 2006. Nomograms. Encyclopedia of Statistical Sciences. 8. \cr
Easy web applications in R. \url{https://shiny.rstudio.com/} \cr
Frank E Harrell Jr (2017). rms: Regression Modeling Strategies. R package version 4.5-0.
\url{https://CRAN.R-project.org/package=rms/}
}
\author{
Amirhossein Jalali, Davood Roshan, Alberto Alvarez-Iglesias, John Newell

Maintainer: Amirhossein Jalali <a.jalali2@nuigalway.ie>
}
\section{Please cite as:}{
Jalali, A., Roshan, D., Alvarez-Iglesias, A., Newell, J. (2019). Visualising statistical models using dynamic nomograms. R package version 5.0.
}
\seealso{
\code{\link{DynNom}}, \code{\link{getpred.DN}}
}
\examples{
\dontrun{
# Simple linear regression models
fit1 <- lm(uptake ~ Plant + conc + Plant * conc, data = CO2)
DNbuilder(fit1)

t.data <- datadist(swiss)
options(datadist = 't.data')
ols(Fertility ~ Agriculture + Education + rcs(Catholic, 4), data = swiss) \%>\%
  DNbuilder(clevel = 0.9, m.summary="formatted")

# Generalized regression models
fit2 <- glm(Survived ~ Age + Class + Sex,
  data = as.data.frame(Titanic), weights = Freq, binomial("probit"))
DNbuilder(fit2, DNtitle="Titanic", DNxlab = "Probability of survival")

counts <- c(18, 17, 15, 20, 10, 20, 25, 13, 12)
outcome <- gl(3, 1, 9)
treatment <- gl(3, 3)
d <- datadist(treatment, outcome)
options(datadist = "d")
Glm((2 * counts) ~ outcome + treatment, family = poisson(),
  data = data.frame(counts, outcome, treatment)) \%>\%
  DNbuilder()

# Proportional hazard models
coxph(Surv(time, status) ~ age + strata(sex) + ph.ecog, data = lung) \%>\%
  DNbuilder()

data.kidney <- kidney
data.kidney$sex <- as.factor(data.kidney$sex)
levels(data.kidney$sex) <- c("male", "female")
coxph(Surv(time, status) ~ age + strata(sex) + disease, data.kidney) \%>\%
  DNbuilder(ptype = "1-st")

d <- datadist(veteran)
options(datadist = "d")
fit3 <- cph((Surv(time/30, status)) ~ rcs(age, 4) * strat(trt) + diagtime +
  strat(prior) + lsp(karno, 60), veteran)
DNbuilder(fit3, DNxlab = "Survival probability",
  KMtitle="Kaplan-Meier plot", KMxlab = "Time (Days)", KMylab = "Survival probability")

# Generalized additive models
mgcv::gam(Fertility ~ s(Agriculture) + Education + s(Catholic), data=swiss) \%>\%
  DNbuilder(DNlimits = c(0, 110), m.summary="formatted")

fit4 <- gam::gam(Fertility ~ Education + Catholic + s(Agriculture), fit=FALSE, data=swiss)
DNbuilder(fit4)
}
if (interactive()) {
  data(rock)
  lm(area~I(log(peri)), data = rock) \%>\%
    DNbuilder()
}
}
\keyword{ dynamic nomogram }
\keyword{ model visualisation }
\keyword{ shiny }
