\name{DynNom}
\alias{DynNom}
\alias{DynNom.core}
\alias{DynNom.surv}
\title{
  Dynamic nomogram to visualise statistical models
}
\description{
  \code{DynNom} is a generic function to display the results of statistical model objects as a dynamic nomogram in an 'RStudio' panel or web browser. \code{DynNom} supports a large number of model objects from a variety of packages.
}
\usage{
DynNom(model, data = NULL, clevel = 0.95, m.summary = c("raw", "formatted"),
  covariate = c("slider", "numeric"), ptype = c("st", "1-st"),
  DNtitle = NULL, DNxlab = NULL, DNylab = NULL, DNlimits = NULL,
  KMtitle = NULL, KMxlab = NULL, KMylab = NULL)

DynNom.core(model, data, clevel, m.summary, covariate, DNtitle, DNxlab, DNylab, DNlimits)

DynNom.surv(model, data, clevel, m.summary, covariate,
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
  A dynamic nomogram in a shiny application providing individual predictions which can be used as a model visualisation or decision-making tools.

The individual predictions with a relative confidence interval are calculated using the \code{predict} function, displaying either graphically as an interactive plot in the \code{Graphical Summary} tab or a table in the \code{Numerical Summary} tab. A table of model output is also available in the \code{Model Summary} tab. In the case of the Cox proportional hazards model, an estimated survivor/failure function will be additionally displayed in a new tab.
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
Jalali, A., Roshan, D., Alvarez-Iglesias, A., Newell, J. (2019). Visualising statistical models using dynamic nomograms. R package version 5.0.
}
\seealso{
\code{\link{DNbuilder}}, \code{\link{getpred.DN}}
}
\examples{
\dontrun{
# Simple linear regression models
fit1 <- lm(uptake ~ Plant + conc + Plant * conc, data = CO2)
DynNom(fit1)

t.data <- datadist(swiss)
options(datadist = 't.data')
ols(Fertility ~ Agriculture + Education + rcs(Catholic, 4), data = swiss) \%>\%
  DynNom(clevel = 0.9, m.summary="formatted")

# Generalized regression models
fit2 <- glm(Survived ~ Age + Class + Sex,
  data = as.data.frame(Titanic), weights = Freq, family = binomial("probit"))
DynNom(fit2, DNtitle="Titanic", DNxlab = "Probability of survival")

counts <- c(18, 17, 15, 20, 10, 20, 25, 13, 12)
outcome <- gl(3, 1, 9)
treatment <- gl(3, 3)
d <- datadist(treatment, outcome)
options(datadist = "d")
Glm((2 * counts) ~ outcome + treatment, family = poisson(),
  data = data.frame(counts, outcome, treatment)) \%>\%
  DynNom()

# Proportional hazard models
coxph(Surv(time, status) ~ age + strata(sex) + ph.ecog, data = lung) \%>\%
  DynNom()

data.kidney <- kidney
data.kidney$sex <- as.factor(data.kidney$sex)
levels(data.kidney$sex) <- c("male", "female")
coxph(Surv(time, status) ~ age + strata(sex) + disease, data.kidney) \%>\%
  DynNom(ptype = "1-st")

d <- datadist(veteran)
options(datadist = "d")
fit3 <- cph((Surv(time/30, status)) ~ rcs(age, 4) * strat(trt) + diagtime +
  strat(prior) + lsp(karno, 60), veteran)
DynNom(fit3, DNxlab = "Survival probability",
       KMtitle="Kaplan-Meier plot", KMxlab = "Time (Days)", KMylab = "Survival probability")

# Generalized additive models
mgcv::gam(Fertility ~ s(Agriculture) + Education + s(Catholic), data=swiss) \%>\%
  DynNom(DNlimits = c(0, 110), m.summary="formatted")

fit4 <- gam::gam(Fertility ~ Education + Catholic + s(Agriculture), fit=FALSE, data=swiss)
DynNom(fit4)
}
if (interactive()) {
  data(rock)
  lm(area~I(log(peri)), data = rock) \%>\%
    DynNom()
}
}
\keyword{ dynamic nomogram }
\keyword{ model visualisation }
\keyword{ shiny }
