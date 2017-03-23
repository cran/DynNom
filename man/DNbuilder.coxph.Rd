\name{DNbuilder.coxph}
\alias{DNbuilder.coxph}
\title{
Shiny code for Dynamic Nomograms
}
\description{
\code{DNbuilder.coxph} provides required scripts to deploy an \code{lm} model object as a dynamic nomogram on a server on the web such as the \url{http://shinyapps.io}.
}
\usage{
DNbuilder.coxph(model, data, clevel = 0.95, m.summary = c("raw", "formatted"),
                   covariate = c("slider", "numeric"), ptype = c("st", "1-st"))
}
\arguments{
  \item{model}{
  a coxph model object
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
  If "st" (the default) is chosen, a plot of the estimated survivor function, S(t),  is displayed.  If "1-st" is chosen a plot of 1- S(t) is displayed.
}
}
\value{
A new folder in the current working directory called DynNomapp which contains all the required scripts to deploy this dynamic nomogram on a server on the web such as the \url{http://shinyapps.io}. This folder includes \code{ui.R}, \code{server.R} and \code{global.R} script files needed to build the application and \code{dataset.rds} which is the accompanying dataset and a user guide text file called \code{README.txt} which explains how to deploy the app using all these objects.
}
\author{
Amirhossein Jalali, Alberto Alvarez-Iglesias, John Newell
}
\seealso{
\link{coxph}, \link{DynNom}, \link{DynNom.coxph}
}
\examples{
\dontrun{
data.kidney <- kidney
# always make sure that the categorical variables are in a factor class
data.kidney$sex <- as.factor(data.kidney$sex)
levels(data.kidney$sex) <- c("male", "female")

model1 <- coxph(Surv(time, status) ~ age + sex + disease, data.kidney)
DNbuilder(model1, data.kidney)
DNbuilder(model1, data.kidney, ptype = "1-st")

# a cox model including a strata term
data(lung)
model2 <- coxph(Surv(time, status) ~ age + strata(sex) + ph.ecog , data = lung)
DNbuilder(model2, lung)
}

if (interactive()) {
data.ovary <- ovarian
data.ovary$resid.ds <- as.factor(data.ovary$resid.ds)
levels(data.ovary$resid.ds) <- c("no", "yes")
data.ovary$rx <- as.factor(data.ovary$rx)
data.ovary$ecog.ps <- as.factor(data.ovary$ecog.ps)

model3 <- coxph(Surv(futime, fustat) ~ age + resid.ds * rx + ecog.ps, data = data.ovary)
DNbuilder.coxph(model3, data.ovary)
}
}
\keyword{ dynamic nomograms }
\keyword{ Cox Proportional Hazards Model }
\keyword{ individual prediction }
\keyword{ shiny }
