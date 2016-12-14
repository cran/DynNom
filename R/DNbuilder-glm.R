if (getRversion() >= "2.15.1") utils::globalVariables(c("input.data", "old.d"))

DNbuilder.glm <- function(model, data,
                      clevel = 0.95, covariate = c("slider", "numeric")) {
  if (length(dim(data)) > 2)
    stop("Error in data format: dataframe format required")

  if (attr(model$terms, "dataClasses")[[1]] == "logical")
    stop("Error in model syntax: logical form for response not supported")

  if (tail(names(attr(model$terms,"dataClasses")), n = 1) == "(weights)") {
    n.terms <- length(attr(model$terms,"dataClasses"))
    attr(model$terms,"dataClasses") <- attr(model$terms,"dataClasses")[1:n.terms - 1]
  }

  for(i in 1:length(names(attr(model$terms, "dataClasses")))) {
    com1 = numeric(length(names(data)))
    for(j in 1:length(names(data))) {
      if (names(attr(model$terms, "dataClasses"))[i] == names(data)[j]) com1[j] = 1
    }
    if (sum(com1) == 0)
      stop("Error in model syntax: some of model's terms do not match to variables' name in dataset")
  }

  covariate <- match.arg(covariate)
  linkF <- model$family$linkinv
  mfamily <- model$family$family

  wdir <- getwd()
  app.dir <- paste(wdir, "DynNomapp", sep="/")
  message(paste("creating new directory: ", app.dir, sep=""))
  dir.create(app.dir)
  setwd(app.dir)

  message(paste("Export dataset: ", app.dir, "/dataset.rds", sep=""))
  saveRDS(data, "dataset.rds")

  #################################
  y <- model$model[[1]]
  mterms <- attr(model$terms, 'dataClasses')
  n.mterms <- names(mterms)
  xlevels <- model$xlevels
  df <- model$df.residual
  model.call <- paste('Linear Regression:', model$call[2], sep = ' ')
  plot.title <- paste(clevel * 100, '% ', 'Confidence Interval for Response', sep = '')
  mlfamily <- paste(model$family[[1]], "('", model$family[[2]], "')", sep="")
  if(tail(n.mterms,n=1)=="(weights)"){
    callm = paste(paste(model$call)[1],"(",paste(model$call)[2],", ","family = ",mlfamily,", ","data = data",", ","weights = ", paste(model$call)[length(paste(model$call))] ,")", sep="")
  } else{
    callm = paste(paste(model$call)[1],"(",paste(model$call)[2],", ","family=",mlfamily,", ","data = data",")", sep="")
  }

  if (mfamily == "gaussian" |
      mfamily == "inverse.gaussian" |
      mfamily == "quasi") {
    lubd <- paste("lwb <- pred$fit - (",qt(1 - (1 - clevel)/2, df)," * pred$se.fit)
   upb <- pred$fit + (",qt(1 - (1 - clevel)/2, df)," * pred$se.fit)")
  } else {
    lubd <- paste("lwb <- pred$fit - (",qnorm(1 - (1 - clevel)/2)," * pred$se.fit)
   upb <- pred$fit + (",qnorm(1 - (1 - clevel)/2)," * pred$se.fit)")
  }

  datname <- paste(substitute(data))
  if(length(datname) > 1){
    datname <- datname[1]
    cat("\n Warning messages:
        The data frame name might be incorrect due to its complicated structure.
        You need to edit the following line in global script to calling your data correctly
        data <- ....", "\n")
  }

  #### global.R generator
  GLOBAL=paste("library(ggplot2)
library(shiny)
library(stargazer)
library(compare)

##################################################################
#### You may need to edit the following lines
#### if data or the model are not defined correctly
##################################################################

data <- readRDS('dataset.rds')
model <- ",callm,"
covariate <- '", covariate,"'
", sep="")

  #### server.R generator
  SERVER=paste("input.data <- NULL
old.d <- NULL
xlevels <- model$xlevels
n.mterms <- names(attr(model$terms, 'dataClasses'))
mterms <- attr(model$terms, 'dataClasses')
linkF <- model$family$linkinv
mfamily <- model$family$family
if (tail(n.mterms, n = 1) == '(weights)') {
n.terms <- length(mterms)
mterms <- mterms[1:n.terms - 1]
n.mterms <- names(mterms)
}

server = function(input, output){
q <- observe({
if (input$quit == 1)
stopApp()
})
limits0 <- c(",mean(as.numeric(y)) - 3 * sd(y),", ",mean(as.numeric(y)) + 3 * sd(y),")
limits <- reactive({
if (as.numeric(input$limits) == 1) {
limits <- c(input$lxlim, input$uxlim)
} else {
limits <- limits0
}
})
neededVar <- n.mterms
data <- data[, neededVar]
input.data <<- data[1, ]
input.data[1, ] <<- NA
b <- 1
i.factor <- NULL
i.numeric <- NULL
for (j in 2:length(mterms)) {
for (i in 1:length(data)) {
if (n.mterms[j] == names(data)[i]) {
if (mterms[[j]] == 'factor' |
mterms[[j]] == 'ordered' |
mterms[[j]] == 'logical') {
i.factor <- rbind(i.factor, c(n.mterms[j], j, i, b))
(break)()
}
if (mterms[[j]] == 'numeric') {
i.numeric <- rbind(i.numeric, c(n.mterms[j], j, i))
b <- b + 1
(break)()
}}}
}
nn <- nrow(i.numeric)
if (is.null(nn)) {
nn <- 0
}
nf <- nrow(i.factor)
if (is.null(nf)) {
nf <- 0
}
if (nf > 0) {
output$manySliders.f <- renderUI({
slide.bars <- list(lapply(1:nf, function(j) {
selectInput(paste('factor', j, sep = ''),
names(mterms[as.numeric(i.factor[j, 2])]),
xlevels[[as.numeric(i.factor[j, 2]) - as.numeric(i.factor[j, 4])]], multiple = FALSE)
}))
do.call(tagList, slide.bars)
})
}
if (nn > 0) {
output$manySliders.n <- renderUI({
if (covariate == 'slider') {
slide.bars <- list(lapply(1:nn, function(j) {
sliderInput(paste('numeric', j, sep = ''),
names(mterms[as.numeric(i.numeric[j, 2])]),
min = as.integer(min(na.omit(data[, as.numeric(i.numeric[j, 3])]))),
max = as.integer(max(na.omit(data[, as.numeric(i.numeric[j, 3])]))) + 1,
value = as.integer(mean(na.omit(data[, as.numeric(i.numeric[j, 3])]))))
}))
}
if (covariate == 'numeric') {
slide.bars <- list(lapply(1:nn, function(j) {
numericInput(paste('numeric', j, sep = ''),
names(mterms[as.numeric(i.numeric[j, 2])]),
value = as.integer(mean(na.omit(data[, as.numeric(i.numeric[j, 3])]))))
}))
}
do.call(tagList, slide.bars)
})
}
a <- 0
new.d <- reactive({
if (nf > 0) {
input.f <- vector('list', nf)
for (i in 1:nf) {
input.f[[i]] <- local({
input[[paste('factor', i, sep = '')]]
})
names(input.f)[i] <- i.factor[i, 1]
}
}
if (nn > 0) {
input.n <- vector('list', nn)
for (i in 1:nn) {
input.n[[i]] <- local({
input[[paste('numeric', i, sep = '')]]
})
names(input.n)[i] <- i.numeric[i, 1]
}
}
if (nn == 0) {
out <- data.frame(do.call('cbind', input.f))
}
if (nf == 0) {
out <- data.frame(do.call('cbind', input.n))
}
if (nf > 0 & nn > 0) {
out <- data.frame(do.call('cbind', input.f), do.call('cbind', input.n))
}
if (a == 0) {
wher <- match(names(out), names(input.data)[-1])
out <- out[wher]
input.data <<- rbind(input.data[-1], out)
}
if (a > 0) {
wher <- match(names(out), names(input.data))
out <- out[wher]
input.data <<- rbind(input.data, out)
}
a <<- a + 1
out
})
p1 <- NULL
old.d <- NULL
data2 <- reactive({
if (input$add == 0)
return(NULL)
if (input$add > 0) {
if (isTRUE(compare(old.d, new.d())) == FALSE) {
OUT <- isolate({
pred <- predict(model, newdata = new.d(), type = 'link', conf.int = ",clevel,", se.fit = TRUE)
",lubd,"
lubound <- sort(c(linkF(lwb), linkF(upb)))
d.p <- data.frame(Prediction = linkF(pred$fit),
Lower.bound = lubound[1], Upper.bound = lubound[2])
old.d <<- new.d()
data.p <- cbind(d.p, counter = 1)
p1 <<- rbind(p1, data.p)
p1$count <- seq(1, dim(p1)[1])
p1
})
} else {
p1$count <- seq(1, dim(p1)[1])
OUT <- p1
}
}
OUT
})
if (mfamily == 'gaussian' |
mfamily == 'inverse.gaussian' |
mfamily == 'quasi') {
output$plot <- renderPlot({
if (input$add == 0)
return(NULL)
OUT <- isolate({
if (is.null(new.d()))
return(NULL)
if (is.na(input$lxlim) | is.na(input$uxlim)) {
lim <- limits0
} else {
lim <- limits()
}
yli <- c(0 - 0.5, 10 + 0.5)
if (dim(input.data)[1] > 11)
yli <- c(dim(input.data)[1] - 11.5, dim(input.data)[1] - 0.5)
p <- ggplot(data = data2(), aes(x = Prediction, y = 0:(sum(counter) - 1)))
p <- p + geom_point(size = 4, colour = data2()$count, shape = 15)
p <- p + ylim(yli[1], yli[2]) + coord_cartesian(xlim = lim)
p <- p + geom_errorbarh(xmax = data2()$Upper.bound, xmin = data2()$Lower.bound,
size = 1.45, height = 0.4, colour = data2()$count)
p <- p + labs(title = '",plot.title,"', x = 'Response', y = NULL)
p <- p + theme_bw() + theme(axis.text.y = element_blank(), text = element_text(face = 'bold', size = 14))
print(p)
})
})
}
if (mfamily == 'poisson' |
mfamily == 'quasipoisson' |
mfamily == 'Gamma') {
output$plot <- renderPlot({
if (input$add == 0)
return(NULL)
OUT <- isolate({
if (is.null(new.d()))
return(NULL)
if (is.na(input$lxlim) | is.na(input$uxlim)) {
lim <- c(0, limits0[2])
} else {
lim <- limits()
}
yli <- c(0 - 0.5, 10 + 0.5)
if (dim(input.data)[1] > 11)
yli <- c(dim(input.data)[1] - 11.5, dim(input.data)[1] - 0.5)
p <- ggplot(data = data2(), aes(x = Prediction, y = 0:(sum(counter) - 1)))
p <- p + geom_point(size = 4, colour = data2()$count, shape = 15)
p <- p + ylim(yli[1], yli[2]) + coord_cartesian(xlim = lim)
p <- p + geom_errorbarh(xmax = data2()$Upper.bound, xmin = data2()$Lower.bound,
size = 1.45, height = 0.4, colour = data2()$count)
p <- p + labs(title = '",plot.title,"', x = 'Response', y = NULL)
p <- p + theme_bw() + theme(axis.text.y = element_blank(), text = element_text(face = 'bold', size = 14))
print(p)
})
})
}
if (mfamily == 'binomial' |
mfamily == 'quasibinomial') {
output$plot <- renderPlot({
if (input$add == 0)
return(NULL)
OUT <- isolate({
if (is.null(new.d()))
return(NULL)
if (is.na(input$lxlim) | is.na(input$uxlim)) {
lim <- c(0, 1)
} else {
lim <- limits()
}
yli <- c(0 - 0.5, 10 + 0.5)
if (dim(input.data)[1] > 11)
yli <- c(dim(input.data)[1] - 11.5, dim(input.data)[1] - 0.5)
p <- ggplot(data = data2(), aes(x = Prediction, y = 0:(sum(counter) - 1)))
p <- p + geom_point(size = 4, colour = data2()$count, shape = 15)
p <- p + ylim(yli[1], yli[2]) + coord_cartesian(xlim = lim)
p <- p + geom_errorbarh(xmax = data2()$Upper.bound, xmin = data2()$Lower.bound,
size = 1.45, height = 0.4, colour = data2()$count)
p <- p + labs(title = '",plot.title,"', x = 'Response', y = NULL)
p <- p + theme_bw() + theme(axis.text.y = element_blank(), text = element_text(face = 'bold', size = 14))
print(p)
})
})
}
output$data.pred <- renderPrint({
if (input$add > 0) {
OUT <- isolate({
if (nrow(data2() > 0)) {
if (dim(input.data)[2] == 1) {
in.d <- data.frame(input.data[-1, ])
names(in.d) <- ",n.mterms[2],"
data.p <- cbind(in.d, data2()[1:3])
}
if (dim(input.data)[2] > 1) {
data.p <- cbind(input.data[-1, ], data2()[1:3])
}
stargazer(data.p, summary = FALSE, type = 'text')
}
})
}
})
output$summary <- renderPrint({
if (mfamily == 'binomial' |
mfamily == 'quasibinomial') {
coef.c <- exp(model$coef)
summ <- summary(model)
ci.c <- matrix(NA,length(model$coefficients),2)
colnames(ci.c) <- c('2.5 %','97.5 %')
rownames(ci.c) <- names(model$coefficients)

for(i in 1:length(model$coefficients)){
ci.c[i,1] <- exp(summ$coefficients[i,1] - (summ$coefficients[i,2] * ",qnorm(1 - (1 - clevel)/2),"))
ci.c[i,2] <- exp(summ$coefficients[i,1] + (summ$coefficients[i,2] * ",qnorm(1 - (1 - clevel)/2),"))
}
stargazer(model, coef = list(coef.c), ci.custom = list(ci.c), p.auto = F, type = 'text',
omit.stat = c('LL', 'ser', 'f'), ci = TRUE, ci.level = ",clevel,", single.row = TRUE,
title = '",model.call,"')
} else {
stargazer(model, type = 'text', omit.stat = c('LL', 'ser', 'f'),
ci = TRUE, ci.level = ",clevel,", single.row = TRUE,
title = '",model.call,"')
}
})}
", sep = "")

  #### ui.R generator
  UI=paste("ui = bootstrapPage(fluidPage(
titlePanel('Dynamic Nomogram'),
sidebarLayout(sidebarPanel(uiOutput('manySliders.f'),
uiOutput('manySliders.n'),
checkboxInput('limits', 'Set x-axis ranges'),
conditionalPanel(condition = 'input.limits == true',
numericInput('uxlim', 'x-axis lower', NA),
numericInput('lxlim', 'x-axis upper', NA)),
actionButton('add', 'Predict'),
br(), br(),
helpText('Press Quit to exit the application'),
actionButton('quit', 'Quit')
),
mainPanel(tabsetPanel(id = 'tabs',
tabPanel('Graphical Summary', plotOutput('plot')),
tabPanel('Numerical Summary', verbatimTextOutput('data.pred')),
tabPanel('Model Summary', verbatimTextOutput('summary'))
))))
)", sep = "")

  output=list(ui=UI, server=SERVER, global=GLOBAL)

  text <- paste("This guide will describe how to deploy a shiny application using scripts generated by DNbuilder:

1. Run the shiny app by setting your working directory to the DynNomapp folder, and then run: shiny::runApp()

If you are using the RStudio IDE, you can also run it by clicking the Run App button in the editor toolbar after open one of the R scripts.

2. You may want to modify the codes to apply all the necessary changes. Run again to confirm that your application works perfectly.

3. Deploy the application by either clicking on the Publish button in the top right corner of the running app, or use the generated files and deploy it on your server if you host any.

You can find a full guide of how to deploy an application on shinyapp.io server here:
http://docs.rstudio.com/shinyapps.io/getting-started.html#deploying-applications", sep="")
  message(paste("writing file: ", app.dir, "/README.txt", sep=""))
  writeLines(text, "README.txt")
  message(paste("writing file: ", app.dir, "/ui.R", sep=""))
  writeLines(output$ui, "ui.R")
  message(paste("writing file: ", app.dir, "/server.R", sep=""))
  writeLines(output$server, "server.R")
  message(paste("writing file: ", app.dir, "/global.R", sep=""))
  writeLines(output$global, "global.R")

  setwd(wdir)
  }
