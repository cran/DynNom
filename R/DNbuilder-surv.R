
DNbuilder.surv <- function(model, data, clevel, m.summary, covariate,
                           ptype, DNtitle, DNxlab, DNylab, KMtitle, KMxlab, KMylab) {

  mclass <- getclass.DN(model)$model.class
  mlinkF <- function(mu) exp(-mu)

  input.data <- NULL
  old.d <- NULL
  Surv.in <- length(model$terms[[2]]) != 1
  if (!Surv.in)
    stop(paste("Error in model syntax: the Surv (survival object) object is created outside of", mclass))

  coll=rep(c("#0E0000", "#0066CC", "#E41A1C", "#54A552", "#FF8000", "#BA55D3",
             "#006400", "#994C00", "#F781BF", "#00BFFF", "#A9A9A9"), 10)

  if (mclass %in% c("cph")){
    model <- update(model, x=T, y=T, surv=T)
  }

  if (length(attr(model$terms, "term.labels")) == 0)
    stop("Error in model syntax: the model is null")

  if (any(!is.null(attr(model$terms, "specials")$tt)))
    stop("Error in model syntax: coxph models with a time dependent covariate is not supported")

  if (any(attr(model$terms, "dataClasses")[[1]] == "nmatrix.3", dim(model$y)[2]==3))
    stop("Error in model syntax: start/stop notation not supported")

  if (mclass %in% c("coxph")){
    strata.l <- attr(model$terms, "specials")$strata
    n.strata <- length(attr(model$terms, "specials")$strata)
    dim.terms <- length(attr(model$terms, "dataClasses"))
    if ("(weights)" %in% names(attr(model$terms, "dataClasses"))){
      attr(model$terms,"dataClasses") <- attr(model$terms, "dataClasses")[-c(length(attr(model$terms, "dataClasses")))]
    }
    mterms <- attr(model$terms, "dataClasses")[-1]
    names(mterms)=all.vars(model$terms)[-c(1,2)]
  }

  if (mclass %in% c("cph")){
    strata.l <- levels(model$strata)
    n.strata <- length(attr(model$terms, "specials")$strat)
    dim.terms <- length(model$Design$units) + 1
    mterms=model$Design$assume[model$Design$assume!="interaction"]
    names(mterms)=names(model$Design$units)
  }
  mterms[mterms %in% c("numeric", "asis", "polynomial", "integer", "double", "matrx") |
           grepl("nmatrix", mterms, fixed = T) | grepl("spline", mterms, fixed = T)] = "numeric"
  mterms[mterms %in% c("factor", "ordered", "logical", "category", "scored", "strata")] = "factor"
  preds0 <- as.list(mterms)

  tim <- all.vars(model$terms)[1:2]
  ttim <- list(v.min = floor(min(na.omit(data[tim[1]]))),
               v.max = ceiling(max(na.omit(data[tim[1]]))),
               v.mean = zapsmall(median(data[tim[1]][,1], na.rm=T), digits = 4))

  preds <- preds0
  for (i in 1:length(preds0)){
    if (preds0[[i]]  == "numeric"){
      i.dat <- which(names(preds0[i]) == names(data))
      preds[[i]] <- list(dataClasses = preds0[[i]],
                         v.min = floor(min(na.omit(data[, as.numeric(i.dat)]))),
                         v.max = ceiling(max(na.omit(data[, as.numeric(i.dat)]))),
                         v.mean = mean(data[, as.numeric(i.dat)], na.rm=T)
      )
    }
    if (preds0[[i]]  == "factor"){
      i.dat <- which(names(preds0[i]) == names(data))
      if (mclass %in% c("coxph")){
        if (startsWith(names(attr(model$terms, "dataClasses")[-1]), "strata")[i]){
          preds[[i]] <- list(dataClasses = preds0[[i]], IFstrata = T, v.levels = model$xlevels[[which(grepl(names(preds0[i]), names(model$xlevels), fixed=TRUE))]])
        } else{
          preds[[i]] <- list(dataClasses = preds0[[i]], IFstrata = F, v.levels = model$xlevels[[which(names(preds0[i]) == names(model$xlevels))]])
        }
      }
      if (mclass %in% c("cph")){
        preds[[i]] <- list(dataClasses = preds[[i]], IFstrata = model$Design$assume[i] == "strata", v.levels = model$Design$parms[[which(names(preds[i]) == names(model$Design$parms))]])
      }
    }
  }

  if (length(names(preds)) == 1) {
    input.data <- data.frame(data[0, names(preds)])
    names(input.data)[1] <- names(preds)
  } else {
    input.data <- data[0, names(preds)]
  }
  input.data <- data.frame(data.frame(tim1=0)[0,], input.data)
  names(input.data)[1] <- tim[1]

  if (mclass %in% c("coxph")){
    model <- update(model, data=data)
  }

  wdir <- getwd()
  app.dir <- paste(wdir, "DynNomapp", sep="/")
  message(paste("creating new directory: ", app.dir, sep=""))
  dir.create(app.dir)
  setwd(app.dir)

  message(paste("Export dataset: ", app.dir, "/dataset.RData", sep=""))
  save(data, model, preds, preds0, mlinkF, getpred.DN, getclass.DN, ttim, tim, ptype, n.strata, coll, dim.terms, strata.l,
       DNtitle, DNxlab, DNylab, KMtitle, KMxlab, KMylab, terms, input.data, file = "data.RData")

  message(paste("Export functions: ", app.dir, "/functions.R", sep=""))
  dump(c("getpred.DN", "getclass.DN"), file="functions.R")

  #################################

  p1title.bl <- paste(clevel * 100, "% ", "Confidence Interval for Response", sep = "")
  p1msg.bl <- paste("Confidence interval is not available as there is no standard errors available by '", mclass, "' ", sep="")

  sumtitle.bl = paste("Cox model (", model$call[1],"): ", model$call[2], sep = "")

  if (m.summary == "formatted"){
    sum.bi <- paste("coef.c <- exp(model$coef)
          ci.c <- exp(suppressMessages(confint(model, level = clevel)))
          stargazer(model, coef = list(coef.c), ci.custom = list(ci.c), p.auto = F,
                    type = 'text', omit.stat = c('LL', 'ser', 'f'), ci = TRUE, ci.level = clevel, single.row = TRUE,
                    title = '", sumtitle.bl,"')", sep="")
  }
  if (m.summary == "raw"){
    sum.bi <- paste("summary(model)")
  }

  if (mclass == "cph"){
    datadist.bl <- paste(paste(model$call[[3]]), "=data
t.dist <- datadist(", paste(model$call[[3]]),")
options(datadist = 't.dist')", sep="")
  } else{
    datadist.bl <- ""
  }

  if (mclass == "coxph"){
    library.bl <- paste("library(survival)")
  }
  if (mclass == "cph"){
    library.bl <- paste("library(rms)")
  }

  if (mclass == "cph"){
    stratlvl.bl <- paste("if (length(levels(model$strata)) != length(levels(attr(predict(model, new.d(), type='x', expand.na=FALSE), 'strata')))){
    levels(model$strata) <- levels(attr(predict(model, new.d(), type='x', expand.na=FALSE), 'strata'))
}")
  } else{
    stratlvl.bl <- ""
  }

  if (mclass %in% c("coxph")){
    nam0.bl <- paste("nam0=paste(new.d()[[names(preds[i])]], sep = '')")
  }
  if (mclass %in% c("cph")){
    nam0.bl <- paste("nam0 <- paste(names(preds[i]),'=', new.d()[[names(preds[i])]], sep = '')")
  }

  if (mclass %in% c("coxph")){
    nam.bl <- paste("nam <- paste(nam, ', ', nam0, sep = '')")
  }
  if (mclass %in% c("cph")){
    nam.bl <- paste("nam <- paste(nam, '.', nam0, sep = '')")
  }

  if (mclass %in% c("coxph")){
    stcond.bl <- paste("!try.survfit")
  }
  if (mclass %in% c("cph")){
    stcond.bl <- paste("!nam %in% strata.l")
  }

  #### global.R generator
  GLOBAL=paste("library(ggplot2)
library(shiny)
library(plotly)
library(stargazer)
library(compare)
", library.bl ,"

#######################################################
#### Before publishing your dynamic nomogram:
####
#### - You may need to edit the following lines if
#### data or model objects are not defined correctly
#### - You could modify ui.R or server.R for
#### making any required changes to your app
#######################################################

load('data.RData')
source('functions.R')
", datadist.bl,"
m.summary <- '",m.summary,"'
covariate <- '", covariate,"'
clevel <- ", clevel,"

### Please cite the package if used in publication. Use:
# Jalali A, Alvarez-Iglesias A, Roshan D, Newell J (2019) Visualising statistical models using dynamic nomograms. PLOS ONE 14(11): e0225253. https://doi.org/10.1371/journal.pone.0225253
", sep="")

  #### ui.R generator
  UI=paste("ui = bootstrapPage(fluidPage(
      titlePanel('", DNtitle,"'),
           sidebarLayout(sidebarPanel(uiOutput('manySliders'),
           checkboxInput('trans', 'Alpha blending (transparency)', value = TRUE),
           actionButton('add', 'Predict'),
           br(), br(),
           helpText('Press Quit to exit the application'),
           actionButton('quit', 'Quit')
           ),
           mainPanel(tabsetPanel(id = 'tabs',
           tabPanel('Survival plot', plotOutput('plot')),
           tabPanel('Predicted Survival', plotlyOutput('plot2')),
           tabPanel('Numerical Summary', verbatimTextOutput('data.pred')),
           tabPanel('Model Summary', verbatimTextOutput('summary'))
           )
           )
           )))", sep = "")


  #### server.R generator
  SERVER=paste('server = function(input, output){
observe({if (input$quit == 1)
          stopApp()})

output$manySliders <- renderUI({
        slide.bars <- list()
               for (j in 1:length(preds)){
               if (preds[[j]]$dataClasses == "factor"){
               slide.bars[[j]] <- list(selectInput(names(preds)[j], names(preds)[j], preds[[j]]$v.levels, multiple = FALSE))
               }
               if (preds[[j]]$dataClasses == "numeric"){
               if (covariate == "slider") {
               slide.bars[[j]] <- list(sliderInput(names(preds)[j], names(preds)[j],
               min = preds[[j]]$v.min, max = preds[[j]]$v.max, value = preds[[j]]$v.mean))
               }
               if (covariate == "numeric") {
               slide.bars[[j]] <- list(numericInput(names(preds)[j], names(preds)[j], value = zapsmall(preds[[j]]$v.mean, digits = 4)))
               }}}
               if (covariate == "slider") {
               slide.bars[[length(preds) + 1]] <-
               list(br(), checkboxInput("times", "Predicted Survival at this Follow Up:"),
               conditionalPanel(condition = "input.times == true",
               sliderInput("tim", tim[1], min = ttim$v.min, max = ttim$v.max, value = ttim$v.mean)))
               } else {
               slide.bars[[length(preds) + 1]] <-
               list(br(), checkboxInput("times", "Predicted Survival at this Follow Up:"),
               conditionalPanel(condition = "input.times == true",
               numericInput("tim", tim[1], value = zapsmall(ttim$v.mean, digits = 4))))
               }
               do.call(tagList, slide.bars)
})

a <- 0
      old.d <- NULL
               new.d <- reactive({
               input$add
               input.v <- vector("list", length(preds) + 1)
               input.v[[1]] <- isolate({ input[["tim"]] })
               names(input.v)[1] <- tim[1]
               for (i in 1:length(preds)) {
               input.v[[i+1]] <- isolate({
               input[[names(preds)[i]]]
               })
               names(input.v)[i+1] <- names(preds)[i]
               }
               out <- data.frame(lapply(input.v, cbind))
               if (a == 0) {
               wher <- match(names(out), names(input.data))
               out <- out[wher]
               input.data <<- rbind(input.data, out)
               }
               if (a > 0) {
               wher <- match(names(out), names(input.data))
               out <- out[wher]
               if (!isTRUE(compare(old.d, out))) {
               input.data <<- rbind(input.data, out)
               }}
               a <<- a + 1
               out
               })

p1 <- NULL
old.d <- NULL
data2 <- reactive({
               if (input$add == 0)
               return(NULL)
               if (input$add > 0) {
               if (!isTRUE(compare(old.d, new.d()))) {
               OUT <- isolate({
               new.d <- cbind(st.ind = 1, new.d())
               names(new.d)[1] <- tim[2]
               DNpred <- getpred.DN(model, new.d)
               mpred <- DNpred$pred
               se.pred <- DNpred$SEpred
               pred <- mlinkF(mpred)
               if (is.na(se.pred)) {
               lwb <- NULL
               upb <- NULL
               } else {
               lwb <- sort(mlinkF(mpred + cbind(1, -1) * (qnorm(1 - (1 - clevel)/2) * se.pred)))[1]
               upb <- sort(mlinkF(mpred + cbind(1, -1) * (qnorm(1 - (1 - clevel)/2) * se.pred)))[2]
               if (upb > 1) {
               upb <- 1
               }}
               if (ptype == "st") {
                d.p <- data.frame(Prediction = zapsmall(pred, digits = 2),
               Lower.bound = zapsmall(lwb, digits = 2),
               Upper.bound = zapsmall(upb, digits = 2))
               }
               if (ptype == "1-st") {
               d.p <- data.frame(Prediction = zapsmall(1-pred, digits = 2),
               Lower.bound = zapsmall(1-upb, digits = 2),
               Upper.bound = zapsmall(1-lwb, digits = 2))
               }
               old.d <<- new.d[,-1]
               data.p <- cbind(d.p, counter = TRUE)
               if (DNpred$InRange){
               p1 <<- rbind(p1[,-5], data.p)
               } else{
               p1 <<- rbind(p1[,-5], data.frame(Prediction = NA, Lower.bound = NA, Upper.bound = NA, counter = FALSE))
               }
               p1
               })
               } else {
               p1$count <- seq(1, dim(p1)[1])
               }}
               p1
})

s.fr <- NULL
old.d2 <- NULL
b <- 1
dat.p <- reactive({
               if (!isTRUE(compare(old.d2, new.d()))) {
               ', stratlvl.bl,'
               try.survfit <- !any(class(try(survfit(model, newdata = new.d()), silent = TRUE)) == "try-error")
               if (try.survfit){
               fit1 <- survfit(model, newdata = new.d())
               }
               if (n.strata == 0) {
               sff <- data.frame(summary(fit1)[c("time", "n.risk", "surv")])
               sff <- cbind(sff, event=1-sff$surv, part = b)
               if (sff$time[1] != 0){
               sff <- rbind(data.frame(time=0, n.risk=sff$n.risk[1] ,surv=1, event=0, part=sff$part[1]), sff)
               }}
               if (n.strata > 0) {
               nam <- NULL
               new.sub <- T
               for (i in 1:(dim.terms-1)) {
               if (preds[[i]]$dataClasses == "factor"){
               if (preds[[i]]$IFstrata){
               ', nam0.bl,'
               if (new.sub) {
               nam <- paste(nam0)
               new.sub <- F
               } else {
               ', nam.bl,'
               }}}}
               if (try.survfit){
               survfit_df <- as.data.frame(summary(fit1)[c("time", "n.risk", "strata", "surv")])
               survfit_df$strata <- as.factor(survfit_df$strata)
               levels(survfit_df$strata) <-  nam
               sub.fit1 <- survfit_df
               } else{
               sub.fit1 <- data.frame(time=NA, n.risk=NA, strata=NA, surv=NA, event=NA, part=NA)[0,]
               }
               if (', stcond.bl,'){
               message("The strata levels not found in the original")
               sff <- cbind(sub.fit1, event=NULL, part = NULL)
               b <<- b - 1
               } else{
               sff <- cbind(sub.fit1, event=1-sub.fit1$surv, part = b)
               if (sff$time[1] != 0) {
               sff <- rbind(data.frame(time=0, n.risk=sff$n.risk[1], strata=sff$strata[1] ,surv=1, event=0, part=sff$part[1]), sff)
               }
               sff$n.risk <- sff$n.risk/sff$n.risk[1]
               }
               sff$n.risk <- sff$n.risk/sff$n.risk[1]
               }
               s.fr <<- rbind(s.fr, sff)
               old.d2 <<- new.d()
               b <<- b + 1
               }
               s.fr
})

dat.f <- reactive({
        if (nrow(data2() > 0))
          cbind(input.data, data2()[1:3])
})

# KM plot
output$plot <- renderPlot({
               data2()
               if (input$add == 0)
               return(NULL)
               if (input$add > 0) {
               if (ptype == "st") {
               if (input$trans == TRUE) {
               pl <- ggplot(data = dat.p()) +
               geom_step(aes(x = time, y = surv, alpha = n.risk, group = part), color = coll[dat.p()$part])
               }
               if (input$trans == FALSE) {
               pl <- ggplot(data = dat.p()) +
               geom_step(aes(x = time, y = surv, group = part), color = coll[dat.p()$part])
               }}
               if (ptype == "1-st") {
               if (input$trans == TRUE) {
               pl <- ggplot(data = dat.p()) +
               geom_step(aes(x = time, y = event, alpha = n.risk, group = part), color = coll[dat.p()$part])
               }
               if (input$trans == FALSE) {
               pl <- ggplot(data = dat.p()) +
               geom_step(aes(x = time, y = event, group = part), color = coll[dat.p()$part])
               }}
               pl <- pl + ylim(0, 1) + xlim(0, max(dat.p()$time) * 1.05) +
               labs(title = "', KMtitle,'", x = "', KMxlab,'", y = "', KMylab,'") + theme_bw() +
               theme(text = element_text(face = "bold", size = 12), legend.position = "none", plot.title = element_text(hjust = .5))
               }
               print(pl)
})

output$plot2 <- renderPlotly({
        if (input$add == 0)
               return(NULL)
               if (is.null(new.d()))
               return(NULL)
               lim <- c(0, 1)
               yli <- c(0 - 0.5, 10 + 0.5)
               input.data = input.data[data2()$counter,]
               in.d <- data.frame(input.data)
               xx=matrix(paste(names(in.d), ": ",t(in.d), sep=""), ncol=dim(in.d)[1])
               text.cov=apply(xx,2,paste,collapse="<br />")
               if (dim(input.data)[1] > 11)
               yli <- c(dim(input.data)[1] - 11.5, dim(input.data)[1] - 0.5)
               dat2 <- data2()[data2()$counter,]
               dat2$count = seq(1, nrow(dat2))

               p <- ggplot(data = dat2, aes(x = Prediction, y = count - 1, text = text.cov,
               label = Prediction, label2 = Lower.bound, label3=Upper.bound)) +
               geom_point(size = 2, colour = coll[dat2$count], shape = 15) +
               ylim(yli[1], yli[2]) + coord_cartesian(xlim = lim) +
               labs(title = "', p1title.bl,'",
               x = "', DNxlab,'", y = "', DNylab,'") + theme_bw() +
               theme(axis.text.y = element_blank(), text = element_text(face = "bold", size = 10))
               if (is.numeric(dat2$Upper.bound)){
               p <- p + geom_errorbarh(xmax = dat2$Upper.bound, xmin = dat2$Lower.bound,
               size = 1.45, height = 0.4, colour = coll[dat2$count])
               } else{
               message("', p1msg.bl,'")
               }
               if (ptype == "st") {
               p <- p + labs(title = paste(clevel * 100, "% ", "Confidence Interval for Survival Probability", sep = ""),
               x = DNxlab, y = DNylab)
               }
               if (ptype == "1-st") {
               p <- p + labs(title = paste(clevel * 100, "% ", "Confidence Interval for F(t)", sep = ""),
               x = DNxlab, y = DNylab)
               }
               gp=ggplotly(p, tooltip = c("text","label","label2","label3"))
               gp$elementId <- NULL
               dat.p()
               gp
})

output$data.pred <- renderPrint({
        if (input$add > 0) {
               if (nrow(data2() > 0)) {
               stargazer(dat.f(), summary = FALSE, type = "text")
        }}
})


output$summary <- renderPrint({
', sum.bi,'
})
}', sep = "")

  #################################

  output=list(ui=UI, server=SERVER, global=GLOBAL)

  text <- paste("This guide will describe how to deploy a shiny application using scripts generated by DNbuilder:

1. Run the shiny app by setting your working directory to the DynNomapp folder, and then run: shiny::runApp() If you are using the RStudio IDE, you can also run it by clicking the Run App button in the editor toolbar after open one of the R scripts.

2. You could modify codes to apply all the necessary changes. Run again to confirm that your application works perfectly.

3. Deploy the application by either clicking on the Publish button in the top right corner of the running app, or use the generated files and deploy it on your server if you host any.

You can find a full guide of how to deploy an application on shinyapp.io server here:
http://docs.rstudio.com/shinyapps.io/getting-started.html#deploying-applications

Please cite the package if using in publication.", sep="")
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
