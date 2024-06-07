
DNbuilder.core <- function(model, data, clevel, m.summary, covariate, DNtitle, DNxlab, DNylab, DNlimits) {

  mclass <- getclass.DN(model)$model.class
  mfamily <- getclass.DN(model)$model.family

  if (mclass %in% c("lm", "ols")){
    mlinkF <- function(eta) eta
  } else{
    mlinkF <- ifelse(mclass == "lrm", function(mu) plogis(mu), model$family$linkinv)
  }

  input.data <- NULL
  old.d <- NULL

  if (mclass %in% c("ols", "Glm", "lrm")){
    model <- update(model, x=T, y=T)
  }

  allvars <- all.vars(model$terms)
  if (mclass %in% c("ols", "lrm", "Glm")){
    terms <- c("-", model$Design$assume[model$Design$assume!="interaction"])
    names(terms) = allvars
    if (mclass %in% c("Glm")){
      terms <- c(attr(attr(model$model, "terms"), "dataClasses")[1], terms)
    } else{
      terms <- c(attr(model$terms,"dataClasses")[1], terms)
    }
  }
  if (mclass %in% c("lm", "glm", "gam", "Gam")){
    if(length(attr(model$terms, "dataClasses")) == length(allvars)){
      terms <- attr(model$terms, "dataClasses")
      names(terms) = allvars
    } else{
      terms <- attr(model$terms, "dataClasses")[which(names(attr(model$terms, "dataClasses")) %in% allvars)]
    }
  }

  if (terms[[1]] == "logical")
    stop("Error in model syntax: logical form for response not supported")

  terms[terms %in% c("numeric", "asis", "polynomial", "integer", "double", "matrx") |
          grepl("nmatrix", terms, fixed = T) | grepl("spline", terms, fixed = T)] = "numeric"
  terms[terms %in% c("factor", "ordered", "logical", "category", "scored")] = "factor"

  resp <- terms[1]
  names(resp) <- allvars[1]

  if ("(weights)" %in% names(terms)){
    preds0 <- as.list(terms[-c(1, length(terms))])
  } else{
    preds0 <- as.list(terms[-1])
  }
  names(preds0) <- allvars[-1]

  preds <- list()
  for (i in 1:length(preds0)){
    if (preds0[[i]] == "numeric"){
      i.dat <- which(names(preds0[i]) == names(data))
      preds[[i]] <- list(v.min = floor(min(na.omit(data[, as.numeric(i.dat)]))),
                         v.max = ceiling(max(na.omit(data[, as.numeric(i.dat)]))),
                         v.mean = zapsmall(mean(data[, as.numeric(i.dat)], na.rm=T), digits = 4)
      )
    }

    if (preds0[[i]] == "factor"){
      i.dat <- which(names(preds0[i]) == names(data))
      if (mclass %in% c("ols", "Glm", "lrm", "cph")){
        preds[[i]] <- list(v.levels = model$Design$parms[[which(names(preds0[i]) == names(model$Design$parms))]])
      } else{
        preds[[i]] <- list(v.levels = model$xlevels[[which(names(preds0[i]) == names(model$xlevels))]])
      }
    }
  }

  if (!is.null(DNlimits) & !length(DNlimits)==2)
    stop("A vector of 2 is required as 'DNlimits'")

  if (is.null(DNlimits)){
    if ((mclass %in% c("glm") & mfamily %in% c("binomial", "quasibinomial")) | mclass == "lrm"){
      limits0 <- c(0, 1)
    } else{
      if (mclass %in% c("lm", "glm", "gam", "Gam")) {
        limits0 <- c(mean(model$model[,names(resp)]) - 3 * sd(model$model[,names(resp)]),
                     mean(model$model[,names(resp)]) + 3 * sd(model$model[,names(resp)]))
      }
      if (mclass %in% c("ols", "lrm", "Glm")) {
        limits0 <- c(mean(model$y) - 3 * sd(model$y), mean(model$y) + 3 * sd(model$y))
      }
    }
    if (mclass %in% c("glm", "Glm") & mfamily %in% c("poisson", "quasipoisson", "Gamma")){
      limits0[1] <- 0
    }
  } else{
    limits0 <- DNlimits
  }

  neededVar <- c(names(resp), names(preds0))
  data <- data[, neededVar]
  input.data <- data[0, ]
  model <- update(model, data=data)

  wdir <- getwd()
  app.dir <- paste(wdir, "DynNomapp", sep="/")
  message(paste("creating new directory: ", app.dir, sep=""))
  dir.create(app.dir)
  setwd(app.dir)

  message(paste("Export dataset: ", app.dir, "/dataset.RData", sep=""))
  save(data, model, preds, preds0, resp, mlinkF, getpred.DN, getclass.DN,
       DNtitle, DNxlab, DNylab, DNlimits, limits0, terms, input.data, file = "data.RData")

  message(paste("Export functions: ", app.dir, "/functions.R", sep=""))
  dump(c("getpred.DN", "getclass.DN"), file="functions.R")

  #################################
  if (!is.null(DNlimits)) {
    limits.bl <- paste("limits <- reactive({ DNlimits })")
  } else{
    limits.bl <- paste("limits <- reactive({ if (input$limits) { limits <- c(input$lxlim, input$uxlim) } else {
                         limits <- limits0 } })")
    }

  noSE.bl <- paste("by '", mclass, "'", sep="")

  p1title.bl <- paste(clevel * 100, "% ", "Confidence Interval for Response", sep = "")
  p1msg.bl <- paste("Confidence interval is not available as there is no standard errors available by '", mclass, "' ", sep="")

  if (m.summary == "formatted"){
    if (mclass == "lm"){
      sumtitle.bl <- paste("Linear Regression:", model$call[2], sep = " ")
    }
    if (mclass %in% c("glm")){
      sumtitle.bl <- paste(mfamily, " regression (", model$family$link, "): ", model$formula[2],
                           " ", model$formula[1], " ", model$formula[3], sep = "")
    }
    if (mclass %in% c("ols", "lrm", "Glm")){
      sumtitle.bl <- paste("Linear Regression:", model$call[2], sep = " ")
    }
  } else{
    sumtitle.bl = NULL
  }

  if (m.summary == "formatted"){
    if (mclass %in% c("lm", "glm", "ols", "lrm", "Glm")){
      sum.bi <- paste("stargazer(model, type = 'text', omit.stat = c('LL', 'ser', 'f'), ci = TRUE, ci.level = clevel, single.row = TRUE, title = '",
                      sumtitle.bl,"')", sep="")
    }
    if (mclass == "gam"){
      sum.bi <- paste("Msum <- list(summary(model)$formula, summary(model)$p.table, summary(model)$s.table)
                invisible(lapply(1:3, function(i){ cat(sep='', names(Msum)[i], '\n') ; print(Msum[[i]])}))")
    }
    if (mclass == "Gam"){
      sum.bi <- paste("Msum <- list(model$formula, summary(model)$parametric.anova, summary(model)$anova)
                invisible(lapply(1:3, function(i){ cat(sep='', names(Msum)[i], '\n')) ; print(Msum[[i]])}))")
    }}
  if (m.summary == "raw"){
    if (mclass %in% c("ols", "Glm", "lrm")){
      sum.bi <- paste("print(model)")
    } else{
      sum.bi <- paste("summary(model)")
    }}

  if (mclass %in% c("ols", "Glm", "lrm")){
    datadist.bl <- paste("t.dist <- datadist(data)
options(datadist = 't.dist')", sep="")
  } else{
    datadist.bl <- ""
  }

  if (mclass %in% c("lm", "glm")){
    library.bl <- ""
  } else{
    if (mclass %in% c("ols", "Glm", "lrm")){
      library.bl <- paste("library(rms)")
    }
    if (mclass %in% c("Gam")){
      library.bl <- paste("library(gam)")
    }
    if (mclass %in% c("gam")){
      library.bl <- paste("library(mgcv)")
    }
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

### Please cite the package in publications if it is used.
# Read more on:
# Jalali, A., Alvarez-Iglesias, A., Roshan, D., & Newell, J. (2019). Visualising statistical models using dynamic nomograms. PLoS one, 14(11), e0225253.

", sep="")

  #### ui.R generator
  UI=paste("ui = bootstrapPage(fluidPage(
    titlePanel('", DNtitle,"'),
    sidebarLayout(sidebarPanel(uiOutput('manySliders'),
                               uiOutput('setlimits'),
                               actionButton('add', 'Predict'),
                               br(), br(),
                               helpText('Press Quit to exit the application'),
                               actionButton('quit', 'Quit')
    ),
    mainPanel(tabsetPanel(id = 'tabs',
                          tabPanel('Graphical Summary', plotlyOutput('plot')),
                          tabPanel('Numerical Summary', verbatimTextOutput('data.pred')),
                          tabPanel('Model Summary', verbatimTextOutput('summary'))
    )
    )
    )))", sep = "")

  #### server.R generator
  SERVER=paste('server = function(input, output){
observe({if (input$quit == 1)
          stopApp()})

', limits.bl, '

output$manySliders <- renderUI({
  slide.bars <- list()
               for (j in 1:length(preds)){
               if (terms[j+1] == "factor"){
               slide.bars[[j]] <- list(selectInput(paste("pred", j, sep = ""), names(preds0)[j], preds[[j]]$v.levels, multiple = FALSE))
               }
               if (terms[j+1] == "numeric"){
               if (covariate == "slider") {
               slide.bars[[j]] <- list(sliderInput(paste("pred", j, sep = ""), names(preds0)[j],
               min = preds[[j]]$v.min, max = preds[[j]]$v.max, value = preds[[j]]$v.mean))
               }
               if (covariate == "numeric") {
               slide.bars[[j]] <- list(numericInput(paste("pred", j, sep = ""), names(preds0)[j], value = zapsmall(preds[[j]]$v.mean, digits = 4)))
               }}}
               do.call(tagList, slide.bars)
})

output$setlimits <- renderUI({
        if (is.null(DNlimits)){
               setlim <- list(checkboxInput("limits", "Set x-axis ranges"),
               conditionalPanel(condition = "input.limits == true",
               numericInput("uxlim", "x-axis upper", zapsmall(limits0[2], digits = 2)),
               numericInput("lxlim", "x-axis lower", zapsmall(limits0[1], digits = 2))))
        } else{ setlim <- NULL }
        setlim
})

a <- 0
new.d <- reactive({
               input$add
               input.v <- vector("list", length(preds))
               for (i in 1:length(preds)) {
               input.v[[i]] <- isolate({
               input[[paste("pred", i, sep = "")]]
               })
               names(input.v)[i] <- names(preds0)[i]
               }
               out <- data.frame(lapply(input.v, cbind))
               if (a == 0) {
               input.data <<- rbind(input.data, out)
               }
               if (a > 0) {
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
               isolate({
               mpred <- suppressWarnings({ getpred.DN(model, new.d())$pred })
               se.pred <- getpred.DN(model, new.d(), set.rms=T)$SEpred
               if (is.na(se.pred)) {
               lwb <- "No standard errors"
               upb <- "', noSE.bl,'"
               pred <- mlinkF(mpred)
               d.p <- data.frame(Prediction = zapsmall(pred, digits = 3),
               Lower.bound = lwb, Upper.bound = upb)
               } else {
               lwb <- sort(mlinkF(mpred + cbind(1, -1) * (qnorm(1 - (1 - clevel)/2) * se.pred)))[1]
               upb <- sort(mlinkF(mpred + cbind(1, -1) * (qnorm(1 - (1 - clevel)/2) * se.pred)))[2]
               pred <- mlinkF(mpred)
               d.p <- data.frame(Prediction = zapsmall(pred, digits = 3),
               Lower.bound = zapsmall(lwb, digits = 3),
               Upper.bound = zapsmall(upb, digits = 3))
               }
               old.d <<- new.d()
               data.p <- cbind(d.p, counter = 1, count=0)
               p1 <<- rbind(p1, data.p)
               p1$counter <- seq(1, dim(p1)[1])
               p1$count <- 0:(dim(p1)[1]-1) %% 11 + 1
               p1
               })
               } else {
               p1$count <- seq(1, dim(p1)[1])
               }}
               rownames(p1) <- c()
               p1
})

output$plot <- renderPlotly({
  if (input$add == 0)
               return(NULL)
               if (is.null(new.d()))
               return(NULL)
               coll=c("#0E0000", "#0066CC", "#E41A1C", "#54A552", "#FF8000", "#BA55D3",
               "#006400", "#994C00", "#F781BF", "#00BFFF", "#A9A9A9")
               lim <- limits()
               yli <- c(0 - 0.5, 10 + 0.5)
               dat2 <- data2()
               if (dim(data2())[1] > 11){
               input.data = input.data[-c(1:(dim(input.data)[1]-11)),]
               dat2 <- data2()[-c(1:(dim(data2())[1]-11)),]
               yli <- c(dim(data2())[1] - 11.5, dim(data2())[1] - 0.5)
               }
               in.d <- input.data
               xx <- matrix(paste(names(in.d), ": ", t(in.d), sep = ""), ncol = dim(in.d)[1])
               Covariates <- apply(xx, 2, paste, collapse = "<br />")
               p <- ggplot(data = dat2, aes(x = Prediction, y = counter - 1, text = Covariates,
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
               gp <- ggplotly(p, tooltip = c("text", "label", "label2", "label3"))
               gp$elementId <- NULL
               gp
})

output$data.pred <- renderPrint({
  if (input$add > 0) {
               if (nrow(data2()) > 0) {
               if (dim(input.data)[2] == 1) {
               in.d <- data.frame(input.data)
               names(in.d) <- names(terms)[2]
               data.p <- cbind(in.d, data2()[1:3])
               }
               if (dim(input.data)[2] > 1) {
               data.p <- cbind(input.data, data2()[1:3])
               }}
               stargazer(data.p, summary = FALSE, type = "text")
}
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
