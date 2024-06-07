
utils::globalVariables(c("Prediction", "counter", "Lower.bound", "Upper.bound"))

DynNom.core <- function(model, data, clevel, m.summary, covariate, DNtitle, DNxlab, DNylab, DNlimits) {

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
  }

  neededVar <- c(names(resp), names(preds0))
  data <- data[, neededVar]
  input.data <- data[0, ]

  runApp(list(
    ui = bootstrapPage(fluidPage(
      titlePanel(DNtitle),
      sidebarLayout(sidebarPanel(uiOutput("manySliders"),
                                 uiOutput("setlimits"),
                                 actionButton("add", "Predict"),
                                 br(), br(),
                                 helpText("Press Quit to exit the application"),
                                 actionButton("quit", "Quit")
      ),
      mainPanel(tabsetPanel(id = "tabs",
                            tabPanel("Graphical Summary", plotlyOutput("plot")),
                            tabPanel("Numerical Summary", verbatimTextOutput("data.pred")),
                            tabPanel("Model Summary", verbatimTextOutput("summary"))
      )
      )
      ))),

    server = function(input, output){

      observe({
        if (input$quit == 1)
          stopApp()
      })

      limits <- reactive({
        if (!is.null(DNlimits)) {
          limits <- DNlimits
        } else{
          if (input$limits) {
            limits <- c(input$lxlim, input$uxlim)
          } else {
            limits <- limits0
          }
        }
      })

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
            }
          }
        }
        do.call(tagList, slide.bars)
      })

      output$setlimits <- renderUI({
        if (is.null(DNlimits)){
          setlim <- list(checkboxInput("limits", "Set x-axis ranges"),
                         conditionalPanel(condition = "input.limits == true",
                                          numericInput("uxlim", "x-axis upper", zapsmall(limits0[2], digits = 2)),
                                          numericInput("lxlim", "x-axis lower", zapsmall(limits0[1], digits = 2))))
        } else{
          setlim <- NULL
        }
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
          }
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
          if (!isTRUE(compare(old.d, new.d()))) {
            isolate({
              mpred <- suppressWarnings({ getpred.DN(model, new.d())$pred })
              se.pred <- getpred.DN(model, new.d())$SEpred
              if (is.na(se.pred)) {
                lwb <- "No standard errors"
                upb <- paste("by '", mclass, "'", sep="")
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
          }
        }
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
          labs(title = paste(clevel * 100, "% ", "Confidence Interval for Response", sep = ""),
               x = DNxlab, y = DNylab) + theme_bw() +
          theme(axis.text.y = element_blank(), text = element_text(face = "bold", size = 10))
        if (is.numeric(dat2$Upper.bound)){
          p <- p + geom_errorbarh(xmax = dat2$Upper.bound, xmin = dat2$Lower.bound,
                                  size = 1.45, height = 0.4, colour = coll[dat2$count])

        } else{
          message(paste("Confidence interval is not available as there is no standard errors available by '", mclass, "' ", sep=""))
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
            }
          }
          stargazer(data.p, summary = FALSE, type = "text")
        }
      })

      output$summary <- renderPrint({
        if (m.summary == "formatted"){
          if (mclass == "lm"){
            stargazer(model, type = "text", omit.stat = c("LL", "ser", "f"), ci = TRUE, ci.level = clevel,
                      single.row = TRUE, title = paste("Linear Regression:", model$call[2], sep = " "))
          }

          if (mclass %in% c("glm")){
            stargazer(model, type = "text", omit.stat = c("LL", "ser", "f"),
                      ci = TRUE, ci.level = clevel, single.row = TRUE,
                      title = paste(mfamily, " regression (", model$family$link, "): ", model$formula[2],
                                    " ", model$formula[1], " ", model$formula[3], sep = ""))
          }

          if (mclass == "gam"){
            Msum <- list(summary(model)$formula, summary(model)$p.table, summary(model)$s.table)
            invisible(lapply(1:3, function(i){ cat(sep="", names(Msum)[i], "\n") ; print(Msum[[i]])}))
          }

          if (mclass == "Gam"){
            Msum <- list(model$formula, summary(model)$parametric.anova, summary(model)$anova)
            invisible(lapply(1:3, function(i){ cat(sep="", names(Msum)[i], "\n") ; print(Msum[[i]])}))
          }

          if (mclass %in% c("ols", "lrm", "Glm")){
            stargazer(model, type = "text", omit.stat = c("LL", "ser", "f"), ci = TRUE, ci.level = clevel,
                      single.row = TRUE, title = paste("Linear Regression:", model$call[2], sep = " "))
          }
        }

        if (m.summary == "raw"){
          if (mclass %in% c("ols", "Glm", "lrm")){
            print(model)
          } else{
            summary(model)
          }
        }
      })
    }
  )
  )
}
