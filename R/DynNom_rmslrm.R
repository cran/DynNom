
DynNom.lrm <- function(model, data, clevel = 0.95, m.summary = c("raw", "formatted"),
                       covariate = c("slider", "numeric")) {

  data <- data.frame(data)
  model <- update(model,x=T,y=T)

  if (length(dim(data)) > 2)
    stop("Error in data format: dataframe format required")

  if(length(class(model$y))==1){
    if (class(model$y)[1] == "logical")  stop("Error in model syntax: logical form for response not supported")} else{
      if (class(model$y)[2] == "logical")  stop("Error in model syntax: logical form for response not supported")
    }

  vars=c()
  for(i in 1:length(model$Design$name)){
    if(model$Design$assume[i]!="interaction") vars[i] <- as.character(model$Design$name[i]) else vars[i]="inter"
  }

  vars <- subset(vars,vars!="inter")
  vars <- as.character(c(model$terms[[2]],vars))

  cl.vars <-  model$Design$assume[model$Design$assume!="interaction"]
  cvars <- NULL

  if(length(class(model$y))==1){
    cvars[1] <- class(model$y)[1]} else{
      cvars[1] <- class(model$y)[2]
    }

  for(i in 2:length(vars))  cvars[i]=cl.vars[i-1]

  covariate <- match.arg(covariate)
  m.summary <- match.arg(m.summary)
  input.data <- NULL
  old.d <- NULL

  runApp(list(
    ui = bootstrapPage(fluidPage(
      titlePanel("Dynamic Nomogram"),
      sidebarLayout(sidebarPanel(uiOutput("manySliders.f"),
                                 uiOutput("manySliders.n"),
                                 checkboxInput("limits", "Set x-axis ranges"),
                                 conditionalPanel(condition = "input.limits == true",
                                                  numericInput("uxlim", "x-axis lower", NA),
                                                  numericInput("lxlim", "x-axis upper", NA)),
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
      q <- observe({
        if (input$quit == 1)
          stopApp()
      })

      limits0 <- c(0,1)
      limits <- reactive({
        if (as.numeric(input$limits) == 1) {
          limits <- c(input$lxlim, input$uxlim)
        } else {
          limits <- limits0
        }
      })

      neededVar <- vars[-1]
      data <- cbind(resp=model$y,na.omit(data[, neededVar]))
      data <- as.data.frame(data)
      colnames(data) <- c("resp",neededVar)
      input.data <<- data[1, ]
      input.data[1, ] <<- NA

      b <- 1
      i.factor <- NULL
      i.numeric <- NULL
      for (j in 2:length(vars)) {
        for (i in 1:length(data)) {
          if (vars[j] == names(data)[i]) {
            if (cvars[j] == "category" |
                cvars[j] == "scored"|
                cvars[j] == "factor"|
                cvars[j] == "ordered") {
              i.factor <- rbind(i.factor, c(vars[j], j, i, b))
              (break)()
            }
            if (cvars[j] == "rcspline"|
                cvars[j] == "asis"|
                cvars[j] == "lspline"|
                cvars[j] == "polynomial"|
                cvars[j] == "numeric" |
                cvars[j] == "integer"|
                cvars[j] == "double"|
                cvars[j] == "matrx") {
              i.numeric <- rbind(i.numeric, c(vars[j], j, i))
              b <- b + 1
              (break)()
            }
          }
        }
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
            selectInput(paste("factor", j, sep = ""),
                        vars[as.numeric(i.factor[j, 2])],
                        model$Design$parms[[i.factor[j,1]]], multiple = FALSE)
          }))
          do.call(tagList, slide.bars)
        })
      }

      if (nn > 0) {
        output$manySliders.n <- renderUI({
          if (covariate == "slider") {
            slide.bars <- list(lapply(1:nn, function(j) {
              sliderInput(paste("numeric", j, sep = ""),
                          vars[as.numeric(i.numeric[j, 2])],
                          min = floor(min(na.omit(data[, as.numeric(i.numeric[j, 3])]))),
                          max = ceiling(max(na.omit(data[, as.numeric(i.numeric[j, 3])]))),
                          value = mean(na.omit(data[, as.numeric(i.numeric[j, 3])]))
                          )
            }))
          }
          if (covariate == "numeric") {
            slide.bars <- list(lapply(1:nn, function(j) {
              numericInput(paste("numeric", j, sep = ""),
                           vars[as.numeric(i.numeric[j, 2])],
                           value = round(mean(na.omit(data[, as.numeric(i.numeric[j, 3])]))))
            }))
          }
          do.call(tagList, slide.bars)
        })
      }

      a <- 0
      new.d <- reactive({
        input$add
        if (nf > 0) {
          input.f <- vector("list", nf)
          for (i in 1:nf) {
            input.f[[i]] <- isolate({
              input[[paste("factor", i, sep = "")]]
            })
            names(input.f)[i] <- i.factor[i, 1]
          }
        }
        if (nn > 0) {
          input.n <- vector("list", nn)
          for (i in 1:nn) {
            input.n[[i]] <- isolate({
              input[[paste("numeric", i, sep = "")]]
            })
            names(input.n)[i] <- i.numeric[i, 1]
          }
        }
        if (nn == 0) {
          out <- data.frame(do.call("cbind", input.f))
        }
        if (nf == 0) {
          out <- data.frame(do.call("cbind", input.n))
        }
        if (nf > 0 & nn > 0) {
          out <- data.frame(do.call("cbind", input.f), do.call("cbind", input.n))
        }
        if (a == 0) {
          wher <- match(names(out), names(input.data)[-1])
          out <- out[wher]
          input.data <<- rbind(input.data[-1], out)
        }
        if (a > 0) {
          wher <- match(names(out), names(input.data))
          out <- out[wher]
          if (isTRUE(compare(old.d, out)) == FALSE) {
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
          if (isTRUE(compare(old.d, new.d())) == FALSE) {
            OUT <- isolate({

              if(is.error(try(predict(model, newdata = new.d(), se.fit = TRUE)))==T){
                d.p <- data.frame(Prediction = NA, Lower.bound = NA,
                                  Upper.bound = NA)
              }   else{
              pred <- predict(model, newdata = new.d(), se.fit = TRUE)
              lwb <- pred$linear.predictors - (qnorm(1 - (1 - clevel)/2) * pred$se.fit)
              upb <- pred$linear.predictors + (qnorm(1 - (1 - clevel)/2) * pred$se.fit)
              lubound <- sort(c(plogis(lwb), plogis(upb)))

              if(is.infinite(round(plogis(pred$linear.predictors),digits = 4))){
                d.p <- data.frame(Prediction = NA, Lower.bound = NA,
                                  Upper.bound = NA)
              } else{
                d.p <- data.frame(Prediction = round(plogis(pred$linear.predictors),digits = 4),
                                  Lower.bound = round(lubound[1],digits=4), Upper.bound = round(lubound[2],digits = 4))
              }
              old.d <<- new.d()
              data.p <- cbind(d.p, counter = 1)
              p1 <<- rbind(p1, data.p)
              p1$count <- seq(1, dim(p1)[1])
              p1
            }
            })
          } else {
            p1$count <- seq(1, dim(p1)[1])
            OUT <- p1
          }
        }
        OUT
      })

      output$plot <- renderPlotly({
        if (input$add == 0)
          return(NULL)
        if (is.null(new.d()))
          return(NULL)
        if (dim(na.omit(data2()))[1]==0 ) return(NULL)
        if (is.na(input$lxlim) | is.na(input$uxlim)) {
          lim <- c(0, 1)
        } else {
          lim <- limits()
        }
        PredictNO <- 0:(sum(data2()$counter) - 1)
        in.d <- data.frame(input.data[-1,])
        xx=matrix(paste(names(in.d), ": ",t(in.d), sep=""), ncol=dim(in.d)[1])
        Covariates=apply(xx,2,paste,collapse="<br />")
        yli <- c(0 - 0.5, 10 + 0.5)
        if (dim(input.data)[1] > 11)
          yli <- c(dim(input.data)[1] - 11.5, dim(input.data)[1] - 0.5)
        p <- ggplot(data = data2()[!is.na(data2()$Prediction),], aes(x = Prediction, y = PredictNO,
                                                                     text = Covariates, label = Prediction, label2 = Lower.bound, label3=Upper.bound))
        p <- p + geom_point(size = 2, colour = data2()$count[!is.na(data2()$Prediction)], shape = 15)
        p <- p + ylim(yli[1], yli[2]) + coord_cartesian(xlim = lim)
        p <- p + geom_errorbarh(xmax = data2()$Upper.bound[!is.na(data2()$Prediction)], xmin = data2()$Lower.bound[!is.na(data2()$Prediction)],
                                size = 1.45, height = 0.4, colour = data2()$count[!is.na(data2()$Prediction)])
        p <- p + labs(title = paste(clevel * 100, "% ", "Confidence Interval for Response", sep = ""),
                      x = "Probability", y = NULL)
        p <- p + theme_bw() + theme(axis.text.y = element_blank(), text = element_text(face = "bold", size = 10))
        gp=ggplotly(p, tooltip = c("text","label","label2","label3"))
        gp
      })

      output$data.pred <- renderPrint({
        if (input$add > 0) {
          if (nrow(data2() > 0)) {
            if (dim(input.data)[2] == 1) {
              in.d <- data.frame(input.data[-1, ])
              names(in.d) <- vars[2]
              data.p <- cbind(in.d, data2()[1:3])
              data.p$Prediction[is.na(data.p$Prediction)] <- "Not"
              data.p$Lower.bound[is.na(data.p$Lower.bound)] <- "IN"
              data.p$Upper.bound[is.na(data.p$Upper.bound)] <- "RANGE"
            }
            if (dim(input.data)[2] > 1) {
              data.p <- cbind(input.data[-1, ], data2()[1:3])
              data.p$Prediction[is.na(data.p$Prediction)] <- "Not"
              data.p$Lower.bound[is.na(data.p$Lower.bound)] <- "IN"
              data.p$Upper.bound[is.na(data.p$Upper.bound)] <- "RANGE"
            }
            stargazer(data.p, summary = FALSE, type = "text")
          }
        }
      })

      output$summary <- renderPrint({
        if (m.summary == "raw"){
          print(model)
        } else{
          coef.c <- exp(model$coef)
          ci.c <- matrix(NA,length(model$coefficients),2)
          colnames(ci.c) <- c("2.5 %","97.5 %")
          rownames(ci.c) <- names(model$coefficients)

          for(i in 1:length(model$coefficients)){
            ci.c[i,1] <- exp(model$coefficients[[i]] - (sqrt(vcov(model)[i,i])*qnorm(1 - (1 - clevel)/2)))
            ci.c[i,2] <- exp(model$coefficients[[i]] + (sqrt(vcov(model)[i,i])*qnorm(1 - (1 - clevel)/2)))
          }

          if(is.null(model$stat)==T){
            stargazer(model, coef = list(coef.c), ci.custom = list(ci.c), p.auto = F, type = "text",
                      omit.stat = c("LL", "ser", "f"), ci = TRUE, ci.level = clevel, single.row = TRUE,
                      title = paste("logistic", " regression (", "logit", "): ", model$call[[2]][2],
                                    " ", model$call[[2]][1], " ", model$call[[2]][3], sep = ""))
          }
          else{
            stargazer(model,model$stats, coef = list(coef.c), ci.custom = list(ci.c), p.auto = F, type = "text",
                      omit.stat = c("LL", "ser", "f"), ci = TRUE, ci.level = clevel, single.row = TRUE,
                      title = paste("logistic", " regression (", "logit", "): ", model$call[[2]][2],
                                    " ", model$call[[2]][1], " ", model$call[[2]][3], sep = ""))
          }
        }
      })}
  )
  )
}
