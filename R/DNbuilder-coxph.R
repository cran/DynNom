
DNbuilder.coxph <- function(model, data, clevel = 0.95, m.summary = c("raw", "formatted"),
                        covariate = c("slider", "numeric"), ptype = c("st", "1-st")) {

  datname <- paste(substitute(data))
  if(length(datname) > 1){
    datname <- datname[1]
    cat("\n Warning messages:
        The data frame name might be incorrect due to its complicated structure,
        You may need to edit the global.R to calling your data correctly", "\n")
  }

  covariate <- match.arg(covariate)
  m.summary <- match.arg(m.summary)
  ptype <- match.arg(ptype)

  y <- model$y
  mterms <- attr(model$terms, 'dataClasses')
  n.mterms <- names(mterms)
  mclass <- attr(model, 'class')
  mspecials <- attr(model$terms, 'specials')
  mtermslab <- attr(model$terms, 'term.labels')

  n.strata <- length(mspecials$strata)
  dim.terms <- length(n.mterms)

  model.call <- paste('Cox model:', model$call[2], sep = ' ')
  plot.title <- paste(clevel * 100, '% ', 'Confidence Interval for Survival Probability', sep = '')
  if(tail(n.mterms,n=1)=="(weights)"){
    callm = paste(paste(model$call)[1],"(",paste(model$call)[2],", ","data = data",", ","weights = ", paste(model$call)[length(paste(model$call))] ,")", sep="")
  } else{
    callm = paste(paste(model$call)[1],"(",paste(model$call)[2],", ","data = data",")", sep="")
  }

  if (ptype == 'st') {
    d.p <- "d.p <- data.frame(Prediction = exp(-pred$fit), Lower.bound = lwb,
Upper.bound = upb)"
  }
  if (ptype == '1-st') {
    d.p <- "d.p <- data.frame(Prediction = 1 - exp(-pred$fit), Lower.bound = 1 - upb,
Upper.bound = 1 - lwb)"
  }

  if (ptype == 'st') {
    pl <- "pl <- ggplot(data = dat.p()) +
geom_step(aes(x = time, y = surv, alpha = n.risk, group = part), color = coll[dat.p()$part]) +
ylim(0, 1) + xlim(0, max(dat.p()$time) * 1.05) +
labs(title = 'Estimated Survival Probability', x = 'Follow Up Time', y = 'S(t)') + theme_bw() +
theme(text = element_text(face = 'bold', size = 14), legend.position = 'none')"
}
  if (ptype == '1-st') {
    pl <- "pl <- ggplot(data = dat.p()) +
geom_step(aes(x = time, y = event, alpha = n.risk, group = part), color = coll[dat.p()$part]) +
ylim(0, 1) + xlim(0, max(dat.p()$time) * 1.05) +
labs(title = 'Estimated Probability', x = 'Follow Up Time', y = 'F(t)') +
theme_bw() + theme(text = element_text(face = 'bold', size = 14), legend.position = 'none')"
}

  if (ptype == 'st') {
    pl2 <- "pl <- ggplot(data = dat.p()) +
geom_step(aes(x = time, y = surv, group = part), color = coll[dat.p()$part]) +
ylim(0, 1) + xlim(0, max(dat.p()$time) * 1.05) +
labs(title = 'Estimated Survival Probability', x = 'Follow Up Time', y = 'S(t)') + theme_bw() +
theme(text = element_text(face = 'bold', size = 14), legend.position = 'none')"
}
  if (ptype == '1-st') {
    pl2 <- "pl <- ggplot(data = dat.p()) +
geom_step(aes(x = time, y = event, group = part), color = coll[dat.p()$part]) +
ylim(0, 1) + xlim(0, max(dat.p()$time) * 1.05) +
labs(title = 'Estimated Probability', x = 'Follow Up Time', y = 'F(t)') +
theme_bw() + theme(text = element_text(face = 'bold', size = 14), legend.position = 'none')"
}

  if (ptype == 'st') {
    ptitle <- paste("labs(title = '",paste(clevel * 100, '% ', 'Confidence Interval for Survival Probability', sep = ''),"',
x = 'Survival Probability', y = NULL)", sep="")
  }
  if (ptype == '1-st') {
    ptitle <- paste("labs(title = '",paste(clevel * 100, '% ', 'Confidence Interval for F(t)', sep = ''),"',
x = 'Probability', y = NULL)", sep="")
  }

  if (m.summary == 'raw'){
    m.print <- paste("summary(model)", sep="")
  } else{
    m.print <- paste("coef.c <- exp(model$coef)
ci.c <- exp(suppressMessages(confint(model, level = ",clevel,")))
stargazer(model, coef = list(coef.c), ci.custom = list(ci.c), p.auto = F,
type = 'text', omit.stat = c('LL', 'ser', 'f'), ci = TRUE, ci.level = ",clevel,",
single.row = TRUE, title = '",callm,"')", sep='')
  }

  sub.fit1 <- ""
  if (n.strata > 0) {
    sub.fit1 <- "sub.fit1 <- reactive({
nam <- NULL
aa <- 0
fit1 <- survfit(model, newdata = new.d())
l.s <- mspecials$strata
for (i in l.s) {
nam0 <- paste(new.d()[[which(i.factor[, 2] == i)]], sep = '')
if (aa == 0) {
nam <- paste(nam0)
}
if (aa > 0) {
nam <- paste(nam, ', ', nam0, sep = '')
}
aa <- aa + 1
}
sub.fit1 <- subset(as.data.frame(summary(fit1)[c(2:4,6:7)]), strata == nam)
return(sub.fit1)
})"
  }

  tt <- n.mterms[1]
  if (substr(tt,1,5) != "Surv("){
    stop("Error in model syntax: 'time' and 'status' do not find in model's terms")
  }

  if (length(dim(data)) > 2)
    stop("Error in data format: dataframe format required")

  if (mterms[[1]] == 'logical')
    stop("Error in model syntax: logical form for response not supported")

  if (tail(n.mterms, n = 1)=='(weights)') {
    n.terms <- length(mterms)
    mterms <- mterms[1:n.terms - 1]
  }

  if (mclass[1] == 'coxph.null') {
    stop("Error in model syntax: the model is null")
  }

  for (i in 2:dim.terms) {
    if (substr(n.mterms[i], 1, 6) == 'strata') {
      nch <- nchar(n.mterms[i])
      n.mterms[i] <- substr(n.mterms[i], 8, (nch - 1))
    }
  }

  if (!is.null(mspecials$tt)) {
    stop("Error in model syntax: coxph models with a time dependent covariate is not supported")
  }

  for(i in 2:length(n.mterms)) {
    com1=numeric(length(names(data)))
    for(j in 1:length(names(data))) {
      if (n.mterms[i]==names(data)[j]) com1[j]=1
    }
    if (sum(com1)==0)
      stop("Error in model syntax: some of model's terms do not match to variables' name in dataset")
  }

  wdir <- getwd()
  app.dir <- paste(wdir, "DynNomapp", sep="/")
  message(paste("creating new directory: ", app.dir, sep=""))
  dir.create(app.dir)
  setwd(app.dir)

  message(paste("Export dataset: ", app.dir, "/dataset.rds", sep=""))
  saveRDS(data, "dataset.rds")

  coll <- c(1:50)
  input.data <- NULL
  old.d <- NULL

  #### global.R generator
  GLOBAL=paste("library(ggplot2)
library(shiny)
library(plotly)
library(stargazer)
library(compare)
library(survival)

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
coll <- c(1:50)
mterms <- attr(model$terms, 'dataClasses')
n.mterms <- names(mterms)
mclass <- attr(model, 'class')
mspecials <- attr(model$terms, 'specials')
mtermslab <- attr(model$terms, 'term.labels')
n.strata <- length(mspecials$strata)
dim.terms <- length(n.mterms)
tt <- n.mterms[1]
for (i in 2:dim.terms) {
if (substr(n.mterms[i], 1, 6) == 'strata') {
nch <- nchar(n.mterms[i])
n.mterms[i] <- substr(n.mterms[i], 8, (nch - 1))
}
}

server = function(input, output){
observe({
if (input$quit == 1)
stopApp()
})
neededVar <- n.mterms[-1]
if (length(mtermslab) == 1) {
input.data <<- data.frame(data[1, neededVar])
names(input.data)[1] <<- n.mterms[-1]
} else {
input.data <<- data[1, neededVar]
}
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
dd <- unlist(strsplit(substr(tt, 6, nchar(tt) - 1), '[,]'))
tim <- dd[1]
sts <- substr(dd[2], 2, nchar(dd[2]))
if (length(mtermslab) == 1) {
input.data <<- data.frame(cbind(stt = NA, ti = NA, cov = NA), NO=NA)
names(input.data)[3] <<- paste(mtermslab)
names(input.data)[1:2] <<- c(paste(sts), paste(tim))
} else {
data1 <- data[, neededVar]
input.data <<- cbind(stt = NA, ti = NA, data1[1, ], NO=NA)
names(input.data)[1:2] <<- c(paste(sts), paste(tim))
input.data[1, ] <<- NA
}
if (length(i.numeric) == 0) {
i.numeric <- matrix(ncol = 3)
i.numeric <- rbind(i.numeric, V1 = paste(tim))
i.numeric[dim(i.numeric)[1], 3] <- which(names(data) == i.numeric[dim(i.numeric)[1],1])
i.numeric <- rbind(i.numeric, V1 = paste(sts))
i.numeric[dim(i.numeric)[1], 3] <- which(names(data) == i.numeric[dim(i.numeric)[1], 1])
i.numeric <- i.numeric[-1, ]
} else {
i.numeric <- rbind(i.numeric, V1 = paste(tim))
i.numeric[dim(i.numeric)[1], 3] <- which(names(data) == i.numeric[dim(i.numeric)[1], 1])
i.numeric <- rbind(i.numeric, V1 = paste(sts))
i.numeric[dim(i.numeric)[1], 3] <- which(names(data) == i.numeric[dim(i.numeric)[1], 1])
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
model$xlevels[[as.numeric(i.factor[j, 2]) - as.numeric(i.factor[j, 4])]], multiple = FALSE)
}))
do.call(tagList, slide.bars)
})
}
if (nn > 1) {
output$manySliders.n <- renderUI({
if (covariate == 'slider') {
if (nn > 2){
slide.bars <- list(lapply(1:(nn - 2), function(j) {
sliderInput(paste('numeric', j, sep = ''), i.numeric[j, 1],
min = floor(min(na.omit(data[, as.numeric(i.numeric[j, 3])]))),
max = ceiling(max(na.omit(data[, as.numeric(i.numeric[j, 3])]))),
value = mean(na.omit(data[, as.numeric(i.numeric[j, 3])])))
}), br(), checkboxInput('times', 'Predicted Survival at this Follow Up:'),
conditionalPanel(condition = 'input.times == true',
sliderInput(paste('numeric', (nn - 1), sep = ''), i.numeric[(nn - 1), 1],
min = floor(min(na.omit(data[, as.numeric(i.numeric[(nn - 1), 3])]))),
max = ceiling(max(na.omit(data[, as.numeric(i.numeric[(nn - 1), 3])]))),
value = mean(na.omit(data[, as.numeric(i.numeric[(nn - 1), 3])])))))
}
if (nn == 2){
slide.bars <- list(br(), checkboxInput('times', 'Predicted Survival at this Follow Up:'),
conditionalPanel(condition = 'input.times == true',
sliderInput(paste('numeric', (nn - 1), sep = ''), i.numeric[(nn - 1), 1],
min = floor(min(na.omit(data[, as.numeric(i.numeric[(nn - 1), 3])]))),
max = ceiling(max(na.omit(data[, as.numeric(i.numeric[(nn - 1), 3])]))),
value = mean(na.omit(data[, as.numeric(i.numeric[(nn - 1), 3])])))))
}
}
if (covariate == 'numeric') {
if (nn > 2){
slide.bars <- list(lapply(1:(nn - 2), function(j) {
numericInput(paste('numeric', j, sep = ''), i.numeric[j, 1],
value = round(mean(na.omit(data[, as.numeric(i.numeric[j, 3])]))))
}), br(), checkboxInput('times', 'Predicted Survival at this Follow Up:'),
conditionalPanel(condition = 'input.times == true',
numericInput(paste('numeric', (nn - 1), sep = ''), i.numeric[(nn - 1), 1],
value = round(mean(na.omit(data[, as.numeric(i.numeric[(nn - 1), 3])]))))))
}
if (nn == 2){
slide.bars <- list(br(), checkboxInput('times', 'Predicted Survival at this Follow Up:'),
conditionalPanel(condition = 'input.times == true',
numericInput(paste('numeric', (nn - 1), sep = ''), i.numeric[(nn - 1), 1],
value = round(mean(na.omit(data[, as.numeric(i.numeric[(nn - 1), 3])]))))))
}
}
do.call(tagList, slide.bars)
})
}
a <- 0
new.d <- reactive({
input$add
if (nf > 0) {
input.f <- vector('list', nf)
for (i in 1:nf) {
input.f[[i]] <- isolate({ input[[paste('factor', i, sep = '')]] })
names(input.f)[i] <- i.factor[i, 1]
}}
if (nn > 1) {
input.n <- vector('list', (nn - 1))
for (i in 1:(nn - 1)) {
input.n[[i]] <- isolate({ input[[paste('numeric', i, sep = '')]] })
names(input.n)[i] <- i.numeric[i, 1]
}}
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
out2 <- cbind(out[wher], NO=input$add)
input.data <<- rbind(input.data[-1], out2)
}
if (a > 0) {
wher <- match(names(out), names(input.data))
out2 <- cbind(out[wher], NO=input$add)
if (isTRUE(compare(old.d, out)) == FALSE) {
input.data <<- rbind(input.data, out2)
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
OUT <- isolate({
if (isTRUE(compare(old.d, new.d())) == FALSE) {
new.d <- cbind(stat = 1, new.d())
names(new.d)[1] <- paste(sts)
if (n.strata > 0) {
pred <- predict(model, newdata = new.d, se.fit = TRUE,
conf.int = ",clevel,", type = 'expected', reference = 'strata')
}
if (n.strata == 0) {
pred <- predict(model, newdata = new.d, se.fit = TRUE,
conf.int = ",clevel,", type = 'expected')
}
upb <- exp(-(pred$fit - (",qnorm(1 - (1 - clevel)/2)," * pred$se.fit)))
if (upb > 1) {
upb <- 1
}
lwb <- exp(-(pred$fit + (",qnorm(1 - (1 - clevel)/2)," * pred$se.fit)))
",d.p,"
old.d <<- new.d[,-1]
data.p <- cbind(d.p, counter = 1, NO = input$add)
p1 <<- rbind(p1, data.p)
p1$count <- seq(1, dim(p1)[1])
p1
} else {
p1$count <- seq(1, dim(p1)[1])
OUT <- p1
}})
}
OUT
})
s.fr <- NULL
old.d2 <- NULL
b <- 1
St <- TRUE
",sub.fit1,"
dat.p <- reactive({
if (isTRUE(compare(old.d2, new.d())) == FALSE) {
fit1 <- survfit(model, newdata = new.d())
if (n.strata == 0) {
sff <- as.data.frame(summary(fit1)[c(2:4,6:7)])
sff <- cbind(sff, event=1-sff$surv, part = b)
if (sff$time[1] != 0){
sff2 <- sff[1, ]
sff2[1, ] <- NA
sff2$time[1] <- 0
sff2$n.risk[1] <- model$n
sff2$surv[1] <- 1
sff2$event[1] <- 0
sff2$part[1] <- sff$part[1]
s.f <- rbind(sff2, sff)
} else {
s.f <- sff
}}
if (n.strata > 0) {
sff <- cbind(sub.fit1(), part = b)
sff <- cbind(sff, event=1-sff$surv)
if (sff$time[1] != 0) {
sff2 <- sff[1, ]
sff2[1, ] <- NA
sff2$time[1] <- 0
sff2$n.risk[1] <- sff[1,2]
sff2$surv[1] <- 1
sff2$event[1] <- 0
sff2$part[1] <- sff$part[1]
s.f <- rbind(sff2, sff)
} else {
s.f <- sff
}
s.f$n.risk <- s.f$n.risk/s.f$n.risk[1]
}
if (dim(s.f)[1] < 3) {
St <<- FALSE
stop('Error in data structure: There is not enough data in the current strata level')
}
s.fr <<- rbind(s.fr, s.f)
old.d2 <<- new.d()
b <<- b + 1
}
s.fr
})
output$plot <- renderPlot({
if (St == TRUE) {
if (input$add == 0)
return(NULL)
if (input$add > 0) {
if (input$trans == TRUE) {
",pl,"
}
if (input$trans == FALSE) {
",pl2,"
}}
data2()
print(pl)
}
if (St == FALSE) {
print('Restart the application')
}})
output$plot2 <- renderPlotly({
if (input$add == 0)
return(NULL)
if (is.null(new.d()))
return(NULL)
lim <- c(0, 1)
yli <- c(0 - 0.5, 10 + 0.5)
PredictNO <- 0:(sum(data2()$counter) - 1)
in.d <- data.frame(input.data[-1,-dim(input.data)[2]])
xx=matrix(paste(names(in.d), ': ',t(in.d), sep=''), ncol=dim(in.d)[1])
text.cov=apply(xx,2,paste,collapse='<br />')
if (dim(input.data)[1] > 11)
yli <- c(dim(input.data)[1] - 11.5, dim(input.data)[1] - 0.5)
p <- ggplot(data = data2(), aes(x = Prediction, y = PredictNO, text = text.cov,
label = Prediction, label2 = Lower.bound, label3=Upper.bound)) +
geom_point(size = 2, colour = data2()$count, shape = 15) +
ylim(yli[1], yli[2]) + coord_cartesian(xlim = lim) +
geom_errorbarh(xmax = data2()$Upper.bound, xmin = data2()$Lower.bound,
size = 1.45, height = 0.4, colour = data2()$count) +
",ptitle," +
theme_bw() + theme(axis.text.y = element_blank(), text = element_text(face = 'bold', size = 10))


gp=ggplotly(p, tooltip = c('text','label','label2','label3'))
dat.p()
gp
})
output$data.pred <- renderPrint({
if (input$add > 0) {
if (nrow(data2() > 0)) {
di <- ncol(input.data)
data.p <- merge(input.data[-1, ], data2()[1:5], by = 'NO')
data.p <- data.p[, !(colnames(data.p) %in% c('NO', 'counter'))]
stargazer(data.p, summary = FALSE, type = 'text')
}}
})
output$summary <- renderPrint({
",m.print,"
})
}
", sep = "")

  #### ui.R generator
  UI=paste("ui = bootstrapPage(fluidPage(
titlePanel('Dynamic Nomogram'),
sidebarLayout(sidebarPanel(uiOutput('manySliders.f'),
uiOutput('manySliders.n'),
checkboxInput('trans', 'Alpha blending (transparency)', value = TRUE),
actionButton('add', 'Predict'),
br(), br(),
helpText('Press Quit to exit the application'),
actionButton('quit', 'Quit')
),
mainPanel(tabsetPanel(id = 'tabs',
tabPanel('Estimated S(t)', plotOutput('plot')),
tabPanel('Predicted Survival', plotlyOutput('plot2')),
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
