
utils::globalVariables(c("newcat", "."))

getdata.DN <- function(model){
  mclass <- getclass.DN(model)$model.class

  getdata.rms.DN <- function(model){
    model <- update(model, x=T, y=T)
    vars=c()
    for(i in 1:length(model$Design$name)){
      if(model$Design$assume[i]!="interaction") vars[i] <- as.character(model$Design$name[i]) else vars[i]="inter"
    }
    vars <- subset(vars,vars!="inter")
    cl.vars <-  model$Design$assume[model$Design$assume!="interaction"]
    if("category" %in% cl.vars){
      vars.cat <- vars[cl.vars=="category"]
      dummy.join.f <- function(cat.var.name){
        cat.var.levels <- model$Design$parms[names(model$Design$parms)==cat.var.name][[1]]
        cat.var.df <- as_tibble(model$x) %>%
          dplyr::select(starts_with(paste(cat.var.name,"=",sep=""))) %>%
          dplyr::select(-contains("*")) %>%
          mutate(newcat=ifelse(rowSums(.)==0,1,0)) %>%
          dplyr::select(newcat,everything())
        colnames(cat.var.df) <- cat.var.levels
        cat.var.dat <- names(cat.var.df[1:ncol(cat.var.df)])[max.col(cat.var.df[1:ncol(cat.var.df)])]
        return(cat.var.dat)
      }
      vars.cat.dat <- sapply(1:length(vars.cat),function(i) dummy.join.f(vars.cat[i]))
      colnames(vars.cat.dat) <- vars.cat
    } else{
      vars.cat.dat <- data.frame(row.names=1:nrow(model$x))
    }
    if("matrix" %in% cl.vars){
      vars.cont.mtx <- vars[cl.vars=="matrix"]
      mod.dat.cont.mtx <- data.frame(model$x) %>% dplyr::select(which(colnames(model$x)==1))
      colnames(mod.dat.cont.mtx) <- vars.cont.mtx
    } else{
      mod.dat.cont.mtx <- data.frame(row.names=1:nrow(model$x))
    }
    if(sum(!cl.vars %in% c("strata","category","matrix"))!=0){
      vars.cont.no.mtx <- vars[cl.vars!="strata" & cl.vars!="category" & cl.vars!="matrix"]
      mod.dat.cont.no.mtx <- data.frame(model$x) %>% dplyr::select(one_of(vars.cont.no.mtx))
    } else{
      mod.dat.cont.no.mtx <- data.frame(row.names=1:nrow(model$x))
    }
    mod.dat.cont <- cbind(mod.dat.cont.no.mtx,mod.dat.cont.mtx)
    if("strata" %in% cl.vars){
      vars.str <- vars[cl.vars=="strata"]
      st <- as.character(model$strata)
      st.dat <- lapply(1:length(st),function(i) unlist(strsplit(st[i], "[.]")))
      st.dat=do.call("rbind",st.dat)
      st.dat <- as.data.frame(st.dat)
      st.dat <- sapply(1:ncol(st.dat),function(x) as.factor(st.dat[,x]))
      st.dat <- sapply(1:ncol(st.dat),function(i) unlist(strsplit(st.dat[,1], "[=]"))[ c(F,T) ])
      colnames(st.dat) <- vars.str
    } else{
      st.dat <- data.frame(row.names=1:nrow(model$x))
    }
    mod.df.x <- cbind(vars.cat.dat,mod.dat.cont,st.dat)
    if(class(model)[1]=="cph"){
      model <- update(model, surv=T)

      tt <- c(model$terms[[2]])
      if(substring(tt,1,5)=="Surv("){
        dd <- unlist(strsplit(substr(tt, 6, nchar(tt) - 1), "[,]"))
        tim <- dd[1]
        sts <- substr(dd[2], 2, nchar(dd[2]))
      } else{
        tim <- colnames(model$y)[1]
        sts <- colnames(model$y)[2]
      }

      mod.df.y <- data.frame(model$y[,1],model$y[,2])
      colnames(mod.df.y) <- c(tim,sts)
    } else{
      mod.df.y <- data.frame(model$y)
      colnames(mod.df.y) <- c(model$terms[1][[2]])
    }
    mod.df <- cbind(mod.df.y,mod.df.x)

    return(mod.df)
  }

  if (!mclass %in% c("lm", "glm", "coxph", "ols", "lrm", "Glm", "cph", "gam", "Gam"))
    stop("Unrecognized model object type.")

  if (mclass %in% c("ols", "lrm", "Glm", "cph"))
    model <- update(model, x=T, y=T)

  if (any(class(try(find_data(model), silent = TRUE)) == "try-error")){
    if (mclass %in% c("ols", "lrm", "Glm", "cph")){
      dat <- getdata.rms.DN(model)
    } else{
      stop("Dataset cannot be extracted")
    }
  } else{
    dat <- find_data(model)
  }
  dat
}
