logistic_plot <- function(X,Y){
  library(ggplot2)
  n=dim(data.frame(X))[2]
  if(n==1)
  {
    plot=data.frame(X=X,response=as.factor(Y))
  }else{
    plot=data.frame(X,response=as.factor(Y))
  }
  plot.list <- list()
  var.list <- colnames(plot)
  for(i in 1:n)
  {
    temp.data <- data.frame(x=plot[,i],y=as.numeric(plot$response)-1)
    gp <- ggplot(temp.data, aes(x=x, y=y)) + 
      geom_point(alpha=.5) +
      labs(x=var.list[i],y = "Response")+
      stat_smooth(method="glm", se=FALSE, method.args = list(family=binomial),col="red", lty=2)
    plot.list[i] <- list(gp)
  }
  return(plot.list)
}

Beta.init <- function(X,Y){
  Xi=X.format(X,intercept = T)
  Y=as.numeric(Y)-1
  Beta <- solve(t(Xi)%*%Xi)%*%t(Xi)%*%Y
  return(Beta)
}

P.i <- function(beta,X){
  X=X.format(X,intercept = F)
  output <- rep(NA,nrow(X))
  for(i in 1:nrow(X)){
    output[i] <- 1/(1+exp(-t(X[i,])%*%beta))
  }
  return(output)
}

loss_func <- function(beta,y,x){
  p <- P.i(beta,x)
  temp <- -y*log(p)-(1-y)*log(1-p)
  return(sum(temp))
}

Beta.hat <-  function(X,Y,method="BFGS"){
  Beta <- Beta.init(X,Y)
  Xi=X.format(X,intercept = T)
  Y=as.numeric(Y)-1
  Beta.hat <- optim(Beta,loss_func,y=Y,x=Xi,method=method)$par
  return(Beta.hat)
}



boot.confi <- function(X,Y,alpha,B=20){
  X=X.format(X,intercept = F)
  n=nrow(X)
  boot_mat <- matrix(data = NA,nrow = B,ncol = ncol(X)+1)
  colnames(boot_mat) <- c("intercept",colnames(X))
  for(i in 1:B){
    resample <- sample(1:n, replace = TRUE)
    boot_mat[i,]  <- Beta.hat(X[resample,],Y[resample])
  }
  confi.interval <- apply(boot_mat,2,quantile,probs=c(alpha/2,1-alpha/2))
  return(confi.interval)
}

logistic_pred <- function(model,X){
  Xi=X.format(X,intercept = T)
  predic <- 1/(1+exp(-Xi%*%model))
  return(predic)
}

confusion.matrix <- function(pred.value,actual.value,cutoff=0.5){
  Confusion <- data.frame(pred=pred.value,actual=actual.value)
  Confusion[Confusion>=cutoff]=1
  Confusion[Confusion<cutoff]=0
  Confusion[,c(1,2)]=lapply(Confusion[,c(1,2)], as.factor)
  N=nrow(Confusion)
  TP=nrow(Confusion[Confusion$pred=="1"&Confusion$actual=="1",])
  TN=nrow(Confusion[Confusion$pred=="0"&Confusion$actual=="0",])
  FP=nrow(Confusion[Confusion$pred=="1"&Confusion$actual=="0",])
  FN=nrow(Confusion[Confusion$pred=="0"&Confusion$actual=="1",])
  
  Accuracy <- (TP+TN)/N
  Prevalence <- (TP+FN)/N
  Sensitivity <- TP/(TP+FN)
  Specificity <- TN/(TN+FP)
  False.Discovery.Rate <- FP/(TP+FP)
  Diagnostic.Odds.Ratio <- (TP*TN)/(FN*FP)
  metrics <- c(Prevalence,Accuracy,Sensitivity,Specificity,False.Discovery.Rate,Diagnostic.Odds.Ratio)
  names(metrics) <-c("Prevalence","Accuracy","Sensitivity","Specificity","False.Discovery.Rate","Diagnostic.Odds.Ratio")
  matrix <- c(TP,FP,FN,TN)
  names(matrix) <- c("TP","FP","FN","TN")
  output.list <- list("metrics"=metrics,"matrix"=matrix)
  return(output.list)
}

metrics.plot <- function(X,Y,interval=c(0.1,0.9),step=0.1){
  library(ggplot2)
  model=Beta.hat(X,Y)
  predict <- logistic_pred(model,X)
  actual.value <- as.numeric(Y)-1
  cutoff.list <- seq(interval[1],interval[2],step)
  n=length(cutoff.list)
  data <- matrix(NA,nrow=n,ncol=6)
  for(i in 1:n){
    data[i,] <- confusion.matrix(predict,actual.value,cutoff=cutoff.list[i])$metrics
  }
  image <- data.frame(cutoff=rep(cutoff.list,6),value=c(data[,1],data[,2],data[,3],data[,4],data[,5],data[,6]),group=rep(c("Prevalence","Accuracy","Sensitivity","Specificity","False.Discovery.Rate","Diagnostic.Odds.Ratio"),n,each=n))
  plot <- ggplot(data=image,aes(x=cutoff,y=value,color=group))+geom_line(size=1)
  return(plot)
}

logistic.regression <- function(X.temp,Y.temp,method="BFGS",cutoff=0.5,alpha=0.1,B=20){
  beta.initial <- Beta.init(X.temp,Y.temp)
  model <- Beta.hat(X.temp,Y.temp,method)
  CI <- boot.confi(X.temp,Y.temp,alpha,B=20)
  predict <- logistic_pred(model,X.temp)
  actual.value <- as.numeric(Y.temp)-1
  Analysis <- confusion.matrix(predict,actual.value,cutoff=cutoff)
  
  Yi <- as.numeric(Y.temp)-1
  level=as.character(unique(Y.temp))
  names(level) <- unique(Yi)
  matrix <- matrix(Analysis$matrix,nrow=2,ncol=2)
  rownames(matrix) <- c(paste("Actual.",level["1"],sep=""),paste("Actual.",level["0"],sep=""))
  colnames(matrix) <- c(paste("Predicted.",level["1"],sep=""),paste("Predicted.",level["0"],sep=""))
  
  beta.info <- data.frame(model,beta.initial,t(CI))
  colnames(beta.info) <- c("Beta.hat","Beta.initial",paste("CI:",alpha/2,"%",sep=""),paste("CI:",1-alpha/2,"%",sep=""))
  plot <- metrics.plot(X.temp,Y.temp)
  result.list <- list("Level"=level,"Beta"=beta.info,"Confusion.Matrix"=matrix,"Metrics"=Analysis$metrics,"Plot"=plot)
  result.list
  return(result.list)
}

X.format <- function(X,intercept=F)
{
    X=as.matrix(X)
    if(dim(X)[2]==1){colnames(X) <- "Predictor"}
    if(intercept){
      intercept <- rep(1,nrow(X))
      X <- as.matrix(cbind(intercept,X))
    }
    return(X)
}



