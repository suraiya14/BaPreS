
#RF data
predict_resultsP<- function(x,y) {
  
  
library(e1071)
library(caret)

  
  
  
  
  file_n<-x
  
  file_other<-y
  
  data1 <- read.csv(file_n, header = TRUE)
  
  
  cols<-ncol(data1)
  
  data_backup<-data1
  
  nrows_training<-nrow(data1)
  
  data_other <- read.csv(file_other, header = TRUE)
  data_other_backup<-data_other
  data_other2<-data_other
  cols2<-ncol(data_other)
  
  
  
  
  
  DF2<-data_other
  data_other_backup<-data_other
  nrows_testing<-nrow(data_other)
  
  
  #print(nrows_training)
  data1<-data1[-c(cols)]
  DF1<-data1
  
  
  
  
  
  data_norm<-preProcess(data1,method=c("center", "scale"))
  #data_norm$mean
  data1<-predict(data_norm, data1)
  
  data1["Output"]<-data_backup[,cols]
  #data1$Output <- as.factor(data1$Output)
  data_other<-predict(data_norm, data_other)
 
  
  
  
  
  train<-data1
  ncol(train)
  test<-data_other
  
  suraiya<-123
  set.seed(suraiya)

  train$Output<-as.factor(train$Output)
  #tmodel2<-tune(svm, Output~., data = train, ranges = list(epsilon =seq(0,1,0.1), cost=2^(2:7)))
  tmodel2<-tune(svm, Output~., data = train, probability = TRUE , ranges = list(epsilon =seq(0,1,0.1), cost=2^(2:7)))
  #tmodel2<-tune(svm, Output~., data = train, ranges = list(epsilon =seq(0,1,0.1), cost=2^(2:7)))
  
  mymodel2<-tmodel2$best.model
  #mymodel2
  #results<-predict(mymodel2, test)
  #results
  probabilities<-predict(mymodel2, test, probability = TRUE)
  #probabilities
  
  probabilities <- attr(probabilities, "probabilities")
  #probabilities
  combined_vec <- data.frame(probabilities)
    
    

 
  
  return (combined_vec)
 
 
  
}




