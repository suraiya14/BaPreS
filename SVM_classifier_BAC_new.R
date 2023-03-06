#x<-"D:\\PARGT_Windows\\selected_train_test_merged_file.csv"
#y<-"D:\\PARGT_Windows\\input_seq.csv"


predict_results<- function(x,y) {
  
  
  #library(ROSE) 
  library(e1071)
  library(caret)
  #library(randomForest)
  #library(ROCR)
  #library(klaR)
  
  #for(kayes in 21:50){
  #kayes<-123
  
  #list_under<-c()
  #title<-paste0("aac", ",", "bl", ",", "dfr")
  
  #list_under<-c(list_under, title)
  
  
  
  
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
  
  
  
  # if(nrows_testing<100){
  #   data_other <- data_other[-((nrows_testing+1):nrows_testing_scale), ]
  # }
  # nrow(data_other)
  
  
  
  
  train<-data1
  ncol(train)
  test<-data_other
  
  kayes<-123
  set.seed(kayes)
  
  #train$Output[train$Output==1]<-"Yes"
  train$Output<-as.factor(train$Output)
  #tmodel2<-tune(svm, Output~., data = train, ranges = list(epsilon =seq(0,1,0.1), cost=2^(2:7)))
  #tmodel2<-tune(svm, Output~., data = train, ranges = list(epsilon =0, cost=8))
  tmodel2<- tune(svm, Output~., data = train, ranges = list(epsilon =seq(0,1,0.1), cost=2^(2:7)))
  tmodel2$best.parameters
  mymodel2<-tmodel2$best.model
  #mymodel2
  
  #scale(test$pseudo_10)
  #setdiff(colnames(train),colnames(test))
  #results<-predict(tmodel2, test)
  #ind <- colSums(is.na(test)) == nrow(test)
  #ind
  #pos<-which(ind==TRUE)
  #test[,pos]<-data_other_backup[,pos]
  #names(test)[ind]
  
  results<-predict(mymodel2, test)
  #print(results)
  
  #conf2_1<-confusionMatrix(predict(mymodel2, test), test$Output, positive = '1')
  
  
  #acc_under<-conf2_1$overall[1]
  
  res_vec<-c()
  for(i in 1:nrows_testing){
    if(results[i]==1){
      res_vec<-c(res_vec,1)
      
    }
    else{
      res_vec<-c(res_vec,-1)
    }
  }
  return(res_vec)
  
}

#x<-"E:\\PARGT\\training_set_aac.csv"
#y<-"E:\\PARGT\\input_seq.csv"

#predict_results(x,y)