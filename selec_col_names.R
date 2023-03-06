df<-read.csv("D:\\PARGT_Windows\\new_with_8CTD\\new_with_8CTD\\p85\\selected_train_test_merged_file.csv")
nc<-ncol(df)
df<-df[,-nc]
clnm<-colnames(df)


paste0("'",paste0(clnm, sep="", collapse="','"),"'")
