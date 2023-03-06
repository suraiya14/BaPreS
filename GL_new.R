#x<-"DLGPPISLERLDVGTNLGNAIAKLEAKELLESSD"

extractSSF <- function(x) {



#ncrna_na<-read.csv("/Applications/kayes/Antiviral\ peptide\ work/for\ X\ sequences/VEEV/VEEV_4/VEEV_4.csv",header = TRUE)

#l<-nrow(ncrna_na)
#l
ncrna_ss<-x

list.of.packages <- c("BiocManager")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
new.packages
if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')



list.of.packages <- c("Biostrings", "DECIPHER", "BiocGenerics","parallel")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) BiocManager::install(new.packages)


#source("http://bioconductor.org/biocLite.R")
#biocLite("Biostrings")

#BiocManager::install("Biostrings")
#BiocManager::install("DECIPHER")

library(Biostrings)
ncrna<-AAStringSet(ncrna_ss)
ncrna

#Secondary structure features

#BiocManager::install(c("DECIPHER", "BiocGenerics","parallel"))

#library(RSQLite)








library(DECIPHER)


#fas <- system.file("extdata", "50S_ribosomal_protein_L2.fas", package="DECIPHER")

#file = "C:\\Soil ML works\\Training\\AMR seq\\1561417752.fas.1"
#file = "C:\\Soil ML works\\Training\\PATRIC_Non-AMR seq\\1561502488.fas.1_deleted _one_4290_no_seq.1"
#file = "C:\\Soil ML works\\cleaned soil_all.fasta"

#aa = readAAStringSet(file)
#length(aa)





#aa[[1]]
#dna <- readDNAStringSet(fas)
#aa <- translate(dna)
#aa
hec <- PredictHEC(ncrna)
hec
ln<-length (hec)
#t<-hec[1]
#bb<-toString(t)
#bb[1]
# ss<-strsplit(t, "")[[1]]
# len<-length(ss)

# hec[1459]


l<-length(ncrna)

#col_n<-"ss_1"

#for(k in 2:6){
#  col_n<-paste0(col_n," ","ss_",k)
  
#}

#ls<-c(col_n)
ls<-c()

for(i in 1:l){
  t<-hec[i]
  data_psi<-strsplit(t, "")[[1]]
  len<-length(data_psi)
  add_ch<-""
  add_ch_exclude<-""
  flag<-""
  freq<-0
  x<-0
  y<-0
  z<-0
  SH<-0
  SE<-0
  SC<-0
  
  
  
  for(v in 1:len){
    if(toString(data_psi[v])=="H"){
      SH<-SH+v
      SH
      
    }
    if(toString(data_psi[v])=="E"){
      SE<-SE+v
      SE
      
    }
    if(toString(data_psi[v])=="C"){
      SC<-SC+v
      SC
      
    }
    add_ch<-paste0(add_ch,toString(data_psi[v]))
  }
  
  rr <- rle(strsplit(add_ch,"")[[1]])
  rr
  #SH<-sum(rr$lengths[which(rr$values == "H")])
  MH<-max(rr$lengths[which(rr$values == "H")])
  if (MH==-Inf || MH==Inf){
    MH<-0
  }
  
  CMVH<-SH/(len*(len-1))
  NMH<-MH/len
  
  #SE<-sum(rr$lengths[which(rr$values == "E")])
  ME<-max(rr$lengths[which(rr$values == "E")])
  ME
  if (ME==-Inf || ME==Inf){
    ME<-0
  }
  
  CMVE<-SE/(len*(len-1))
  NME<-ME/len
  
  #SC<-sum(rr$lengths[which(rr$values == "C")])
  MC<-max(rr$lengths[which(rr$values == "C")])
  
  CMVC<-SC/(len*(len-1))
  #NMC<-MC/len
  rr_len<-length(rr$values)
  rr_len
  for(v in 1:rr_len){
    if(rr$values[v]!="C")  {
      if(rr$values[v]!= flag){
        add_ch_exclude<-paste0(add_ch_exclude, toString(rr$values[v]))
        flag<-toString(rr$values[v])
        freq<-freq+1
      }
    }
  }
  add_ch_exclude
  freq
  
  #count_EHE<-str_count(add_ch_exclude, fixed("EHE"))
  count_EHE<-gregexpr("(?=EHE)",add_ch_exclude,perl=TRUE)
  count_EHE1<-count_EHE[[1]]
  count_EHE1
  if(count_EHE1[1]<0){
    count_EHE_len<-0
  }
  else{
    count_EHE_len<-length(count_EHE1)
  }
  
  count_EHE_len
  if(freq>2){
    f_EHE<-count_EHE_len/(freq-2)
  }
  else{
    f_EHE<-0
  }
  
  f_EHE
  #add_line<-paste0(toString(CMVH), " ", toString(CMVE), " ", toString(CMVC), " ", toString(NMH), " ", toString(NME), " ", toString(f_EHE))
  #add_line
  
  #ls<-c(ls,add_line)
  ls<-c(ls,CMVH, CMVE, CMVC, NMH, NME, f_EHE)
   
  
}

return (ls)
}
