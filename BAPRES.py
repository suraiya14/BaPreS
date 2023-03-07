#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#import time
#start_time = time.time()
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
import os
import sys
import subprocess
import rpy2
from rpy2 import *
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector
#from rpy2.robjects import FloatVector, r
#from rpy2.robjects import globalenv
from rpy2.robjects import numpy2ri
#from rpy2.robjects.numpy2ri import numpy2ri
from rpy2.robjects.packages import STAP
#numpy2ri.activate()
#from rpy2.robjects import pandas2ri
import shutil
from shutil import copyfile
import math
import re
import fileinput
import tkinter
from tkinter import *
#import tkFont
# import filedialog module
from tkinter import filedialog
#from PIL import Image
#from PIL import ImageTk
from IPython.display import display
from IPython.display import HTML
import IPython.core.display as di 

di.display_html('<script>jQuery(function() {if (jQuery("body.notebook_app").length == 0) { jQuery(".input_area").toggle(); jQuery(".prompt").toggle();}});</script>', raw=True)

di.display_html('''<button onclick="jQuery('.input_area').toggle(); jQuery('.prompt').toggle();">Toggle code</button>''', raw=True)


#generating 621D features

def feature_extraction(file_name,svalue):
    
    fasta_file_name=file_name + ".fasta"
    lines = [line.rstrip('\n') for line in open(fasta_file_name)]    
    save_sequences= ''
    sequences= []
    sequences_name= []
    for l in lines:
        if(l == ""): #blank lines are disregarded
                pass
        elif (l[0] == '>'):
            sequences_name.append(l[1:])
            sequences.append(save_sequences)
            save_sequences= ''
        else:
            save_sequences+= l

    sequences.append(save_sequences)
    del sequences[0]
    
    
    all_generated_features=pd.DataFrame()
    

     
    #amino acid composition
    
    def aa_composition(seq):
        
        r=ro.r
        r.source("AAC_BAC.R") #calling R script to compute amino acid composition
  
        aa_comp = r.extractAAC_BAC(seq)
        return aa_comp
    
    #dipeptide composition
    
    def dc_composition(seq):
        
        r=ro.r
        r.source("DC_BAC.R") #calling R script to compute amino acid composition
  
        dc_comp = r.extractDC_BAC(seq)
        return dc_comp
    
    #pseudo-amino acid composition
    
    def pseaac_composition(seq,path):
        
        r=ro.r
        r.source("file_with_pkgTest.R") #check whether required R package is installed
        r.pkgTest("base")
        r.pkgTest("protr")
    #r.source("SVM_prediction.R") #calling R script to get prediction results from SVM
    #r.source("file_with_pkgTest_decipher.R")
    #r.pkgTest_decipher("DECIPHER")
        r.source("PAAC_BAC.R") #calling R script to compute amino acid composition
  
        pseaac_comp = r.extractPAAC_BAC(seq,path)
        return pseaac_comp
    
    #amphiphilic pseudo-amino acid composition
    
    def apseaac_composition(seq,path):
        
        r=ro.r
        r.source("file_with_pkgTest.R") #check whether required R package is installed
        r.pkgTest("base")
        r.pkgTest("protr")
    #r.source("SVM_prediction.R") #calling R script to get prediction results from SVM
    #r.source("file_with_pkgTest_decipher.R")
    #r.pkgTest_decipher("DECIPHER")
        r.source("APAAC_BAC.R") #calling R script to compute amino acid composition
  
        apseaac_comp = r.extractAPAAC_BAC(seq,path)
        return apseaac_comp
    
    #composition for CTD model 
    
    def ctd_composition(seq):
        
        r=ro.r
        r.source("CTDC.R") #calling R script to compute composition for CTD model
  
        ctd_comp = r.extractCTDC_BAC(seq)
        return ctd_comp
    
    
    #transition for CTD model
    
    def ctd_transition(seq):
        
        r=ro.r
        r.source("CTDT.R") #calling R script to compute transition for CTD model
  
        ctd_trans = r.extractCTDT_BAC(seq)
        return ctd_trans
    
    
    #distribution for CTD model
    
    def ctd_distribution(seq):
        
        r=ro.r
        r.source("CTDD.R") #calling R script to compute distribution for CTD model
  
        ctd_distr = r.extractCTDD_BAC(seq)
        return ctd_distr
    

    
    #generate secondary structure features
        
    def secondstruct_feat(seq):
        
        
        base = importr('base')
        
        r=ro.r
        #r.source("file_with_pkgTest.R") #check whether required R package is installed
        #r.pkgTest("base")
        
        #r.source("file_with_pkgTest_decipher.R")
        #r.pkgTest_decipher("DECIPHER")
        r.source("GL_new.R") #calling R script to compute secondary structure features
  
        ss_ft = r.extractSSF(seq)
        #print(ssfile)
        return ss_ft
    
    
    #####################################################################################
    
    
    aa_dict='ARNDCEQGHILKMFPSTWYV'
    
    for i in range(1,21):
        all_generated_features["aac_%s"%i]=""
    
    for i in range(1,401):
        all_generated_features["dipep_%s"%i]=""
        
    for i in range(1,31):
        all_generated_features["pseudo_%s"%i]=""
        
    for i in range(1,41):
        all_generated_features["amphipseudo_%s"%i]=""
        
    for i in range(1,22):
        all_generated_features["comp_%s"%i]=""
    
    for i in range(1,22):
        all_generated_features["tran_%s"%i]=""
        
    for i in range(1,106):
        all_generated_features["dist_%s"%i]=""
    
    for i in range(1,7):
        all_generated_features["ss_%s"%i]=""
    
    
    # Execute for each sequence
    
    seqs_length=[]
    counter=0
    counting_seq=0

    for seq in sequences:
        seqs_length.append(len(seq))
        counter+=1
        aa_comp=aa_composition(seq)
        dir_path = os.getcwd()
        dir_path = dir_path.replace("\\", "/")
        dir_path+="/AAidx.csv"
        dc_comp=dc_composition(seq)
        pseaac_comp=pseaac_composition(seq,dir_path)
        apseaac_comp=apseaac_composition(seq,dir_path)
        ctd_comp = ctd_composition(seq)
        ctd_trans = ctd_transition(seq)
        ctd_distr = ctd_distribution(seq)
        ss_struct = secondstruct_feat(seq)
        
        all_generated_features.at[counting_seq,"aac_1":"aac_20"]=aa_comp
        all_generated_features.at[counting_seq,"dipep_1":"dipep_400"]=dc_comp
        all_generated_features.at[counting_seq,"pseudo_1":"pseudo_30"]=pseaac_comp
        all_generated_features.at[counting_seq,"amphipseudo_1":"amphipseudo_40"]=apseaac_comp
        all_generated_features.at[counting_seq,"comp_1":"comp_21"]=ctd_comp
        all_generated_features.at[counting_seq,"tran_1":"tran_21"]=ctd_trans
        all_generated_features.at[counting_seq,"dist_1":"dist_105"]=ctd_distr[:105]
        all_generated_features.at[counting_seq,"ss_1":"ss_6"]=ss_struct
        
        x_len=len(seq)
        
        #get path for current directiry
        

        
        
        counting_seq+=1
    
    
    
    eliminate= []
    
    if svalue==0:
        #features_selected=["S5", "D80", "C21", "S4", "T10", "D81"]
        features_selected=["pseudo_11","pseudo_2","dist_3","dist_76","dist_18","dist_62","dist_93","pseudo_13","dist_75","pseudo_14","pseudo_4",
"pseudo_15","pseudo_6","pseudo_20","dist_47","dist_17","pseudo_1","dist_77","pseudo_10","dist_91","pseudo_16","pseudo_9",
"dist_16","dist_2","pseudo_17","dist_72","pseudo_7","dist_78","dist_1","pseudo_19","comp_16","ss_1","dist_26","ss_3",
"tran_18","dist_29","pseudo_8","dist_69","aac_8","dist_23","pseudo_3","comp_15","dist_7","dist_20","pseudo_12","dist_10",
"comp_4","dist_70","comp_5","dist_58","comp_18","dist_4","comp_2","dist_21","dist_56","ss_2","comp_10","dist_96","dist_105",
"dist_90","dist_44","dist_53","tran_6","dist_85","pseudo_18","dist_63","dist_67","dist_50","dist_102","dist_13","dist_79",
"dist_61","dist_82","comp_17","dist_59","comp_11","dist_97","dist_64","dist_14","dist_87","dist_99","ss_6","dist_73",
"tran_16","dist_88","dist_34","dist_24","dist_94","tran_17","dist_6","dist_28","dist_27","dist_68","dist_84","comp_1",
"tran_2","pseudo_5","tran_3","tran_10","tran_21","comp_3","tran_20","dist_81","dist_15","tran_19","ss_4","dist_89",
"dist_19","dist_100","comp_6","tran_1","dist_38","dist_103","aac_1","dist_71","dist_8","dist_22","tran_5","comp_13",
"comp_19","dist_83","dist_66","dist_30","dist_37","dist_49","dist_52","dist_55","dist_5","dist_40","dist_65","dist_9"]
        
    for i in range(1,21):
        a="aac_%s"%i
        if a not in features_selected:
            eliminate.append(a)
            
    for i in range(1,401):
        a="dipep_%s"%i
        if a not in features_selected:
            eliminate.append(a)
    
    for i in range(1,31):
        a="pseudo_%s"%i
        if a not in features_selected:
            eliminate.append(a)        
    
    for i in range(1,41):
        a="amphipseudo_%s"%i
        if a not in features_selected:
            eliminate.append(a)
    
    for i in range(1,22):
        a="comp_%s"%i
        if a not in features_selected:
            eliminate.append(a)
            
    for i in range(1,22):
        a="tran_%s"%i
        if a not in features_selected:
            eliminate.append(a)
    
    for i in range(1,106):
        a="dist_%s"%i
        if a not in features_selected:
            eliminate.append(a)  
    
            
    for i in range(1,7):
        a="ss_%s"%i
        if a not in features_selected:
            eliminate.append(a)
    
    all_generated_features=all_generated_features.drop(eliminate, axis=1)
    
    if svalue==0:
        #all_generated_features=all_generated_features[['S5', 'D80', 'C21', 'S4', 'T10', 'D81']]
        all_generated_features=all_generated_features[["pseudo_11","pseudo_2","dist_3","dist_76","dist_18","dist_62","dist_93","pseudo_13","dist_75","pseudo_14","pseudo_4",
"pseudo_15","pseudo_6","pseudo_20","dist_47","dist_17","pseudo_1","dist_77","pseudo_10","dist_91","pseudo_16","pseudo_9",
"dist_16","dist_2","pseudo_17","dist_72","pseudo_7","dist_78","dist_1","pseudo_19","comp_16","ss_1","dist_26","ss_3",
"tran_18","dist_29","pseudo_8","dist_69","aac_8","dist_23","pseudo_3","comp_15","dist_7","dist_20","pseudo_12","dist_10",
"comp_4","dist_70","comp_5","dist_58","comp_18","dist_4","comp_2","dist_21","dist_56","ss_2","comp_10","dist_96","dist_105",
"dist_90","dist_44","dist_53","tran_6","dist_85","pseudo_18","dist_63","dist_67","dist_50","dist_102","dist_13","dist_79",
"dist_61","dist_82","comp_17","dist_59","comp_11","dist_97","dist_64","dist_14","dist_87","dist_99","ss_6","dist_73",
"tran_16","dist_88","dist_34","dist_24","dist_94","tran_17","dist_6","dist_28","dist_27","dist_68","dist_84","comp_1",
"tran_2","pseudo_5","tran_3","tran_10","tran_21","comp_3","tran_20","dist_81","dist_15","tran_19","ss_4","dist_89",
"dist_19","dist_100","comp_6","tran_1","dist_38","dist_103","aac_1","dist_71","dist_8","dist_22","tran_5","comp_13",
"comp_19","dist_83","dist_66","dist_30","dist_37","dist_49","dist_52","dist_55","dist_5","dist_40","dist_65","dist_9"]]
    
    all_generated_features.to_csv("%s.csv" %file_name, header=True, index=False)
    



#predicting antiviral peptide sequences
def predict_BAC_sequences(svalue):
    Delete()
    #reading training and test datsets 
    
    dir_path = os.getcwd()
    dir_path = dir_path.replace("\\", "/")
    if svalue==0:
        training_file_path=dir_path + "/ML-selected_training_merged_file.csv" #path for aac training set
    
    
    testing_file_path=dir_path + "/input_seq.csv" #path for input sequences
    
    #print(training_file_path)
    #print(testing_file_path)
    
    #base = importr('base')
    #utils = rpackages.importr('utils')
    #utils.chooseCRANmirror(ind=1) # select the first mirror in the list
    #packnames = ('ROSE', 'e1071', 'caret') 
    #utils.install_packages(StrVector(packnames))
    r=ro.r
    r.source("file_with_pkgTest.R") #check whether required R package is installed
    #r.pkgTest("ROSE")
    r.pkgTest("base")
    r.pkgTest("e1071") #install e1071 R package if not installed
    r.pkgTest("caret")
    
    #r.pkgTest("protr")
    #r.source("SVM_prediction.R") #calling R script to get prediction results from SVM
    #r.source("file_with_pkgTest_decipher.R")
    #r.pkgTest_decipher("DECIPHER")
    r.source("SVM_classifier_BAC_new.R") #calling R script to get prediction results from SVM
    predictions = r.predict_results(training_file_path,testing_file_path)
    #print(predictions)
    #print(predictions[1])
    
    
    # #################################################################
    lines = [line.rstrip('\n') for line in open('input_seq.fasta')]
    #print(lines)

    #reading all sequences in sequences and their names in sequences_name
    save_sequences= ''
    sequences_input= []
    sequences_name_input= []


    for l in lines:
        if(l == ""): #blank lines are disregarded
                pass
        elif (l[0] == '>'):
            sequences_name_input.append(l[1:])
            sequences_input.append(save_sequences)
            save_sequences= ''
        else:
            save_sequences+= l

    sequences_input.append(save_sequences)
    del sequences_input[0]
    
   
    
 
   #generate prediction statistics

    predict_file = open("predicted_bacteriocin_sequences.fasta", "w+")

    count_resistance_sequences=0   
    for i in range(len(predictions)):
        if (predictions[i]==1 ):
            predict_file.write(str(sequences_name_input[i])+ '\n'+ '\n')
            count_resistance_sequences+=1
    predict_file.close()
    del predict_file


    one_line="Total number of predicted bacteriocin sequences = "+str(count_resistance_sequences) + "\n"  +"\n" 
    with open("predicted_bacteriocin_sequences.fasta", 'r+') as fp:
        lines = fp.readlines()     
        lines.insert(0, one_line)  
        fp.seek(0)                 
        fp.writelines(lines)  
        
     #prediction probability
     #r.source("svm_with_percent_identity.R") #calling R script to get prediction results from SVM
     #predictionsR = r.predict_resultsP(training_file_path,testing_file_path)
    r.source("svm_with_percent_identity.R") #calling R script to get prediction results from SVM
    predictionsR = r.predict_resultsP(training_file_path,testing_file_path)
    #df = pd.read_csv(predictionsR) 
    prob_values=pd.DataFrame(predictionsR)  
    row = prob_values.transpose()
    #print(row)
    row=row.rename(columns={ 0: "Positive", 1: "Negative"})
    #print(row)
    #print(sequences_name_input)
    df = pd.DataFrame(sequences_name_input)
    df=df.rename(columns={ 0: "NCBI no"})
    #print(df)
    concatenated_df =pd.concat([df, row], axis=1, join="inner")
    #print(concatenated_df)
    #concatenated_df = concatenated_df.rename(columns={ 0: "NCBI", '0': "Positive", 1: "Negative"})
    #print(concatenated_df)

    concatenated_df.to_csv('probability_results.csv', header=True, index=False)



#include new sequences to the training and test datasets 

def add_new_sequences(svalue,sgval): 
    file_new = open("input_seq.fasta", "r")
    #data_new = file_new.read()
    file_new.close()
    
    
        
    feature_extraction('input_seq',svalue)
    df1 = pd.read_csv('input_seq.csv')
    if sgval==1:
        Otpt = [1] * df1.shape[0]
    else:
        Otpt = [-1] * df1.shape[0]
    df1["Output"]= Otpt
    df1.to_csv("seq_excld_header.csv", header=False, index=False)
    file_new_features = open("seq_excld_header.csv", "r")
    data_new_features = file_new_features.read()
    file_new_features.close()
    
    if svalue==0:
        file_all_features= open("ML-selected_training_merged_file.csv","a")
    
    
    file_all_features.write(data_new_features)
    file_all_features.close()

    # reset training and test data sets
def restore_training_data():
    
    copyfile("ML-selected_training_merged_file_actual.csv", "ML-selected_training_merged_file.csv")
    

    

    
#build graphical user interface

root = tkinter.Tk()
root.title("BaPreS")
root.geometry("530x570")
#520x485
root.configure(background="Light blue")
#root.wm_attributes('-alpha', 0.7)
#peach puff

#canvas = Canvas(root, width=235, height = 139)  
#canvas = Canvas(root, width=159, height = 173, bg='Light blue',highlightthickness=0)
canvas = Canvas(root, width=50, height = 5, bg='Light blue',highlightthickness=0)
#canvas = Canvas(root, width=362, height = 352, bg='white smoke',highlightthickness=0)
#canvas.pack(padx=1, pady=1)
#canvas.pack(padx=1, pady=1)

# Create label
l = Label(root, text = "BaPreS: A Software Tool for Predicting \n Bacteriocin Protein Sequences")
l.config(font =("Courier", 14))
l.pack(padx=5, pady=5)

#T = Text(root, font="none 12 bold",bd=0,height=2, width=37, padx=0, pady=0)


#T.pack(padx=0, pady=0)


#T.insert(END, "      BacPred: A Software Tool for Predicting \n                Bacteriocin Peptide Sequence","center")

#for file upload


def openFile():
    filepath = filedialog.askopenfilename(initialdir="C:\\Users\\Cakow\\PycharmProjects\\Main",
                                          title="Open Input Sequence File?",
                                          filetypes= (("fasta files","*.fasta"),
                                          ("all files","*.*")))
    
    #print(filepath)
    filepath = open(filepath,'r')
    content= filepath.read()
    my_text.insert(END,content)
    saveFile()
    filepath.close()
    
#end file upload and save

#file save 
def saveFile():
    text_file=open('input_seq.fasta','w')
    text_file.write(my_text.get(1.0,END))
    my_text.delete('1.0', END)


#end file save

def Delete():
    text_read.delete('1.0', END)
    
def predictResult():
    Delete()
    #filepath = filedialog.askopenfilename(initialdir="predicted_AVP_sequences.fasta")
    filepath = open("predicted_bacteriocin_sequences.fasta",'r')
    content= filepath.read()
    text_read.insert(END,content)
    filepath.close()
    
def predictProbabiliteResult():    
    Delete()
    df = pd.read_csv('probability_results.csv')
    #filepath = open("probability_results.csv",'r')
    #content= filepath.read()
    #table = df.to_markdown(tablefmt="pipe", index=False)
    table = df.to_string(index=False, float_format='%.4f', justify='center')
    content=table
    text_read.insert(END,content)
    

    
def setup_window(soption):
    window = Toplevel(root)
    window.geometry("320x60")
    window.configure(background="tan")
    if soption not in choices:
        Lbl = Label(window, bg="tan",fg="black",text="Please Select an Option")
        Lbl.config(font=('Helvetica', 8, 'bold'))
        Lbl.pack( )
    else:
        Lbl = Label(window, bg="tan",fg="black",text="Done!")
        Lbl.config(font=('Helvetica', 8, 'bold'))
        Lbl.pack( )
    
    Btn=Button(window, text="OK", command=window.destroy) 
    Btn.config(font=('Helvetica', 8, 'bold'))
    Btn.pack()
  

    
#peform operation based on choice    
def Send():
   
    sf = "%s" % var.get()
    if var.get()==choices[0]:
        op_value=0
        feature_extraction('input_seq',op_value)
        predict_BAC_sequences(op_value)
        
    
    elif var.get()==choices[1]:
        op_value=0
        signvalue=1
        add_new_sequences(op_value,signvalue)
    
    elif var.get()==choices[2]:
        op_value=0
        signvalue=-1
        add_new_sequences(op_value,signvalue)
    
    
    elif var.get()==choices[3]:
        restore_training_data()
    
    
    setup_window(sf)
   
    
#file open 
button = Button(text="Open and Save Input Sequences",command=openFile)
button.config(font=('Helvetica', 8, 'bold'))
button.pack(padx=10, pady=10)



#data save
#button = Button(text='Save',command=saveFile)
#button.config(font=('Helvetica', 8, 'bold'))
#button.pack(padx=10, pady=10)
#data save eb=nd

var = tkinter.StringVar(root)
# initial value
var.set('< Please Select Option >')
choices = ['Predict Bacteriocin Sequences',            'Add New  Bacteriocin Sequences',             
           'Add New Non-bacteriocin Sequences',  \
           
           'Restore Training Set']
option = tkinter.OptionMenu(root, var, *choices)
#option.config(bg = "GREEN")
#helv35=font.Font(family='Helvetica', size=36)
option.config(font=('Helvetica', 8, 'bold')) 
#option["menu"].config(bg="GREEN")
option["menu"].config(font=('Helvetica', 8, 'bold'))
option.pack( padx=10, pady=40)
button = tkinter.Button(root, text="Submit", command=Send)
button.config(font=('Helvetica', 8, 'bold'))
button.pack(padx=10, pady=10)


#text box for file save data

frame = Frame(root)
text_read=Text(
    frame,
    font=("Helvetical",8),
    height=15,
    width=45,
    wrap='word',
)
text_read.pack(side=LEFT,expand=True)

button1 = tkinter.Button(text='Prediction Result',command=predictResult)


button1.config(font=('Helvetica', 8, 'bold'))
button1.pack(side=tkinter.LEFT, padx=5, pady=5)

button2 = tkinter.Button(text='Probability Result',command=predictProbabiliteResult)
button2.config(font=('Helvetica', 8, 'bold'))
button2.pack(side=tkinter.RIGHT, padx=5, pady=5)

#button2 = tkinter.Button(text='Detele Text Data',command=Delete)
#button2.config(font=('Helvetica', 8, 'bold'))
#button2.pack(side=tkinter.RIGHT, padx=10, pady=10)

sb = Scrollbar(frame)
sb.pack(side=RIGHT, fill=BOTH)

text_read.config(yscrollcommand=sb.set)
sb.config(command=text_read.yview)

frame.pack(expand=True)



my_text = Text(
    #frame,
    #wrap='word',
    font=("Helvetical",8)
)

#my_text.pack(side=LEFT,expand=True)
my_text.pack_forget()




#my_text = Text(root, width=50, height=5, font=("Helvetical",8))
#my_text.pack()


#end saved opened file

# Python program to create
# a file explorer in Tkinter

# import all components
# from the tkinter library
#from tkinter import *

# import filedialog module
#from tkinter import filedialog

# Function for opening the
# file explorer window
#data save



root.mainloop()
#end_time=time.time()
#print("--- %s seconds ---" % (end_time-start_time))


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




