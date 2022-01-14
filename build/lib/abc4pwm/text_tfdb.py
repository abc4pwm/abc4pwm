#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 11:20:59 2019

@author: omerali
"""
import sys, os
from abc4pwm.cleandatabase_for_classification import cleanTFDB
import numpy as np


class make_txt_of_clean_databse():
    def __init__(self, path_to_input_pwms, output_path):
        """

        :param path_to_input_pwms: folder which contain pwms of database interies
        """
        print("\nPreparing Database file in text.")
        if not os.path.exists(output_path):
            os.makedirs(output_path, exist_ok=True)
        self.creat_database(path_to_input_pwms, output_path)
        print("Task Completed. \nDatabase file in text can be seen in data/", output_path )

    def creat_database(self,path_to_input_pwms,output_path):

        pwms = sorted(os.listdir(path_to_input_pwms))


        obj_cleanTFDB = cleanTFDB(path_to_input_pwms,0)
        tfsdb, for_txt_tfsEW, for_txt_tfsFromHTF_NP = obj_cleanTFDB.get_tfsdb(path_to_input_pwms)

        tfs = obj_cleanTFDB.make_TFsList_from_filenames(pwms)

        # obj_text = codecs.open('tfsEW.json', 'r', encoding='utf-8').read()
        # b_new = json.loads(obj_text)
        # for_txt_tfsEW = np.array(b_new)
        #
        #
        #
        # obj_text = codecs.open('tfsHTF.json', 'r', encoding='utf-8').read()
        # b_new = json.loads(obj_text)
        # for_txt_tfsFromHTF_NP = np.array(b_new, dtype='object')



        for_txt_replacmentsHTF = {
                       "Basic helix-loop-helix factors (bHLH)": "bHLH",
                       "Basic leucine zipper factors (bZIP)" : "bZIP",
                       "Nuclear receptors with C4 zinc fingers" : "Nuclear receptor",
                       "T-Box factors" : "T-box",
                       "Homeo domain factors" : "Homeodomain",
                       "High-mobility group (HMG) domain factors" : "HMG",
                       "Fork head / winged helix factors" : "Forkhead",
                       "Grainyhead domain factors" : "Grainyhead",
                       "p53 domain factors" : "p53",
                       "Fork head - winged helix factors" : "Forkhead",
                       "zinc factor" : "ZF",
                       "zinc finger-type factors" : "ZF" ,
                       "SMAD/NF-1 DNA-binding domain factors" : "SMAD",
                       "SAND domain factors" : "SAND",
                       "C2H2 zinc finger factors" : "C2H2 ZF",
                       "Paired box factors" : "Paired box",
                       "Tryptophan cluster factors" : "Tryptophan",
                       "C2CH THAP-type zinc finger factors" : "C2CH_THAP-type ZF",
                       "Runt domain factors" : "Runt",
                       "Rel homology region (RHR) factors" : "RHR",
                       "SMAD-NF-1 DNA-binding domain factors" : "SMAD",
                       "Other C4 zinc finger-type factors" : "Other C4 ZF",
                       "DM-type intertwined zinc finger factors" : "DM-type intertwined ZF",
                       "TEA domain factors" : "TEA",
                       "Basic helix-span-helix factors" : "bHSH"


                       }



        for k,v in for_txt_replacmentsHTF.items():
            for l in range(len(for_txt_tfsFromHTF_NP[:,1])):
                if v == for_txt_tfsFromHTF_NP[l,1]:
                    for_txt_tfsFromHTF_NP[l,1] = for_txt_tfsFromHTF_NP[l,1].replace(v,k)


        for k,v in for_txt_replacmentsHTF.items():
            for l in range(len(for_txt_tfsEW[:,1])):
                if v == for_txt_tfsEW[l,1]:
                    for_txt_tfsEW[l,1] = for_txt_tfsEW[l,1].replace(v,k)

        for_txt_tfsFromJaspar = [['EWSR1::FLI1','Tryptophan cluster factors','JP'],
                         ['FLI1','Tryptophan cluster factors','HTF'],
                         ['HLTF','Tryptophan cluster factors','JP'],
                         ['SPZ1','Basic helix-loop-helix factors (bHLH)','JP'],
                         ['TATA','TATA','HTF'],
                        ]
        for_txttfsdb = np.concatenate((for_txt_tfsEW, for_txt_tfsFromHTF_NP,for_txt_tfsFromJaspar), axis=0)
        txttfsdb = []


        for tf in tfs:                             #iterate through list
            tf = str(tf)
            for j in range(len(tfsdb[:,0])):               #Iterate through DB to compare TF with
                if tf in tfsdb[j,0]:                       #if matches
                    txttfsdb.append([tf, tfsdb[j,2] , tfsdb[j,1] , for_txttfsdb[j,1] ])          #prepare destination file name
        #            print tfsdb[j,1], for_txttfsdb[j,1][7:]
                    break                                                #don't need to check remaining (optimization)








        with open(os.path.join(output_path,"txttfsdb.txt"), "w") as txt_file:
            txt_file.writelines("TF Name |\tSource Database |\tDBD Short name |\tDBD Full name \n")
            for line in txttfsdb:
                txt_file.writelines("\t".join(line) + "\n")


# if __name__ == "__main__":
#     # if len(sys.argv) < 2:
#     #     print("Usage: python text_tfdb.py <path_to_input_pwms>")
#     #     exit()
#
    # textDbObj = make_txt_of_clean_databse('../data/in/in_pwms/')