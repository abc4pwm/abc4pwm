#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 14 22:45:13 2019

@author: omerali
"""
from pathlib import Path
import numpy as np
import os
from glob import glob
import shutil
from collections import Counter
from abc4pwm.cleandatabase_for_classification import cleanTFDB
import codecs
import json

class ClassificationPwm(object):
    def __init__(self, in_pwm_folder, out_pwm_folder, original_pwm_folder, new_database = 0 ):
        """

        :param in_pwm_folder: input folder which contain pwms, (transcription factors)
        :param out_pwm_folder: output folder which contain pwms, (transcription factors)
        :param original_pwm_folder: keep copy of input pwms in this folder as input files are renamed and moved to
        output folder
        :param new_database: if you want to load a new database for classification or use old one
        """
        print("\nTask: Classification of TFs based on their DNA Binding Domain")
        self.input_folder_path = in_pwm_folder
        self.output_folder_path = out_pwm_folder
        if not os.path.exists(self.output_folder_path):
            os.makedirs(self.output_folder_path, exist_ok=True)
        self.original_pwm_folder_path = original_pwm_folder
        self.load_new_db = new_database

        leaf_folder = Path(self.input_folder_path)
        folder_for_dbs = leaf_folder.parent
        if self.load_new_db:
            print("Reading new database")
            # pwms = sorted(os.listdir(self.input_folder_path))
            caller = cleanTFDB(self.input_folder_path,1)
        else:
            print("Reading old database")

        self.tfsdb = np.array(json.loads(codecs.open(os.path.join(str(folder_for_dbs), "tfsdb.json"), 'r').read()))

        self.empty_dir(self.output_folder_path)
        self.classify(self)

    @staticmethod
    def folderize_DBDs(self):
    # this function looks for the DBD of file and put files in their respective DBD folder

        mlpfiles = [i.split('/')[-1] for i in glob(os.path.join(self.output_folder_path, "*.mlp"))]
        for pwm in mlpfiles:
            folder_name = pwm.split('-')[-1].split('.mlp')[0]

            path_to_dbd_folder = os.path.join(self.output_folder_path, folder_name)
            if not os.path.exists(path_to_dbd_folder):
                os.makedirs(path_to_dbd_folder)

            old_mlp_path = os.path.join(self.output_folder_path, pwm)
            new_mlp_path = os.path.join(path_to_dbd_folder, pwm)
            shutil.move(old_mlp_path, new_mlp_path)

    @staticmethod
    def moving_files_input(self):
        #moving files to input folder from original files for the next run
        dest_copy = self.input_folder_path
        src_copy = self.original_pwm_folder_path

        pwms = [i.split('/')[-1] for i in glob(os.path.join(src_copy, "*.mlp"))]

        for file_name in pwms:
            full_file_name = os.path.join(src_copy, file_name)
            if os.path.isfile(full_file_name):
                shutil.copy(full_file_name, dest_copy)

    @staticmethod
    def classify(self):
        #this funtion classify the files received in input directory according to the database received in tfsdb
        #also prepares a summary of classifcation in Summary_Classification.txt
        unknowns = []

        pwms = [i.split('/')[-1] for i in glob(os.path.join(self.input_folder_path, "*.mlp"))]
        counterHTF = 0
        counterEW = 0
        counterJP = 0
        total_pwms_in_each_dbd = []

        for i in pwms:
            tf = i.split('_')[0]
            flag = 0

            for j in range(len(self.tfsdb[:,0])):
                if tf in self.tfsdb[j,0]:
                    tf_name_with_dbd = i[:len(i)-4] +'-'+ str(self.tfsdb[j,1]).replace(' ','_').replace('-','_').replace('/','_') + '.mlp'
                    flag = 1
                    total_pwms_in_each_dbd.append(self.tfsdb[j,1])
                    if self.tfsdb[j,2] == 'HTF':
                        counterHTF = counterHTF + 1
                    elif self.tfsdb[j,2] == 'EW':
                        counterEW = counterEW + 1
                    else:
                        pass
                    break


            if (len(i.split('::')))>1:
                double_tfs = []
                double_tfs.append(i.split('::')[0])
                double_tfs.append(i.split('::')[1])
                double_tfs = set(double_tfs)
                for l in double_tfs:
                    for j in range(len(self.tfsdb[:,0])):
                        if l in self.tfsdb[j,0]:
                            tf_name_with_dbd = i[:len(i)-4] +'-'+ str(self.tfsdb[j,1]).replace(' ','_').replace('-','_') +'.mlp'
                            flag = 1
                            if self.tfsdb[j,2] == 'HTF':
                                counterHTF = counterHTF + 1
                            elif self.tfsdb[j,2] == 'EW':
                                counterEW = counterEW + 1
                            elif self.tfsdb[j,2] == 'JP':
                                counterJP = counterJP + 1
                            else:
                                pass

                            break

            if flag ==0:
                unknowns = np.append(unknowns, tf)
                tf_name_with_dbd = i[:len(i)-4] +'-Unknown.mlp'
            tf_name_without_dbd = os.path.join(self.input_folder_path, str(i))
            tf_name_with_dbd = os.path.join(self.output_folder_path, tf_name_with_dbd)
            os.rename(tf_name_without_dbd, tf_name_with_dbd)

        leaf_folder = Path(self.output_folder_path)
        out_dir = leaf_folder.parent
        path_to_text_reports = os.path.join(out_dir,'reports_in_text/')
        if not os.path.exists(path_to_text_reports):
            os.makedirs(path_to_text_reports, exist_ok=True)

        if os.path.exists(os.path.join(path_to_text_reports,'classificationSummary.txt')):
            os.remove(os.path.join(path_to_text_reports,'classificationSummary.txt'))
        with open(os.path.join(path_to_text_reports,'classificationSummary.txt'), 'w') as f:
            f.writelines("Summary of classification \n" +
            "Total Found PWMs " + str(counterHTF+counterEW+counterJP) + "\n" +
            "edgar-wingender : "  + str(counterEW)+ "\n" +
            "HumanTF : "  + str(counterHTF)+ "\n" +
            "JASPAR : " + str(counterJP)+ "\n" +

            str(len(set(self.tfsdb[:,0]))) + ' are number of Unique tfs in all pwms' + "\n" +
            str(len(set(self.tfsdb[:,1]))) + ' are number of Unique DBDs' + "\n" +
            str(Counter(total_pwms_in_each_dbd)))

        self.moving_files_input(self)

        self.folderize_DBDs(self)
        print("Task Completed. \nClassification summary can be seen in text file SummaryClassfication.txt in: ",  path_to_text_reports)
        print("Please see Output files in ", self.output_folder_path)

    @staticmethod
    def empty_dir(folder):
        # function for deleting files from a folder

        for filename in os.listdir(folder):
            file_path = os.path.join(folder, filename)
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                print(('Failed to delete %s. Reason: %s' % (file_path, e)))


if __name__ == "__main__":
    classObj = ClassificationPwm('../data/in/in_pwms','../data/out/classification_out/','../data/in/pwm_original', 0)
#     if len(sys.argv) < 4:
#        print("Usage: python classification.py <input_directory> <output_directory> <original_directory>")
#        exit()
#
#     in_dir = sys.argv[1]
#     out_dir = sys.argv[2]
#     original_dir = sys.argv[3]
#
#     try:
#         new_db_read = int(sys.argv[4])
#     except:
#         new_db_read = 0
#
#     if new_db_read:
#         print("Reading new database")
#         pwms = sorted(os.listdir(in_dir))
#         caller = cleanTFDB(pwms)
#         tfsdb, db_ew, db_htf = caller.get_tfsdb()
#
#     else:
#         tfsdb = np.array(json.loads(codecs.open(os.path.join('../data/in/',"tfsdb.json"), 'r').read()))
#         print("Reading old database")
#
#
# empty_dir(out_dir)
# classify(in_dir, out_dir, original_dir, tfsdb)