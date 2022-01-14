#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 11:30:31 2019

@author: omerali
"""
from pathlib import Path
import numpy as np
import pandas as pd
import os
import glob
from copy import deepcopy
from bs4 import BeautifulSoup as BS
import requests
import json
import codecs
import urllib

class edgar_wingender_scraper():
    #this class scrapes the webpage http://www.edgar-wingender.de/huTF_classification.html
    # and returns the data in tfdb
    tfdb = []
    currentSclass = ''
    currentclass = ''
    currentfamily = ''
    currentsubfamily = ''

    def __init__(self, url):
        website_url = requests.get(url).text
        soup = BS(website_url, 'lxml')
        My_table = soup.find('table', {'class': ''})
        self.all_trs = My_table.findAll('tr')[4:]

    def checktr(self, tr):
        #this fucntion check every row if it is subgamily or class or subclass etc
        #    print(CN)
        CN = list(tr.attrs.values())[0][0]
        td_list = tr.findAll('td')
        if CN != 'genus_tr':
            if CN == 'family_tr':
                self.currentfamily = td_list[2].get_text()
            if CN == 'subfamily_tr':
                self.currentsubfamily = td_list[2].get_text()
            if CN == 'class_tr':
                self.currentclass = td_list[2].get_text().replace('/', '-')
            if CN == 'superclass_tr':
                self.currentSclass = td_list[1].get_text()
            return False
        else:
            return True

    def scrape(self):
        for tr in self.all_trs:
            try:
                if self.checktr(tr):
                    td_list = tr.findAll('td')
                    self.tfdb.append([td_list[2].get_text().replace('\\u03b1', ''),
                                      self.currentSclass,
                                      self.currentclass,
                                      self.currentfamily,
                                      self.currentsubfamily
                                      ])
            except:
                pass
        return self.tfdb


class cleanTFDB():
    # this class prepares a clean database from two sources
    # 1. edgar wingender 2. Human transcription factor
    # return 1. complete TF database, 2. TF databse form source edgar wingender, 3. from Human TF


    replacmentsEW = {}
    replacmentsHTF = {}
    tfsFromEW = []
    tfsUnknown = []
    tfsFromHTF = []

    def __init__(self, input_folder_path, read_new=0):
        """
        :param input_folder_path: folder which contain input files (pwms) in  .mlp
        :param read_new: if want new database files downloaded from internet. (internet connection required)
        """
        print("\nTask: Uniformly named database for Classification")


        pwms = glob.glob(os.path.join(input_folder_path,'*.mlp'))

        leaf_folder = Path(input_folder_path)
        folder_for_dbs = leaf_folder.parent

        tfs = self.make_TFsList_from_filenames(pwms)

        DB_HTF, DB_EW = self.read_source_databases(folder_for_dbs,read_new)
        self.checkdbtf(tfs, DB_EW, DB_HTF)

        self.replacmentsEW = {"Basic helix-loop-helix factors (bHLH)": "bHLH",
                              "Basic leucine zipper factors (bZIP)": "bZIP",
                              "Nuclear receptors with C4 zinc fingers": "Nuclear receptor",
                              "T-Box factors": "T-box",
                              "Homeo domain factors": "Homeodomain",
                              "High-mobility group (HMG) domain factors": "HMG",
                              "Fork head / winged helix factors": "Forkhead",
                              "Grainyhead domain factors": "Grainyhead",
                              "p53 domain factors": "p53",
                              "Fork head - winged helix factors": "Forkhead",
                              "zinc factor": "ZF",
                              "zinc finger-type factors": "ZF",
                              "SMAD/NF-1 DNA-binding domain factors": "SMAD",
                              "SAND domain factors": "SAND",
                              "C2H2 zinc finger factors": "C2H2 ZF",
                              "Paired box factors": "Paired box",
                              "Tryptophan cluster factors": "Tryptophan",
                              "C2CH THAP-type zinc finger factors": "C2CH_THAP-type ZF",
                              "Runt domain factors": "Runt",
                              "Rel homology region (RHR) factors": "RHR",
                              "SMAD-NF-1 DNA-binding domain factors": "SMAD",
                              "Other C4 zinc finger-type factors": "Other C4 ZF",
                              "Basic helix-span-helix factors": "bHSH",
                              "C2HC zinc finger factors": "C2HC ZF",
                              "STAT domain factors": "STAT",
                              "MADS box factors": "MADS",
                              "DM-type intertwined zinc finger factors": "DM-type intertwined ZF",
                              "TEA domain factors": "TEA"}

        self.replacmentsHTF = {"E2F": "Forkhead",
                               "Ets": "Tryptophan",
                               "AP-2": "bHSH",
                               "Tryptophan clust": "Tryptophan",
                               "CBF-NF-Y": "Heteromeric CCAAT-binding factors",
                               "Ets; AT hook": "Tryptophan",
                               "C2H2 zinc finger": "C2H2 ZF",
                               "Myb-SANT": "Tryptophan",
                               "Rel": "RHR",
                               "HMG/Sox":"HMG",
                               "MADS box": "MADS",
                               "CUT; Homeodomain": "Homeodomain",
                               "Homeodomain; Paired box": "Homeodomain",




                               }

        self.replacements()
        self.get_tfsdb(folder_for_dbs)
        print("Task Completed. \nPlease see new created uniformly named database in tfsdb.json in ", folder_for_dbs)


    def read_source_databases(self, folder_for_dbs, read_new = 0):
        if read_new:
            print("Preparing new dataset from online resources.")
            url = 'http://www.edgar-wingender.de/huTF_classification.html'
            driver = edgar_wingender_scraper(url)
            df_EW = pd.DataFrame(driver.scrape())

            DB_EW = deepcopy(df_EW[[0, 2]])
            DB_EW.loc[:, '3'] = 'EW'
            DB_EW.to_json(os.path.join(folder_for_dbs,'ew_db.json'))
            DB_EW = DB_EW.to_numpy()



            #Human TF: as HTF       Egdar Wingender = EW

            url = 'http://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.txt'

            testfile = urllib.request.urlretrieve(url, os.path.join(folder_for_dbs,'DB_HumanTF.txt'))
        else:
            if not os.path.exists(os.path.join(folder_for_dbs,'ew_db.json')):
                print('Any input file of Edgar Wingender classification cannot be found'
                      'Please select True for read_new to create a new db.')
                exit()
            else:
                print("Preparing dataset from offline files in ", folder_for_dbs)
                DB_EW = pd.read_json(os.path.join(folder_for_dbs,'ew_db.json'))
                DB_EW = DB_EW.to_numpy()

        df_HTF = pd.read_csv(os.path.join(folder_for_dbs,'DB_HumanTF.txt'), header=None, delimiter='\t')
        df_HTF = df_HTF.drop_duplicates()  # drop if there are any duplicates

        DB_HTF = deepcopy(df_HTF[[2, 3]])

        DB_HTF.loc[:, '4'] = 'HTF'
        DB_HTF = DB_HTF.to_numpy()

        return DB_HTF, DB_EW

    def make_TFsList_from_filenames(self, pwms):
        tfs = []

        for i in pwms:
            tfs.append(i.split('/')[-1].split('_')[0])
        return tfs

    def checkdbtf(self, tfs, DB, DB_HTF):
        for i in tfs:
            flag = 0
            #            i = i.decode('utf-8')
            for j in range(len(DB[:, 0])):
                tfstr = DB[j, 0].replace('\\u03b1', '')
                if i in tfstr:
                    self.tfsFromEW.append([tfstr, DB[j, 1][7:], DB[j, 2]])

                    flag = 1
                    break
            if not flag:
                self.tfsUnknown.append(i)

        for i in self.tfsUnknown:
            for j in range(len(DB_HTF[:, 0])):
                tfstr = DB_HTF[j, 0].replace('\\u03b1', '')
                if i in tfstr:
                    self.tfsFromHTF.append(([tfstr, DB_HTF[j, 1], DB_HTF[j, 2]]))
                    break

        self.tfsFromEW_NP = np.asarray(self.tfsFromEW, dtype=None, order=None)
        self.tfsFromHTF_NP = np.asarray(self.tfsFromHTF, dtype=None, order=None)

        return self.tfsFromEW_NP, self.tfsFromHTF_NP

    def replacements(self):

        for k, v in list(self.replacmentsEW.items()):
            for l in range(self.tfsFromEW_NP.shape[0]):
                if k == self.tfsFromEW_NP[l, 1]:
                    self.tfsFromEW_NP[l, 1] = self.tfsFromEW_NP[l, 1].replace(k, v)

        for k, v in list(self.replacmentsHTF.items()):
            for l in range(self.tfsFromHTF_NP.shape[0]):
                if k == self.tfsFromHTF_NP[l, 1]:
                    self.tfsFromHTF_NP[l, 1] = self.tfsFromHTF_NP[l, 1].replace(k, v)



    def get_tfsdb(self, folder_for_dbs):

        tfsFromJaspar = [['EWSR1::FLI1', 'Tryptophan', 'JP'],
                         ['FLI1', 'Tryptophan', 'JP'],
                         ['HLTF', 'Tryptophan', 'JP'],
                         ['SPZ1', 'bHLH', 'JP'],
                         ['TATA', 'TATA', 'HTF'],
                         ]
        tfsFromJaspar = np.array(tfsFromJaspar)

        self.tfsdb = np.concatenate((self.tfsFromEW_NP, self.tfsFromHTF_NP, tfsFromJaspar))

        json.dump(self.tfsdb.tolist(), codecs.open( os.path.join(folder_for_dbs,"tfsdb.json"), "w"))

        return self.tfsdb, self.tfsFromEW_NP, self.tfsFromHTF_NP


if __name__ == "__main__":
    # if len(sys.argv) < 2:
    #     print("Usage: python cleandatabase_for_classification_pwm.py <inputdir>")
    #     exit()
    #
    # in_dir = sys.argv[1]
    # #    in_dir = 'in/'
    #
    # # url = 'http://www.edgar-wingender.de/huTF_classification.html'
    # # driver = edgar_wingender_scraper(url)

    caller = cleanTFDB('../data/in/in_pwms', 0)