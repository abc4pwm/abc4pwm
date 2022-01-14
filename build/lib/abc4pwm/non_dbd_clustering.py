#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 14:43:24 2019

@author: omerali
"""

import numpy as np
import os, shutil
from pathlib import Path
from time import gmtime, strftime
from glob import glob
import json
from distutils.dir_util import copy_tree
import multiprocessing as mp
from abc4pwm.convert_count_to_pwm import motif_weight2p
from abc4pwm.similarity_score import compute_similarity_score4alignment
from abc4pwm.energy_to_p import read_energy_matrix

from functools import partial


from sklearn.cluster import AffinityPropagation

class non_dbd_ClusteringPwm():

    def __init__(self, input_folder_path, output_folder_path, max_no_processors = 5):
        """

        :param input_folder_path: this should point to folder which contain DBD folders
        :param output_folder_path: this should point to folder which contain clustered DBD folders
        :param max_no_processors: for parallel processing, select maximum number of processors. Default is 5


        """
        if not os.path.exists(output_folder_path):
            os.makedirs(output_folder_path, exist_ok=True)

        self.output_folder_path = output_folder_path
        leaf_folder = Path(self.output_folder_path)
        out_dir = leaf_folder.parent
        # if not os.path.exists(os.path.join(out_dir, 'non_dbd_uncertain_pwms')):
        #     os.makedirs(os.path.join(out_dir, 'non_dbd_uncertain_pwms'), exist_ok=True)
        # if not os.path.exists(output_folder_path):
        #     os.makedirs(output_folder_path, exist_ok=True)

        # self.empty_dir(os.path.join(out_dir, 'non_dbd_uncertain_pwms'))
        self.empty_dir(output_folder_path)
        self.total_clusters = 0
        self.unclustered_dbds = 0
        self.unclustered_pwms = 0
        self.total_bad_pwms = 0
        self.max_processors = max_no_processors

        self.empty_dir(output_folder_path)
        copy_tree(input_folder_path,output_folder_path)
        input_folder_path = output_folder_path

        print("\nClustering TFs (Non DBD) ", input_folder_path)
        self.non_dbd_drive_clustering(self, input_folder_path)


        path_to_text_reports = os.path.join(out_dir, 'reports_in_text/')

        if not os.path.exists(path_to_text_reports):
            os.makedirs(path_to_text_reports)
        if os.path.exists(os.path.join(path_to_text_reports, "non_dbd_clusterSummary.txt")):
            os.remove(os.path.join(path_to_text_reports, "non_dbd_clusterSummary.txt"))
            os.makedirs(path_to_text_reports, exist_ok=True)

        with open(path_to_text_reports+"non_dbd_clusterSummary.txt", 'w') as cs:
            cs.writelines("Clustering Time: " + str(strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime()))+ "\n"
                          "Clustering Technique: Affinity Propagation Clustering \n"
                          "Total Clusters made : " + str(self.total_clusters) + "\n")


        print("Task completed. \n "
              "Please see clusters in : ", input_folder_path, "\n"
              "Clustering summary in data/out/reports_in_text")


    @staticmethod
    def non_dbd_drive_clustering(self, inputdir):
        #this function prepares a similiarity matrix in parallel and send to a function for clustering
        #rename in files according to their clusters and put in respective cluster folder
        #also call representative motif function at the end

        leaf_folder = Path(self.output_folder_path)
        out_dir = leaf_folder.parent
        # dst_for_bad_pwms = os.path.join(out_dir, 'non_dbd_uncertain_pwms/')

        # pwms = sorted([i for i in os.listdir(inputdir) if not os.path.isdir(i)])
        temp = []
        pwms = glob(os.path.join(inputdir, "*.mlp"))
        for i in pwms:
            temp.append(i.split('/')[-1])
        pwms = temp
        n_processors = int(np.ceil(len(pwms)/30)+1)
        if n_processors > self.max_processors:
            n_processors = self.max_processors

        pool = mp.Pool(processes=n_processors)

        start = 0
        processor_capacity = int(np.ceil(len(pwms) / n_processors))
        end = processor_capacity
        chunks_of_pwms = []
        for i in range(n_processors):
            chunks_of_pwms.append(pwms[start:end])
            start = end
            end = end + processor_capacity
            # print(pwms[start:end])
        calculate_similarity = partial(self.calculate_similarity_matrix, inputdir=inputdir)
        chunks_similarity_matrix = pool.map(calculate_similarity, chunks_of_pwms)
        similarity_matrix = np.concatenate((chunks_similarity_matrix), axis=0)

        clusters_labels = self.clustering(similarity_matrix)



        self.renaming_mlp(self, inputdir, clusters_labels,pwms)
        self.total_clusters += len(np.unique(clusters_labels))

        self.folderizeclusters(os.path.join(inputdir,'out/'))

        # self.total_bad_pwms = make_representative_pwm(os.path.join(inputdir,'out/'), False, self.mean_threshold, self.z_score_threshold, self.top_occurrence, self.occurrence_threshold)
        # return self.total_bad_pwms

    @staticmethod
    def folderizeclusters(folder_path):
        #    folder_path = "test_out/"

        mlpfiles = [f for f in sorted(os.listdir(folder_path)) if os.path.isfile(os.path.join(folder_path, f))]

        for pwm in mlpfiles:
            folder_name = pwm.split('_')[0]

            new_path = os.path.join(folder_path, folder_name)
            if not os.path.exists(new_path):
                os.makedirs(new_path)

            old_mlp_path = os.path.join(folder_path, pwm)
            new_mlp_path = os.path.join(new_path, pwm)
            shutil.move(old_mlp_path, new_mlp_path)

    @staticmethod
    def calculate_similarity_matrix(pwms_full, inputdir):
        #this function extracts matrices from files, convert them to probablity and return similarity matrix

        temp = []
        full_pwms = glob(os.path.join(inputdir, "*.mlp"))
        for i in full_pwms:
            temp.append(i.split('/')[-1])
        full_pwms = temp

        df = np.zeros((len(pwms_full), len(full_pwms)))
        for index1, i in enumerate(pwms_full):
            matrix1, matrix_string1, maximum_feq1, total_maximum1, info1 = read_energy_matrix(os.path.join(inputdir, i))
            matrix1 = motif_weight2p(matrix1)

            for index2, j in enumerate(full_pwms):
                matrix2, matrix_string2, maximum_feq2, total_maximum2, info2 = read_energy_matrix(os.path.join(inputdir, j))
                matrix2 = motif_weight2p(matrix2)

                df[index1, index2] = compute_similarity_score4alignment(matrix1, matrix2)


        similarityMatrix = np.asarray(df)
        return similarityMatrix



    @staticmethod
    def read_similarity_matrix():
        #funciton for reading similariy matrix stores in a json file
        with open('similarityMatrix.json') as f:
            similarityMatrix = np.array(json.load(f))
            return similarityMatrix

    @staticmethod
    def clustering(similarityMatrix):
        #function for clustering algorithm
        clusters = AffinityPropagation(damping=0.5, max_iter=400, convergence_iter=30, preference=None, affinity='precomputed', verbose=False).fit(similarityMatrix)
        return clusters.labels_


    @staticmethod
    def renaming_mlp(self, inputdir, clusters,pwms):
        #this function add cluster number to every file


        if not os.path.exists(os.path.join(inputdir,'out/')):
            os.mkdir(os.path.join(inputdir,'out/'))
        outputdir = os.path.join(inputdir,'out/')
        self.empty_dir(outputdir)

        for ind, i in enumerate(pwms):
            src = os.path.join(inputdir,i)  # renaming
            dst = str(clusters[ind]) + '_' + str(i)
            dst = os.path.join(outputdir,dst)

            os.rename(src, dst)

    @staticmethod
    def empty_dir(folder):
        #function for deleting files from a folder

        for filename in os.listdir(folder):
            file_path = os.path.join(folder, filename)
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                print(('Failed to delete %s. Reason: %s' % (file_path, e)))


    @staticmethod
    def moving_files_input(src_copy, dest_copy):
        #moving files to input folder from original files for the next run


        if not os.path.exists(dest_copy):
            os.mkdir(dest_copy,mode=0o777)

        pwms = os.listdir(src_copy)

        for indexpwm, i in enumerate(pwms):
            if i.startswith('.'):
                pwms.pop(indexpwm)

        for file_name in pwms:
            full_file_name = os.path.join(src_copy, file_name)
            if os.path.isfile(full_file_name):
                shutil.move(full_file_name, dest_copy)



if __name__ == "__main__":
    # clusteringClassobj = non_dbd_ClusteringPwm('../data/out/uncertain_pwms/','../data/out/uncertain_clustering_out/')
    clusteringClassobj = non_dbd_ClusteringPwm('ensemble/omer_ensemble/demo_out_ensemble/swi4/demo_in/','ensemble/omer_ensemble/out/clustering_out/')

    # if len(sys.argv) < 2:
    #     print(
    #         "Usage: python clustering.py <input_directory> ")
    #     exit()
    #
    # input_directory = sys.argv[1]
    # empty_dir(os.path.join('../data/out/','bad_pwms'))
    #
    # dbds = sorted(os.listdir(input_directory))
    # for ind, i in enumerate(dbds):
    #     if i.startswith('.DS'):
    #         dbds.pop(ind)
    # print("Clustering...")
    # for i in dbds:
    #     in_directory = os.path.join(input_directory,i)
    #     drive_clustering(in_directory)
    # # in_directory = 'c_out/HMG'
    # # drive_clustering(in_directory)
    #
    # print("See clusters in :  ../data/out/c_out/<dbd_folder>/out/")


