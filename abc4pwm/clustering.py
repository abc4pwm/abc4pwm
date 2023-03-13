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
import json
from glob import glob
import warnings
warnings.filterwarnings("ignore")
from distutils.dir_util import copy_tree
import multiprocessing as mp
from abc4pwm.convert_count_to_pwm import motif_weight2p
from abc4pwm.similarity_score import compute_similarity_score4alignment
from abc4pwm.energy_to_p import read_energy_matrix
from abc4pwm.non_dbd_clustering import non_dbd_ClusteringPwm

from functools import partial


from sklearn.cluster import AffinityPropagation

class ClusteringPwm():

    def __init__(self, input_folder_path,
                 output_folder_path,
                 path_to_text_report='default',
                 in_dbd = True,
                 minimum_pwms_in_dbd=5,
                 max_no_processors =5,
                 seed = 0,
                 damp=0.5,
                 max_iter=400,
                 convergence_iter=30,
                 preference=None):
        """

        :param input_folder_path: this should point to folder which contain DBD folders
        :param output_folder_path: this should point to folder which contain clustered DBD folders
        :param in_dbd: this should be true if you want clustering inside dbd. Otherwise false.
        :param minimum_pwms_in_dbd: a dbd having less than this number of pwm will not be clustered
        :param max_no_processors: for parallel processing, select maximum number of processors. Default is 5
        :param damp: Damping factor (between 0.5 and 1) is the extent to which the current value is maintained relative to incoming values (weighted 1 - damping).
        This in order to avoid numerical oscillations when updating these values (messages).
        :param max_iter: default=200, Maximum number of iterations.
        :param convergence_iter: default=15 Number of iterations with no change in the number of estimated clusters that stops the convergence.
        :param preference: default=None Preferences for each point - points with larger values of preferences are more likely to be chosen as exemplars.
        The number of exemplars, ie of clusters, is influenced by the input preferences value. If the preferences are not passed as arguments, they will be set to the median of the input similarities.
        """

        if in_dbd:
            print("\nTask: Clustering of TFs based on their DNA Binding Domain")
        else:
            print("\nTask: Clustering of TFs")

        input_folder_path = os.path.join(input_folder_path, '')
        self.output_folder_path = os.path.join(output_folder_path, '')
        leaf_folder = Path(self.output_folder_path)
        out_dir = leaf_folder.parent

        if not os.path.exists(output_folder_path):
            os.makedirs(output_folder_path, exist_ok=True)

        if 'default' in path_to_text_report:
            path_to_text_reports = os.path.join(out_dir, "reports_in_text")
        else:
            path_to_text_reports = path_to_text_report
        if not in_dbd:
            clusteringClassobj = non_dbd_ClusteringPwm(input_folder_path, output_folder_path, path_to_text_reports, seed, damp, max_iter, convergence_iter, preference)
            exit()
        self.empty_dir(output_folder_path)
        self.minimum_pwms_in_dbd = minimum_pwms_in_dbd
        self.max_processors = max_no_processors
        self.total_clusters = 0
        self.unclustered_dbds = 0
        self.unclustered_pwms = 0

        copy_tree(input_folder_path,output_folder_path)
        input_folder_path = output_folder_path

        dbds = sorted(os.listdir(input_folder_path))
        for ind, i in enumerate(dbds):
            if i.startswith('.DS'):
                dbds.pop(ind)
        for i in dbds:
            path_to_dbd = os.path.join(input_folder_path, i)
            self.drive_clustering(self, path_to_dbd, seed, damp, max_iter, convergence_iter, preference)



        if os.path.exists(os.path.join(path_to_text_reports,"clusterSummary.txt")):
            os.remove(os.path.join(path_to_text_reports,"clusterSummary.txt"))
            os.makedirs(path_to_text_reports, exist_ok=True)
        if not os.path.exists(path_to_text_reports):
            os.makedirs(path_to_text_reports, exist_ok=True)
        with open(os.path.join(path_to_text_reports,"clusterSummary.txt"),'w') as cs:
            cs.writelines("Clustering Time: " + str(strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime()))+ "\n"
                          "Clustering Technique: Affinity Propagation Clustering \n"
                          "Minimum PWMs in a DBD considered for Clustering: " + str(self.minimum_pwms_in_dbd) + "\n\n"
                          "Total DBDs : " + str(len(dbds)) + "\n"
                          "Total Clusters made : " +  str(self.total_clusters) + "\n"
                          "UnClustered DBDs due to less than threshold : " +  str(self.unclustered_dbds) + "\n"
                          "UnClustered pwms in DBDs less than threshold: " +  str(self.unclustered_pwms) + "\n"
                          "Damping: " + str(damp) + "\n"
                          "Preference: " + str(preference) + "\n"
                          "Max Iterations: " + str(max_iter) + "\n"
                          "Convergence Iterations: " + str(convergence_iter) + "\n")


        print("Task completed. \n "
              "Please see clusters in : ", input_folder_path, "<dbd_folder>/out \n"
              "Clustering summary in , ", path_to_text_reports)


    @staticmethod
    def drive_clustering(self, inputdir, seed, damp, max_iter, convergence_iter, preference):
        #this function prepares a similiarity matrix in parallel and send to a function for clustering
        #rename in files according to their clusters and put in respective cluster folder
        #also call representative motif function at the end



        pwms = [i.split('/')[-1] for i in glob(os.path.join(inputdir, "*.mlp"))]
        pwms = sorted(pwms)
        print(self.minimum_pwms_in_dbd)
        print(type(self.minimum_pwms_in_dbd))
        if len(pwms) < int(self.minimum_pwms_in_dbd):

            self.unclustered_dbds+=1
            self.unclustered_pwms+=len(pwms)

            return 1
        else:
            n_processors = int(np.ceil(len(pwms)/30))
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
            calculate_similarity = partial(self.calculate_similarity_matrix, inputdir=inputdir)
            chunks_similarity_matrix = pool.map(calculate_similarity, chunks_of_pwms)
            similarity_matrix = np.concatenate((chunks_similarity_matrix), axis=0)

        clusters_labels = self.clustering(similarity_matrix, seed, damp, max_iter, convergence_iter, preference)



        self.renaming_mlp(self, inputdir, clusters_labels,pwms)
        self.total_clusters += len(np.unique(clusters_labels))

        self.folderizeclusters(os.path.join(inputdir,'out/'))

    @staticmethod
    def folderizeclusters(folder_path):

        mlpfiles = [i.split('/')[-1] for i in glob(os.path.join(folder_path, "*.mlp"))]

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



        full_pwms = [i.split('/')[-1] for i in glob(os.path.join(inputdir, "*.mlp"))]
        full_pwms = sorted(full_pwms)

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
    def clustering(similarityMatrix, seed, damp, max_iter, convergence_iter, preference):
        #function for clustering algorithm
        clusters = AffinityPropagation(damping=damp, max_iter=max_iter, convergence_iter=convergence_iter, preference=preference, random_state=seed , affinity='precomputed', verbose=False).fit(similarityMatrix)
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
    clusteringClassobj = ClusteringPwm('../data/out/classification_out', '../data/out/clustering_out/')
    # clusteringClassobj = ClusteringPwm('../../../abc4pwm_working/AffinityPropogation_Clustering/to_omer/true_peaks/motif_out/neil2-ko/',
    #                                    '../../../abc4pwm_working/AffinityPropogation_Clustering/to_omer/out_cluster/neil2-ko/test1',
    #                                    path_to_text_report='../../../abc4pwm_working/AffinityPropogation_Clustering/to_omer/out_cluster/neil2-ko/test1',
    #                                    in_dbd=0)
