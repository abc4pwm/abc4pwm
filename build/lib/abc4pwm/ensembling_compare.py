import os
import shutil
import glob
import numpy as np
from collections import Counter
from pathlib import Path
from abc4pwm.non_dbd_clustering import non_dbd_ClusteringPwm
from abc4pwm.searching import searching
from abc4pwm.pwm_representative_test import make_representative_pwm

class PredictedCompare():
    def __init__(self, path_to_predicted_files, output_folder, db_folder, min_pwms_in_cluster=3, top_n=2, db_type='folder',input_count=1,db_count=0, db_file_type='.txt',
                     input_file_type='.mlp', input_prob=0, db_prob=1, tf_name=''):
        """
        :param path_to_predicted_files: path to folders of predicted pwms
        :param output_folder: output folder path where results of analysis file and matched TF will be stored
        :param min_pwms_in_cluster: minimum number of pwms in a cluster at clustering step
        :param db_folder: path to folder of database which contains files with which you want to match predicted files.
        :param top_n: define top n(top matches from database) for your predicted file in searching step
        :param db_type: type of database: Options are folder, path(in case of abc4pwm setting)
        :param input_count: True if input file contains counts
        :param db_count: True if database file contains counts
        :param db_file_type: mention the extension of file type. e.g., .mlp, .txt
        :param input_file_type: mention the extension of file type. e.g., .mlp, .txt
        :param input_prop: True if input matrix is already in probability
        :param db_prob: True if database matrix is already in probability
        :param tf_name: TF name for which results are compared with.

        """

        temp_folder = Path(path_to_predicted_files)
        temp2_folder = temp_folder.parent
        dst_folder = os.path.join(temp2_folder, 'in/')
        if not os.path.exists(dst_folder):
            os.makedirs(dst_folder, exist_ok=True)
        self.empty_dir(dst_folder)
        self.prepare_for_clustering(path_to_predicted_files, dst_folder)
        print("step 1: Prepare for clustering Done")
        cluster_out_folder = os.path.join(temp2_folder,'clustering_out/')
        self.do_clustering(dst_folder, cluster_out_folder)
        print("step 2: Clustering Done")

        self.make_rep(os.path.join(cluster_out_folder, 'out/'), dbd='selected', clusters='all', ic=0.4,
                 best_match_initial_motif=1, mean_threshold=0.75, z_score_threshold=-0.9,
                 top_occurrence=0.35, occurrence_threshold=0.25)
        predicted_folder = os.path.join(temp2_folder,'predicted')
        self.ensembling_compare_investigate(os.path.join(cluster_out_folder, 'out/'),
                                       predicted=predicted_folder, min_threshold_for_number_of_pwms=min_pwms_in_cluster)
        print("step 3: Representative Motifs Done")

        self.do_searching(in_folder=predicted_folder,
                     db_folder=db_folder,
                     out_folder=output_folder,
                     db_type=db_type,
                     input_count=input_count,
                     n=top_n,
                     db_count=db_count,
                     db_file_type=db_file_type,
                     input_file_type=input_file_type,
                     input_prob=input_prob,
                     db_prob=db_prob,
                     tf_name=tf_name)
        print("step 4: Searching Done")


    def prepare_for_clustering(self,src, dst):
        if not os.path.exists(dst):
            os.makedirs(dst, exist_ok=True)
        self.empty_dir(dst)
        pwms = glob.glob(src+"/*.mlp")
        for pwm in pwms:
            shutil.copy(pwm,dst)
    def do_clustering(self, in_folder, out_folder):
        non_dbd_ClusteringPwm(in_folder,out_folder)

    def make_rep(self, path_to_clusters,  dbd = 'selected', clusters = 'all' , ic = 0.4 , best_match_initial_motif = 1, mean_threshold= 0.75, z_score_threshold = -0.9,
                     top_occurrence = 0.35, occurrence_threshold=0.25):
        make_representative_pwm(path_to_clusters,  dbd, clusters, ic , best_match_initial_motif, mean_threshold, z_score_threshold ,
                     top_occurrence, occurrence_threshold)

    def empty_dir(self, folder):
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
    def ensembling_compare_investigate(self, path_of_clusters= 'out/clustering_out/out/',predicted = '', min_threshold_for_number_of_pwms= 3):


        clusters_list = [i for i in sorted(os.listdir(path_of_clusters)) if os.path.isdir(os.path.join(path_of_clusters, i))]
        clusters_list = [int(x) for x in clusters_list]
        clusters_list.sort()
        if not os.path.exists(predicted):
            os.makedirs(predicted)
        self.empty_dir(predicted)
        for cluster in clusters_list:
            pwms = glob.glob(os.path.join(path_of_clusters,str(cluster))+"/*.mlp")
            if len(pwms) < int(min_threshold_for_number_of_pwms):
                if not os.path.exists(predicted):
                    os.makedirs(predicted)

            else:
                src = os.path.join(path_of_clusters,str(cluster),"repres",str(cluster)+'_rep.mlp')
                dst = os.path.join(predicted,pwms[0].split('/')[-1]+'_'+str(cluster)+'_rep.mlp')
                os.rename(src,dst)
    def do_searching(self, in_folder,db_folder, out_folder, tf_name ,n,db_type ,input_count,db_count, db_file_type,
                     input_file_type, input_prob, db_prob):
        pwms = glob.glob(in_folder+"/*.mlp")
        all_matched_list = []
        all_score_list = []
        if not os.path.exists(out_folder):
            os.makedirs(out_folder, exist_ok=True)
        with open(os.path.join(out_folder, "ensemble_summary.txt"), 'w') as f:
            f.writelines("Summary of Matched PWMs: \n")
            for ind,pwm in enumerate(pwms):
                f.writelines("\nInput PWM: ")
                f.writelines(str(pwm).split('/')[-1] + "\n")

                new_out_folder=os.path.join(out_folder,str(ind))
                if not os.path.exists(new_out_folder):
                    os.makedirs(new_out_folder, exist_ok=True)
                temp_item_list, temp_list_score = searching(pwm=pwm,out=new_out_folder,tf_name=tf_name,db_path=db_folder,db_type=db_type,
                         input_count=input_count,
                         db_count=db_count,
                         db_file_type=db_file_type,
                         input_file_type=input_file_type,
                         input_prop=input_prob,n=n,
                         db_prob=db_prob)
                f.writelines("\nMatched PWMs : \n")
                for file in temp_item_list:
                    f.writelines(str(file)+"\n")
                all_matched_list.append(temp_item_list)
                all_score_list.append(temp_list_score)

            flat_list_items = [item for sublist in all_matched_list for item in sublist]
            flat_list_scores = [item for sublist in all_score_list for item in sublist]

            f.writelines(str(Counter(flat_list_items))+"\n")
            f.writelines("Median: " + str(np.median(flat_list_scores)))


            print(Counter(flat_list_items))
            # print(flat_list_scores)
            print("Median: " , np.median(flat_list_scores))
            print("See summary of ensemble and search in ", os.path.join(out_folder,"ensemble_summary.txt"))

if __name__ == '__main__':
    src_folder = '../data/out/ensemble/swi4/'
    out= '../data/out/ensemble/search_out/'
    db_folder = '../data/in/SGD/SGD_yeast_pwm/'
    # dst_folder = '../data/out/ensemble/in/'
    # prepare_for_clustering(src_folder, dst_folder)
    # print("step 1")
    # do_clustering(dst_folder, '../data/out/ensemble/clustering_out/')
    # print("step 2")
    #
    # make_rep('../data/out/ensemble/clustering_out/out/', dbd='selected', clusters='all', ic=0.4,
    #          best_match_initial_motif=1, mean_threshold=0.75, z_score_threshold=-0.9,
    #          top_occurrence=0.35, occurrence_threshold=0.25)
    #
    # ensembling_compare_investigate(path_of_clusters='../data/out/ensemble/clustering_out/out/',
    #                                predicted='../data/out/ensemble/predicted')
    # print("step 3")
    #
    # do_searching(in_folder='../data/out/ensemble/predicted/',
    #              db_folder='../data/in/SGD/SGD_yeast_pwm/',
    #              out_folder='../data/out/search_out/',
    #              db_type='folder',
    #              input_count=1,
    #              n=2,
    #              db_count=0,
    #              db_file_type='.txt',
    #              input_file_type='.mlp',
    #              input_prob=0,
    #              db_prob=1,
    #              tf_name='')
    PredictedCompare(path_to_predicted_files=src_folder,
                     output_folder=out,
                     db_folder=db_folder, tf_name='')