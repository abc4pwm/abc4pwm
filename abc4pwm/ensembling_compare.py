import os
import shutil
from glob import glob
from datetime import datetime
import numpy as np
from scipy import stats
from statistics import median
from collections import Counter
from pathlib import Path
from abc4pwm.non_dbd_clustering import non_dbd_ClusteringPwm
from abc4pwm.searching import searching
from abc4pwm.pwm_representative_test import make_representative_pwm
from abc4pwm.convert_count_to_pwm import motif_weight2p
from abc4pwm.energy_to_p import read_energy_matrix
from abc4pwm.similarity_score import compute_similarity_score4alignment
from itertools import takewhile
from distutils.dir_util import copy_tree


class PredictedCompare():
    def __init__(self, path_to_predicted_files, output_folder, db_folder, dst_for_bad_pwms, min_pwms_in_cluster=3, mean_threshold=0.75, z_score_threshold=-0.9,
                 top_occurrence=0.35, occurrences_threshold=0.25, ic_for_rep=0.4,top_n=2, db_type='folder',input_count=1,db_count=0, db_file_type='.txt',
                     input_file_type='.mlp', input_prob=0, db_prob=1, tf_name='', qa = 1, seed = 0, damp=0.5, max_iter=400, convergence_iter=30, preference=None):
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
        ;param qa: quality assesment paramter. qa = 1 for doing quality assesment, otherwise zero.
        :param mean_threshold: minimum mean of similarities to qualify as good cluster
        :param z_score_threshold: minimum negative z score for similarity of a pwm from mean
        :param top_occurrence: top occurrences of pwms which qualify z score criteria
        :param occurrence_threshold: pick occurrences threshold to declare unknown/wrongly clustered pwms
        :param damp: Damping factor (between 0.5 and 1) is the extent to which the current value is maintained relative to incoming values (weighted 1 - damping).
        This in order to avoid numerical oscillations when updating these values (messages).
        :param max_iter: default=200, Maximum number of iterations.
        :param convergence_iter: default=15 Number of iterations with no change in the number of estimated clusters that stops the convergence.
        :param preference: default=None Preferences for each point - points with larger values of preferences are more likely to be chosen as exemplars.
        The number of exemplars, ie of clusters, is influenced by the input preferences value. If the preferences are not passed as arguments, they will be set to the median of the input similarities.

        """


        output_folder = os.path.join(output_folder,'')
        path_to_predicted_files = os.path.join(path_to_predicted_files,'')
        db_folder = os.path.join(db_folder,'')
        temp2_folder = Path(output_folder)
        # temp2_folder = temp_folder.parent
        path_to_text_report = os.path.join(temp2_folder,'reports_in_text/')
        if not os.path.exists(path_to_text_report):
            os.makedirs(path_to_text_report, exist_ok=True)



        with open(os.path.join(path_to_text_report,'parameters.txt'),  'w' ) as fp:
            now = datetime.now()

            current_time = now.strftime("%H:%M:%S")
            fp.writelines("Following parameters are used for insemble investigate run at " + str(current_time) + '\n\n')
            fp.writelines(str(path_to_predicted_files) + " : Path to predicted files \n")
            fp.writelines(str(output_folder) + " : Output folder \n")
            fp.writelines(str(db_folder) + " : Database folders \n")
            fp.writelines(str(dst_for_bad_pwms) + " : Bad pwms folder \n")
            fp.writelines(str(min_pwms_in_cluster) + " : Minimum Pwms in a cluster for results of ensembling learning and quality cluster \n")

            fp.writelines(str(z_score_threshold) + " : Z_score \n")
            fp.writelines(str(mean_threshold) + " : Mean\n")
            fp.writelines(str(top_occurrence) + " : Top Occurances\n")
            fp.writelines(str(occurrences_threshold) + " : Occurances Threshold\n")

            fp.writelines(str(ic_for_rep) + " : ic for rep \n")
            fp.writelines(str(top_n) + " : Number of top mathces for search results \n")
            fp.writelines(str(db_type) + " : Type of database \n")
            fp.writelines(str(input_count) + " : Input is count data? \n")
            fp.writelines(str(db_count) + " : Database is count data? \n")

            fp.writelines(str(seed) + " : Seed \n")
            fp.writelines(str(damp) + " : Damp \n")
            fp.writelines(str(max_iter) + " : Maximum iterations for Affinity Propagation \n")
            fp.writelines(str(convergence_iter) + " : Convergence Iterations \n")
            fp.writelines(str(preference) + " : Preference for Affinity Propagation \n")



        dst_folder = os.path.join(temp2_folder, 'in/')
        out_path_for_qa_clusters = os.path.join(temp2_folder, 'quality_assessed_out/')
        if not os.path.exists(out_path_for_qa_clusters):
            os.makedirs(out_path_for_qa_clusters, exist_ok=True)
        self.empty_dir(out_path_for_qa_clusters)
        if not os.path.exists(dst_for_bad_pwms):
            os.makedirs(dst_for_bad_pwms, exist_ok=True)
        if not os.path.exists(dst_folder):
            os.makedirs(dst_folder, exist_ok=True)
        self.empty_dir(dst_folder)
        self.empty_dir(dst_for_bad_pwms)
        self.prepare_for_clustering(path_to_predicted_files, dst_folder)
        print("step 1: Prepare for clustering Done")
        cluster_out_folder = os.path.join(temp2_folder,'clustering_out/')
        #jbw 2024
        self.do_clustering(dst_folder, cluster_out_folder,path_to_text_report,  seed,  damp, max_iter, convergence_iter, preference)
        print("step 2: Clustering Done")

        if qa:
            print("Going through qa now")
            self.quality_assessment(cluster_out_folder,
                                    out_path_for_qa_clusters,
                                    dst_for_bad_pwms,
                                    path_to_text_reports=path_to_text_report,
                                    mean_threshold=mean_threshold,
                                    z_score_threshold=z_score_threshold,
                                    top_occurrence=top_occurrence,
                                    occurrence_threshold=occurrences_threshold)

        #jbw 2024
        tmp_path_for_clusters=os.path.join(out_path_for_qa_clusters, 'out/')
        if not os.path.exists(tmp_path_for_clusters):
              os.makedirs(tmp_path_for_clusters)

        #self.make_rep(os.path.join(out_path_for_qa_clusters, 'out/'), dbd='selected', clusters='all', ic=ic_for_rep,
        #         best_match_initial_motif=1, mean_threshold=mean_threshold, z_score_threshold=z_score_threshold,
        #         top_occurrence=top_occurrence, occurrences_threshold=occurrences_threshold)

        self.make_rep(tmp_path_for_clusters, dbd='selected', clusters='all', ic=ic_for_rep,
                 best_match_initial_motif=1, mean_threshold=mean_threshold, z_score_threshold=z_score_threshold,
                 top_occurrence=top_occurrence, occurrences_threshold=occurrences_threshold)
        #end jbw

        predicted_folder = os.path.join(temp2_folder,'predicted')
        self.ensembling_compare_investigate(os.path.join(out_path_for_qa_clusters, 'out/'),
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
                     tf_name=tf_name,
                     path_to_text_reports=path_to_text_report,
                      mean_threshold=mean_threshold,
                      z_score_threshold=z_score_threshold,
                      top_occurrence=top_occurrence,
                      occurrence_threshold=occurrences_threshold,
                      min_threshold_for_number_of_pwms=min_pwms_in_cluster,
                          seed=seed)
        print("step 4: Searching Done")


    def prepare_for_clustering(self,src, dst):
        if not os.path.exists(dst):
            os.makedirs(dst, exist_ok=True)
        self.empty_dir(dst)
        pwms = glob(src+"/*.mlp")
        pwms = sorted(pwms)
        for pwm in pwms:
            shutil.copy(pwm,dst)
    #jbw 2024
    def do_clustering(self, in_folder, out_folder, path_to_text_report, seed, damp, max_iter, convergence_iter, preference):
        non_dbd_ClusteringPwm(in_folder,out_folder, path_to_text_report, seed, damp, max_iter, convergence_iter, preference)

    def make_rep(self, path_to_clusters,  dbd = 'selected', clusters = 'all' , ic = 0.4 , best_match_initial_motif = 1, mean_threshold= 0.75, z_score_threshold = -0.9,
                     top_occurrence = 0.35, occurrences_threshold=0.25):
        make_representative_pwm(path_to_clusters,  dbd, clusters, ic , best_match_initial_motif, mean_threshold, z_score_threshold ,
                     top_occurrence, occurrences_threshold)

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
            pwms = glob(os.path.join(path_of_clusters,str(cluster))+"/*.mlp")
            pwms = sorted(pwms)
            if len(pwms) < int(min_threshold_for_number_of_pwms):
                if not os.path.exists(predicted):
                    os.makedirs(predicted)

            else:
                src = os.path.join(path_of_clusters,str(cluster),"repres",str(cluster)+'_rep.mlp')
                dst = os.path.join(predicted,pwms[0].split('/')[-1]+'_'+str(cluster)+'_rep.mlp')
                shutil.copy(src,dst)


    def do_searching(self, in_folder,db_folder, out_folder, tf_name ,n,db_type ,input_count,db_count, db_file_type,
                     input_file_type, input_prob, db_prob, path_to_text_reports, mean_threshold,
                                    z_score_threshold, top_occurrence, occurrence_threshold,min_threshold_for_number_of_pwms , seed):
        pwms = glob(in_folder+"/*.mlp")
        pwms = sorted(pwms)
        all_matched_list = []
        all_score_list = []
        if not os.path.exists(out_folder):
            os.makedirs(out_folder, exist_ok=True)

        with open(os.path.join(path_to_text_reports, "ensemble_summary.txt"), 'w') as f:
            f.writelines("Quality assessment is done for parameters provided below: \n\n")
            f.writelines(str(min_threshold_for_number_of_pwms) + " : Minimum threshold of PWMs in a cluster \n")
            f.writelines(str(n) + " : Top N matches\n")
            f.writelines(str(in_folder) + " : Input Folder for search..\n")
            f.writelines(str(out_folder) + " : Out folder for search..\n")
            f.writelines(str(db_folder) + " : Database folder for search results..\n")

            f.writelines(str(seed) + " : Seed \n")
            f.writelines(str(z_score_threshold) + " : Z_score \n")
            f.writelines(str(mean_threshold) + " : Mean\n")
            f.writelines(str(top_occurrence) + " : Top Occurances\n")
            f.writelines(str(occurrence_threshold) + " : Occurances Threshold\n\n")
            f.writelines("Summary of Matched PWMs: \n")
            for ind,pwm in enumerate(pwms):
                f.writelines("\nSearch out folder " +str(ind) )
                f.writelines("\nInput PWM: ")
                f.writelines(str(pwm).split('/')[-1] + "\n")

                new_out_folder=os.path.join(out_folder,str(ind))
                if not os.path.exists(new_out_folder):
                    os.makedirs(new_out_folder, exist_ok=True)
                temp_item_list, temp_list_score = searching(pwm=pwm, out=new_out_folder, tf_name=tf_name, db_path=db_folder, db_type=db_type,
                         input_count=input_count,
                         db_count=db_count,
                         db_file_type=db_file_type,
                         input_file_type=input_file_type,
                         input_prop=input_prob, n=n,
                         db_prob=db_prob)

                f.writelines("\nMatched PWMs : \n")
                for ind, file in enumerate(temp_item_list):
                    f.writelines(str(file)+"\t"+ str(temp_list_score[ind]) +"\n")
                all_matched_list.append(temp_item_list)
                all_score_list.append(temp_list_score)

            flat_list_items = [item for sublist in all_matched_list for item in sublist]
            flat_list_scores = [item for sublist in all_score_list for item in sublist]

            f.writelines(str(Counter(flat_list_items))+"\n")
            f.writelines("Median: " + str(np.median(flat_list_scores)))


            print(Counter(flat_list_items))
            # print(flat_list_scores)
            print("Median: " , np.median(flat_list_scores))
            print("See summary of ensemble and search in ", os.path.join(path_to_text_reports,"ensemble_summary.txt"))

    def get_items_upto_count(self, dct, n,top_percentage):
      data = dct.most_common(top_percentage)
      #Now collect all items whose value is greater than `n`.
      return list(takewhile(lambda x: x[1] > n, data))


    def quality_assessment(self, cluster_path, out_path_for_qa_clusters ,dst_for_bad_pwms, path_to_text_reports , mean_threshold, z_score_threshold, occurrence_threshold, top_occurrence):

        copy_tree(cluster_path, out_path_for_qa_clusters)

        out_path_for_qa_clusters = os.path.join(out_path_for_qa_clusters,'out')
        clusters_list = [i for i in sorted(os.listdir(out_path_for_qa_clusters)) if os.path.isdir(os.path.join(out_path_for_qa_clusters, i))]
        clusters_list = [int(x) for x in clusters_list]
        clusters_list.sort()
        if not os.path.exists(path_to_text_reports):
            os.makedirs(path_to_text_reports, exist_ok=True)
        with open(os.path.join(path_to_text_reports, 'ensemble_cluster_quality.txt'), 'w') as f:




            f.writelines(
                'Cluster Number '
                '\t Mean '
                '\t Median '
                '\t Minimum '
                '\t Standard Deviation '
                '\t Total No. of PWMs '
                '\t No. of Wrong PWMs '
                '\t Percentage of Quality '
                '\t Wrongly clustered PWMS \n')
            for cluster in clusters_list:

                pwms_full = glob(os.path.join(out_path_for_qa_clusters, str(cluster) + "/*.mlp"))
                pwms_full = sorted(pwms_full)


                similarityMatrix = np.zeros((len(pwms_full), len(pwms_full)))

                for index1, i in enumerate(pwms_full):
                    matrix1, matrix_string1, maximum_feq1, total_maximum1, info1 = read_energy_matrix(i)
                    normalized_matrix1 = motif_weight2p(matrix1)
                    for index2, j in enumerate(pwms_full):
                        matrix2, matrix_string2, maximum_feq2, total_maximum2, info2 = read_energy_matrix(j)

                        normalized_matrix2 = motif_weight2p(matrix2)
                        similarityMatrix[index1, index2] = compute_similarity_score4alignment(normalized_matrix1,
                                                                                              normalized_matrix2)

                uper_triangle_similarity_matrix = similarityMatrix[np.triu_indices(similarityMatrix.shape[0], k=1)]
                uper_triangle_similarity_matrix_indeces = np.triu_indices(similarityMatrix.shape[0], k=1)

                # minimum_similarity_pwm_index = np.unravel_index(np.argmin(sum(similarityMatrix), axis=None),
                #                                                 sum(similarityMatrix).shape)

                if len(pwms_full) > 1:
                    if uper_triangle_similarity_matrix.mean() > mean_threshold:
                        potential_bad_pwms = []
                        for index, k in enumerate(stats.zscore(uper_triangle_similarity_matrix)):
                            if k < z_score_threshold:
                                potential_bad_pwms.append(uper_triangle_similarity_matrix_indeces[0][index])
                                potential_bad_pwms.append(uper_triangle_similarity_matrix_indeces[1][index])

                        bad_pwms_indexes = self.get_items_upto_count(Counter(potential_bad_pwms),
                                                                     int(np.ceil(len(
                                                                         pwms_full) * occurrence_threshold)),
                                                                     int(np.ceil(
                                                                         len(pwms_full) * top_occurrence)))

                        del_from_pwms_list = [i[0] for i in bad_pwms_indexes]
                        bad_pwms = [pwms_full[x] for x in del_from_pwms_list]


                        f.writelines(str(cluster) +
                                     "\t" + str(uper_triangle_similarity_matrix.mean()) +
                                     "\t" + str(median(uper_triangle_similarity_matrix)) +
                                     "\t" + str(min(uper_triangle_similarity_matrix)) + " " +
                                     "\t" + str(uper_triangle_similarity_matrix.std()) +
                                     "\t" + str(len(pwms_full)) +
                                     "\t" + str(len(bad_pwms)) +
                                     "\t" + str(int(((len(pwms_full) - len(bad_pwms)) / len(pwms_full)) * 100)))
                        if len(bad_pwms) > 0:
                            f.writelines("%\t" + str([i[0] for i in bad_pwms_indexes]) + "\n")
                            # total_bad_pwms += len(bad_pwms)
                            for x in bad_pwms:
                                # pwms.remove(x)
                                shutil.move(x, dst_for_bad_pwms)
                        else:
                            f.writelines("\n")
                            #print("Nothing")



if __name__ == '__main__':


    src_folder = '../data/out/ensemble/time_series/t0/0strand/'
    out = '../data/out/ensemble/time_series/t0/0strand/search_out/'
    dest_4_bd_pwms = '../data/out/ensemble/time_series/t0/0strand/bad_pwms/'
    db_folder = '../data/in/in_pwms/'
    db_file_type = '.mlp'
    db_count = 1
    db_prob = 0
    input_count = 1
    input_file_type = '.mlp'
    input_prob = 0
    min_pwms_in_cluster= 5
    top_n=2
    seed= 1
    ic_for_rep = 0
    mean_threshold=0.5
    z_score_threshold=0.0
    damp=0.7
    max_iter=1000
    convergence_iter=30
    top_occurrences=0.5
    occurrences_threshold=0.5
    db_type = 'folder'



    PredictedCompare(path_to_predicted_files=src_folder,
                     dst_for_bad_pwms=dest_4_bd_pwms,
                     output_folder=out,
                     db_folder=db_folder,
                     db_file_type = db_file_type,
                     db_count=db_count,
                     db_prob=db_prob,
                     tf_name='',
                     qa = 1,
                     seed=seed ,
                     ic_for_rep=ic_for_rep,
                     min_pwms_in_cluster=min_pwms_in_cluster, mean_threshold=mean_threshold, z_score_threshold=z_score_threshold,
                     top_occurrence=top_occurrences,
                     occurrences_threshold=occurrences_threshold,
                     top_n=top_n,
                     db_type=db_type,
                     input_count=input_count,
                     input_prob=input_prob,
                     input_file_type=input_file_type,
                     damp=damp,
                     max_iter=max_iter,
                     convergence_iter=convergence_iter,
                     )
