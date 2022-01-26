
import os, sys
from abc4pwm.energy_to_p import read_energy_matrix
from abc4pwm.convert_count_to_pwm import motif_weight2p
from abc4pwm.similarity_score import compute_similarity_score4alignment
import numpy as np
import shutil
from scipy import stats
from statistics import median
from glob import glob
from pathlib import Path
from collections import Counter
from itertools import takewhile
from distutils.dir_util import copy_tree
import pandas as pd


class ClusterQuality():
    def __init__(self, input_path_to_folder_of_DBDs, out_path_for_qa_clusters, output_folder_for_text_report,
                 output_path_for_quality_assessment_file='default',
                 load_new_assesment=0,
                 minimum_pwms_in_dbd=5,
                 mean_threshold=0.77,
                 z_score_threshold=-1.2,
                 top_occurrence=0.05,
                 occurrence_threshold=0.05):
        """
        :param input_path_to_folder_of_DBDs: path to the folder which contained clustered DBD
        ;param out_path_for_qa_clusters: path where quality assessment clusters should go.
        :param output_folder_for_text_report: path to folder where you want your report
        ;param load_new_assessment: write 1 if new assessment is needed, default is 0.
        ;param output_path_for_quality_assessment_file: folder where quality assessment .json file will be saved
        :param mean_threshold: minimum mean of similarities to qualify as good cluster
        :param z_score_threshold: minimum negative z score for similarity of a pwm from mean
        :param top_occurrence: top occurrences of pwms which qualify z score criteria
        :param occurrence_threshold: pick occurrences threshold to declare unknown/wrongly clustered pwms

        """
        input_path_to_folder_of_DBDs = os.path.join(input_path_to_folder_of_DBDs,'')
        out_path_for_qa_clusters = os.path.join(out_path_for_qa_clusters,'')
        output_folder_for_text_report = os.path.join(output_folder_for_text_report,'')
        output_path_for_quality_assessment_file = os.path.join(output_path_for_quality_assessment_file,'')


        self.minimum_pwms_in_dbd = minimum_pwms_in_dbd
        self.mean_threshold = mean_threshold
        self.z_score_threshold = z_score_threshold
        self.top_occurrence = top_occurrence
        self.occurrence_threshold = occurrence_threshold
        self.total_clusters = 0
        self.output_folder_for_text_report = output_folder_for_text_report


        if 'default' in output_path_for_quality_assessment_file:
            leaf_folder = Path(self.output_folder_for_text_report)
            out_dir = leaf_folder.parent
            data_dir = out_dir.parent
            path_for_quality_assessment_file = os.path.join(data_dir,'in/')
        else:
            path_for_quality_assessment_file = output_path_for_quality_assessment_file

        if not os.path.exists(output_folder_for_text_report):
            os.makedirs(output_folder_for_text_report, exist_ok=True)
        if not os.path.exists(out_path_for_qa_clusters):
            os.makedirs(out_path_for_qa_clusters)
        self.empty_dir(out_path_for_qa_clusters)

        copy_tree(input_path_to_folder_of_DBDs, out_path_for_qa_clusters)
        input_path_to_folder_of_DBDs = out_path_for_qa_clusters

        if load_new_assesment:
            print("Preparing summary of clusters...")
            quality_df = self.cluster_quality(self, input_path_to_folder_of_DBDs, output_folder_for_text_report, path_for_quality_assessment_file)
        else:
            if os.path.exists(os.path.join(path_for_quality_assessment_file,'quality_df.json')):
                quality_df = pd.read_json(os.path.join(output_path_for_quality_assessment_file,'quality_df.json'))


            else:
                print("Quality assesment file is not available. Please generate assesment by setting load_new_assesment = True")
                exit()




        print("Now you can see clustering summary and quality in ", self.output_folder_for_text_report, " in cluster_quality.txt \n"
            "Qaulity assessment file is saved in .json formate in ", path_for_quality_assessment_file)

    @staticmethod
    def get_items_upto_count(dct, n,top_percentage):
      data = dct.most_common(top_percentage)
      #Now collect all items whose value is greater than `n`.
      return list(takewhile(lambda x: x[1] > n, data))

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

    @staticmethod
    def cluster_quality(self, input_folder,output_folder,path_for_quality_assessment_file ):
        # calculates quality of any given cluster by different measures.
        # also do some statistical analysis like median, mean etc

        cols = ['dbd','cluster_number', 'mean', 'median','min','max', 'std','total_pwms','uncertain_pwms', 'similarities']
        quality_df_data = []
        leaf_folder = Path(output_folder)
        out_dir = leaf_folder.parent
        path_to_text_reports = output_folder
        dst_for_bad_pwms = os.path.join(out_dir,'uncertain_pwms/')
        if not os.path.exists(dst_for_bad_pwms):
            os.makedirs(dst_for_bad_pwms,exist_ok=True)
        self.empty_dir(dst_for_bad_pwms)
        if os.path.exists(os.path.join(path_to_text_reports, 'cluster_quality.txt')):
            os.remove(os.path.join(path_to_text_reports, 'cluster_quality.txt'))
        if not os.path.exists(path_to_text_reports):
            os.makedirs(path_to_text_reports, exist_ok=True)
        with open(os.path.join(path_to_text_reports, 'cluster_quality.txt'), 'w') as f:
            total_number_of_clusters = 0
            total_bad_pwms = 0
            cluster_size_list = []
            similarityMatrices = []

            dbds = sorted(os.listdir(input_folder))
            for ind, i in enumerate(dbds):
                if i.startswith('.DS'):
                    dbds.pop(ind)

            for dbd in dbds:
                dst_for_bad_pwms_n = os.path.join(dst_for_bad_pwms,dbd,'')
                if not os.path.exists(dst_for_bad_pwms_n):
                    os.makedirs(dst_for_bad_pwms_n, exist_ok=True)
                if not os.path.exists(os.path.join(input_folder , dbd,'out/')):
                    # pwms_full = glob(os.path.join(input_folder,dbd,"/*.mlp"))
                    pwms_full = os.listdir(os.path.join(input_folder,dbd))
                    total_bad_pwms += len(pwms_full)
                    for x in pwms_full:
                        shutil.move(os.path.join(input_folder,dbd,x) , dst_for_bad_pwms)
                    continue

                if 'Unknown' in dbd:
                    for i in os.listdir(os.path.join(input_folder , dbd,'out/')):
                        unknown_pwms = glob(str(os.path.join(input_folder , dbd,'out/'))+'*/*.mlp')
                        for pwm in unknown_pwms:
                            shutil.move(pwm,dst_for_bad_pwms_n)
                    shutil.rmtree(os.path.join(input_folder , dbd))
                    continue
                cluster_path = os.path.join(input_folder, dbd, 'out/')
                f.writelines(str(dbd) + '\n')
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

                clusters_list = [i for i in sorted(os.listdir(cluster_path)) if os.path.isdir(os.path.join(cluster_path,i))]
                clusters_list = [int(x) for x in clusters_list]
                clusters_list.sort()


                total_number_of_clusters += len(clusters_list)
                total_pwms_in_dbd = 0
                for cluster in clusters_list:

                    pwms_full = glob(os.path.join(input_folder,dbd,'out', str(cluster)+ "/*.mlp"))
                    pwms_full = sorted(pwms_full)

                    cluster_size_list.append(len(pwms_full))
                    total_pwms_in_dbd = len(pwms_full) + total_pwms_in_dbd
                    similarityMatrix = np.zeros((len(pwms_full), len(pwms_full)))

                    for index1, i in enumerate(pwms_full):
                        matrix1, matrix_string1, maximum_feq1, total_maximum1, info1 = read_energy_matrix(i)
                        normalized_matrix1 = motif_weight2p(matrix1)
                        for index2, j in enumerate(pwms_full):
                            matrix2, matrix_string2, maximum_feq2, total_maximum2, info2 = read_energy_matrix(j)

                            normalized_matrix2 = motif_weight2p(matrix2)
                            similarityMatrix[index1, index2] = compute_similarity_score4alignment(normalized_matrix1, normalized_matrix2)

                    uper_triangle_similarity_matrix = similarityMatrix[np.triu_indices(similarityMatrix.shape[0], k=1)]
                    uper_triangle_similarity_matrix_indeces = np.triu_indices(similarityMatrix.shape[0], k=1)

                    minimum_similarity_pwm_index = np.unravel_index(np.argmin(sum(similarityMatrix), axis=None), sum(similarityMatrix).shape)


                    if len(pwms_full) > 1:
                        if uper_triangle_similarity_matrix.mean() > self.mean_threshold:
                            potential_bad_pwms = []
                            if uper_triangle_similarity_matrix.std() > 0.0:
                                for index, k in enumerate(stats.zscore(uper_triangle_similarity_matrix)):
                                    if k < self.z_score_threshold:
                                        potential_bad_pwms.append(uper_triangle_similarity_matrix_indeces[0][index])
                                        potential_bad_pwms.append(uper_triangle_similarity_matrix_indeces[1][index])

                                bad_pwms_indexes = self.get_items_upto_count(Counter(potential_bad_pwms),
                                                                        int(np.ceil(len(pwms_full) * self.occurrence_threshold)),
                                                                        int(np.ceil(len(pwms_full) * self.top_occurrence)))

                                del_from_pwms_list = [i[0] for i in bad_pwms_indexes]
                                bad_pwms = [pwms_full[x] for x in del_from_pwms_list]



                                f.writelines(str(cluster) +
                                      "\t" + str(uper_triangle_similarity_matrix.mean()) +
                                      "\t" + str(median(uper_triangle_similarity_matrix)) +
                                      "\t" + str(min(uper_triangle_similarity_matrix)) + " " +
                                      str(minimum_similarity_pwm_index).split(',')[0].split('(')[1] +
                                      "\t" + str(uper_triangle_similarity_matrix.std()) +
                                      "\t" + str(len(pwms_full)) +
                                      "\t" + str(len(bad_pwms)) +
                                      "\t" + str(int(((len(pwms_full)-len(bad_pwms))/len(pwms_full))*100)) )
                                if len(bad_pwms) > 0:
                                    f.writelines("%\t" + str([i[0] for i in bad_pwms_indexes]) + "\n")
                                    total_bad_pwms += len(bad_pwms)
                                    for x in bad_pwms:
                                        # pwms.remove(x)
                                        shutil.move(x, dst_for_bad_pwms_n)
                                else:
                                    f.writelines("\n")
                        else:
                            f.writelines(str(cluster) +
                                         "\t" + str(uper_triangle_similarity_matrix.mean()) +
                                         "\t" + str(median(uper_triangle_similarity_matrix)) +
                                         "\t" + str(min(uper_triangle_similarity_matrix)) + " " +
                                         str(minimum_similarity_pwm_index).split(',')[0].split('(')[1] +
                                         "\t" + str(uper_triangle_similarity_matrix.std()) +
                                         "\t" + str(len(pwms_full)) +
                                         "\t Mean less than "+ str(self.mean_threshold)+ ", bad cluster " +
                                         "\n")
                        try:
                            quality_df_data.append([dbd,
                                                   cluster,
                                                   uper_triangle_similarity_matrix.mean(),
                                                   median(uper_triangle_similarity_matrix),
                                                   min(uper_triangle_similarity_matrix),
                                                   max(uper_triangle_similarity_matrix),
                                                   uper_triangle_similarity_matrix.std(),
                                                   len(pwms_full),
                                                   len(bad_pwms), uper_triangle_similarity_matrix])
                        except:
                            quality_df_data.append([dbd,
                                                    cluster,
                                                    uper_triangle_similarity_matrix.mean(),
                                                    median(uper_triangle_similarity_matrix),
                                                    min(uper_triangle_similarity_matrix),
                                                    max(uper_triangle_similarity_matrix),
                                                    uper_triangle_similarity_matrix.std(),
                                                    len(pwms_full),
                                                    0, uper_triangle_similarity_matrix])
                    else:
                        f.writelines(str(cluster) +
                                  "\t\t\t Single pwm cluster, no need of analysis" +
                                  "\n")



                f.writelines("\nTotal PWMs : " + str(total_pwms_in_dbd)+
                             "\t Total Clusters in DBD: " + str(len(clusters_list)) +
                             "\n\n\n")
            f.writelines("Total DBDs : " + str(len(sorted(os.listdir(input_folder)))) +
                         "\t Total Number of clusters " + str(total_number_of_clusters ) +
                         "\t Total bad PWMs " + str(total_bad_pwms))
            f.writelines("\n"
                         "Clusters size summary : " + str(Counter(cluster_size_list)))

            quality_df = pd.DataFrame(data=quality_df_data, columns=cols)
            leaf_folder = Path(out_dir)
            data_path = leaf_folder.parent
            if not os.path.exists(path_for_quality_assessment_file):
                os.makedirs(path_for_quality_assessment_file, exist_ok=True)
            quality_df.to_json(os.path.join(path_for_quality_assessment_file,'quality_df.json'))
            # print(quality_df['dbd'])
            return quality_df



if __name__ == "__main__":

    # if len(sys.argv) < 2:
    #     "Usage python quality_assessment.py <path_to_folder_of_DBDs>"
    #     exit()
    # clusterQualityobj = ClusterQuality(input_path_to_folder_of_DBDs='../data/out/clustering_out/',
    #                                    out_path_for_qa_clusters='../data/out/quality_assessed_out',
    #                                    output_folder_for_text_report='../data/out/reports_in_text/',
    #                                    mean_threshold=0.80,
    #                                    occurrence_threshold=0.05,
    #                                    z_score_threshold=-1.2,
    #                                    top_occurrence=0.05,
    #                                    load_new_assesment=1)

    # abc4pwm_new qualty assement for comparison
    clusterQualityobj = ClusterQuality(
        input_path_to_folder_of_DBDs='/Users/omerali/Desktop/PhD/Publication/Next/Clustering/Comparison/in_abc4pwm_motif_clustering/',
        out_path_for_qa_clusters='/Users/omerali/Desktop/PhD/Publication/Next/Clustering/Comparison/quality_assessed_out_abc4pwm_new',
        output_folder_for_text_report='/Users/omerali/Desktop/PhD/Publication/Next/Clustering/Comparison/reports_in_text/',
        mean_threshold=0.80,
        occurrence_threshold=0.05,
        z_score_threshold=-1.2,
        top_occurrence=0.05,
        load_new_assesment=1)

    # stamp qualty assement for comparison

    # clusterQualityobj = ClusterQuality(
    #     input_path_to_folder_of_DBDs='/Users/omerali/Desktop/PhD/Publication/Next/Clustering/Comparison/stamp_in_abc4pwm_motif/',
    #     out_path_for_qa_clusters='/Users/omerali/Desktop/PhD/Publication/Next/Clustering/Comparison/quality_assessed_out_stamp',
    #     output_folder_for_text_report='/Users/omerali/Desktop/PhD/Publication/Next/Clustering/Comparison/reports_in_text/',
    #     mean_threshold=0.65,
    #     occurrence_threshold=0.05,
    #     z_score_threshold=-1.2,
    #     top_occurrence=0.05,
    #     load_new_assesment=1)


