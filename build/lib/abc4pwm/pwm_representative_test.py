import os, sys
import numpy as np
from abc4pwm.pwm_representative import PWMRepresentative
from abc4pwm.similarity_score import compute_similarity_score4alignment
from scipy import stats
from glob import glob
from collections import Counter
from itertools import takewhile
from abc4pwm.energy_to_p import read_energy_matrix
from pathlib import Path
import shutil

#this is a modified version of omer's test script , which is using wang's new functions of similarity score and represenative motif alignment
#June 2020 Wang


def motif_weight2p(motif_energy):
    w = np.transpose(motif_energy)
    new_w = (w + 0.25) / (w.sum(axis=0) + 1)
    psam = new_w / np.sum(new_w, axis=0)
    return psam

def get_char(num):
    if num == 0:
        return 'a'
    elif num == 1:
        return 'c'
    elif num == 2:
        return 'g'
    elif num == 3:
        return 't'
    return ''

def trim_edges_of_motif(matrix_1, cut_off = 1.8):

    temp = matrix_1.copy()
    for i in range(2):
        rows_to_delete = []
        start = True
        for ind, line in enumerate(temp):
            t2 = [x for x in line if x > cut_off]
            if not len(t2)>0 and start:
                rows_to_delete.append(ind)
            else:
                break
        temp = np.delete(temp, rows_to_delete, axis=0)
        temp = temp[::-1]
    return temp

def get_items_upto_count(dct, n,top_percentage):
  data = dct.most_common(top_percentage)
  #Now collect all items whose value is greater than `n`.
  return list(takewhile(lambda x: x[1] > n, data))


def read_matrix_file(path):
    with open(path) as bb:
        lines = bb.readlines()
    lst = list()
    for line in lines:
        if '#' in line:
            continue
        splited = line.rstrip().split('\t')[1:]
        lst.append([float(x) for x in splited])
    return np.asarray(lst)

def make_representative_pwm(path_to_clusters,  dbd = 'selected', clusters = 'all' , ic = 0.4 , best_match_initial_motif = 1, mean_threshold= 0.75, z_score_threshold = -0.9,
                 top_occurrence = 0.35, occurrence_threshold=0.25):
    """

    :param path_to_clusters: this path should be path to clusters for which you want to make a representative motif
    :param dbd: dbd for which representative motif is to be calculated. If 'all' then 'path_to_clusters' should be path
    to folder which contain clustered dbds
    :param clusters: cluster numbers which you want to plot, For example, 0,1,2,3,4 if you want to plot these clusters
    ;param ic: information content cutt off for triming edges. Default is 0.4
    :param best_match_initial_motif: True if you want best match inital motif , False if random motif
    :param mean_threshold: minimum mean of similarities to qualify as good cluster
    :param z_score_threshold: minimum negative z score for similarity of a pwm from mean
    :param top_occurrence: top occurences of pwms which qualify z score criteria
    :param occurrence_threshold: pick occurences threshold to declare unknown/wrongly clustered pwms

    :return:
    """
    print('\nPreparing representative motifs...')

    if 'all' in dbd or '*' in dbd:
        dbds = [i for i in sorted(os.listdir(path_to_clusters)) if not i.endswith('.DS_Store')]
        for dbd in dbds:
            if os.path.exists(os.path.join(path_to_clusters, dbd, 'out/')):
                clusters = [i for i in sorted(os.listdir(os.path.join(path_to_clusters, dbd, 'out/'))) if
                        not i.endswith('.DS_Store')]
                representaive_for_clusters(os.path.join(path_to_clusters, dbd, 'out/'),clusters, ic, best_match_initial_motif, mean_threshold, z_score_threshold,
                 top_occurrence, occurrence_threshold)
        exit()

    else:
        if 'all' in clusters or '*' in clusters:
            clusters = [i for i in sorted(os.listdir(path_to_clusters)) if not i.endswith('.DS_Store')]

        else:
            clusters = clusters.split(',')


        for ind, i in enumerate(clusters):
            if not os.path.isdir(os.path.join(path_to_clusters, i)) or i.endswith('.txt'):
                clusters.pop(ind)
        representaive_for_clusters(path_to_clusters, clusters, ic, best_match_initial_motif, mean_threshold,
                                   z_score_threshold,
                                   top_occurrence, occurrence_threshold)

    print('Task Completed')
def representaive_for_clusters(in_folder, clusters, ic, best_match_initial_motif, mean_threshold, z_score_threshold,
                 top_occurrence, occurrence_threshold):
    for cluster in clusters:

        addition = ''
        output_folder = os.path.join(in_folder, cluster, 'repres/')
        file_name = cluster + '_rep.mlp'
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        # pwms: List[Union[Union[bytes, str], Any]] = glob(os.path.join(in_folder,cluster) + '/*.mlp')
        pwms = glob(os.path.join(in_folder, cluster) + '/*.mlp')

        if len(pwms)>=2:
            pwms_key = np.asarray([x.split('/')[-1].split('_')[0].lower() for x in pwms])


            similarityMatrix = np.zeros((len(pwms), len(pwms)))
            for index1, i in enumerate(pwms):
                matrix1, matrix_string1, maximum_feq1, total_maximum1, info1 = read_energy_matrix(i)
                matrix1 = motif_weight2p(matrix1)

                for index2, j in enumerate(pwms):
                    matrix2, matrix_string2, maximum_feq2, total_maximum2, info2 = read_energy_matrix(j)
                    matrix2 = motif_weight2p(matrix2)

                    similarityMatrix[index1, index2] = compute_similarity_score4alignment(matrix1, matrix2)
                    # print(index1, index2, similarityMatrix[index1, index2] )
            #
            # print(similarityMatrix)
            # print("Median: ",np.median(similarityMatrix, axis=0))
            # print("Sum: ",sum(similarityMatrix))

            potential_rep_index = np.unravel_index(np.argmax(sum(similarityMatrix), axis=None),sum(similarityMatrix).shape)
            potential_rep_index = int(str(potential_rep_index).split(',')[0].split('(')[1])
            matrix1, matrix_string1, maximum_feq1, total_maximum1, info1 = read_energy_matrix(pwms[potential_rep_index])

            #raw energey matrix  - Wang
            if best_match_initial_motif:
                #preselected initial matrix may have better results than the ranomdly selected median length of matrix
              potential_rep = matrix1
              #delete initinal matrix from pwm list because it will be added in the first round alignment.
              # del pwms[potential_rep_index]
              # print(info1[0])
            else:
               print('Random initial matrix')
               #use median sized motif as inital value has better results
               potential_rep = None

            uper_triangle_similarity_matrix = similarityMatrix[np.triu_indices(similarityMatrix.shape[0], k=1)]
            uper_triangle_similarity_matrix_indeces = np.triu_indices(similarityMatrix.shape[0], k=1)

            if uper_triangle_similarity_matrix.mean() > mean_threshold:  # if mean is > than 0.80, do analysis
                potential_bad_pwms = []
                if uper_triangle_similarity_matrix.std() > 0.0:
                    # print(stats.zscore(uper_triangle_similarity_matrix))
                    for index, k in enumerate(stats.zscore(uper_triangle_similarity_matrix)):
                        if k < z_score_threshold:
                            potential_bad_pwms.append(uper_triangle_similarity_matrix_indeces[0][index])
                            potential_bad_pwms.append(uper_triangle_similarity_matrix_indeces[1][index])

                    bad_pwms_indexes = get_items_upto_count(Counter(potential_bad_pwms),
                                                            int(np.ceil(len(pwms) * occurrence_threshold)),
                                                            int(np.ceil(len(pwms) * top_occurrence)))

                    del_from_pwms_list = [i[0] for i in bad_pwms_indexes]
                    bad_pwms = [pwms[x] for x in del_from_pwms_list]
                    # print("uncertain pwms removed: ", len(bad_pwms_indexes))
                    # print(bad_pwms_indexes)
                    if len(bad_pwms) > 0:
                        # total_bad_pwms_in_dbd += len(bad_pwms)
                        for x in bad_pwms:
                            pwms.remove(x)
                            # shutil.move(x, dst_for_bad_pwms)


            pwm_list = []
            #add wang
            pwm_name_list=[]
            # print("length After: ", len(pwms))
            for index, item in enumerate(pwms):
                # converted_file = convert_count_to_pwm(item, os.path.join(output_folder_for_text_report, 'count'))
                #add wang
                pwm_name_list.append(os.path.basename(item))
                pwm_list.append(read_matrix_file(item))

            #added  test without initianl selected PWM - Wang
            out_rep_1, out_count = PWMRepresentative().calculate_representative(pwm_list, potential_rep, pwm_name_list)
            ic_cutt_off = sum(out_rep_1[1,:])*ic
            out_rep = trim_edges_of_motif(out_rep_1, ic_cutt_off)
            max_index = list(map(get_char, np.argmax(out_rep, axis=1)))
            with open(str(os.path.join(output_folder, file_name)), 'w') as ff:
                ff.write(('').join(max_index) + '\t#\n')
                ff.write('0	# p-value (N/A)\n')
                ff.write('1000	# bonferroni (N/A)\n')
                ff.write('0	1	0	1	1	0	1	1	1	1	# b1, w2, b2, T-value,'
                         ' R2, alpha-w1, alpha-b1, alpha-w2, alpha-b2, beta (set to constant)\n')
                ff.write('#	a	c	g	t	#	P=>P*4\n')
                for line in range(out_rep.shape[0]):
                    ff.write('\t'.join([max_index[line]] + [str(x) for x in out_rep[line, :]]) + '\n')

if __name__ == "__main__":
    make_representative_pwm('../data/out/quality_assessed_out/', dbd='all')

#     if len(sys.argv) < 2:
#         print("Usage: python pwm_representative_test.py <input folder>")
#         exit()
#     input_folder = sys.argv[1]
#     make_representative_pwm(input_folder)
