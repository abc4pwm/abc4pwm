import glob
import os
import numpy as np
from PIL import Image
from abc4pwm.energy_to_p import read_energy_matrix, read_energy_matrix_jaspar, read_energy_matrix_transfac
from abc4pwm.similarity_score import compute_similarity_score4alignment
import shutil
from fpdf import FPDF
from abc4pwm.make_pwm_logo import make_logo
from abc4pwm.convert_count_to_pwm import motif_weight2p, read_bayesianReduce_mlp
from abc4pwm.format_conversion import Conversion
from itertools import islice

def take(n, iterable):
    "Return first n items of the iterable as a list"
    return list(islice(iterable, n))

def searching(pwm, db_path, out, tf_name = '', db_type = 'path'  , db_format = 'abc4pwm' ,n = 5,
              input_count = False, db_count = False, db_file_type = '.mlp', input_file_type = '.mlp', input_prop = False, db_prob = False):

    """

    :param pwm: input file for searching
    :param db_path: database folder path of target files with which you want to search and match
    :param out: output folder where resulted pdf file wil store
    :param db_type: type of database: Options are folder, path(in case of abc4pwm setting),
    :param db_format: database file formate, options are abc4pwm, transfac, jaspar
    :param n: number of top matches
    :param input_count: True if input file contains counts
    :param db_count: True if database file contains counts
    :param db_file_type: mention the extension of file type. e.g., .mlp, .txt
    :param input_file_type: mention the extension of file type. e.g., .mlp, .txt
    :param input_prop: True if input matrix is already in probability
    :param db_prob: True if database matrix is already in probability
    :return:
    """


    print("Searching Task started. . . .  \n")


    if not os.path.exists(out):
        os.makedirs(out, exist_ok=True)
    empty_dir(out)

    matrix1, matrix_string1, maximum_feq1, total_maximum1, info1 = read_energy_matrix(pwm)
    if not input_count:
        if not input_prop:
            matrix1 = motif_weight2p(matrix1)
    else:
        normalized_matrix = (matrix1.transpose() / maximum_feq1).transpose()
        normalized_matrix = (normalized_matrix.transpose() + 0.25 / (np.sum(normalized_matrix, axis=1) + 1)).transpose()
        normalized_matrix = (normalized_matrix.transpose() / np.sum(normalized_matrix, axis=1)).transpose()
        matrix1 = normalized_matrix



    score_list= {}

    if 'path' in db_type:
        dbds = [i for i in os.listdir(db_path) if not i.endswith('.DS_Store')]

        for dbd in dbds:
            if os.path.exists(os.path.join(db_path, dbd, 'out/')):
                clusters = [i for i in sorted(os.listdir(os.path.join(db_path, dbd, 'out/'))) if
                            not i.endswith('.DS_Store')]

                for cluster in clusters:
                    rep = os.path.join(db_path, dbd, 'out/', cluster, 'repres/_'+ str(cluster) +'_rep.mlp')
                    if os.path.exists(rep):
                        matrix2, matrix_string2, maximum_feq2, total_maximum2, info2 = read_energy_matrix(rep)
                        if not db_count:
                            if not db_prob:
                                matrix2 = motif_weight2p(matrix2)
                        else:
                            normalized_matrix = (matrix2.transpose() / maximum_feq2).transpose()
                            normalized_matrix = (normalized_matrix.transpose() + 0.25 / (
                                        np.sum(normalized_matrix, axis=1) + 1)).transpose()
                            normalized_matrix = (
                                        normalized_matrix.transpose() / np.sum(normalized_matrix, axis=1)).transpose()
                            matrix2 = normalized_matrix
                        score = compute_similarity_score4alignment(matrix1, matrix2)
                        score_list[rep] = score

    elif db_type in 'folder':

        pwms = glob.glob(db_path+ "*"+tf_name+"*"+db_file_type)
        if 'abc4pwm' in db_format:
            for motif in pwms:
                if os.path.exists(motif):
                    matrix2, matrix_string2, maximum_feq2, total_maximum2, info2 = read_energy_matrix(motif)
                    if not db_count:
                        if not db_prob:
                            matrix2 = motif_weight2p(matrix2)
                    else:
                        normalized_matrix = (matrix2.transpose() / maximum_feq2).transpose()
                        normalized_matrix = (normalized_matrix.transpose() + 0.25 / (
                                np.sum(normalized_matrix, axis=1) + 1)).transpose()
                        normalized_matrix = (
                                normalized_matrix.transpose() / np.sum(normalized_matrix, axis=1)).transpose()
                        matrix2 = normalized_matrix
                    score = compute_similarity_score4alignment(matrix1, matrix2)
                    score_list[motif] = score
        elif 'transfac' in db_format:
            for motif in pwms:
                if os.path.exists(motif):
                    matrix2, matrix_string2, maximum_feq2, total_maximum2, info2 = read_energy_matrix_transfac(motif)
                    matrix2 = motif_weight2p(matrix2)

                    score = compute_similarity_score4alignment(matrix1, matrix2)
                    score_list[motif] = score
        elif 'jaspar' in db_format:
            for motif in pwms:
                if os.path.exists(motif):
                    matrix2, matrix_string2, maximum_feq2, total_maximum2, info2 = read_energy_matrix_jaspar(motif)
                    matrix2 = motif_weight2p(matrix2)
                    score = compute_similarity_score4alignment(matrix1, matrix2)
                    score_list[motif] = score
        else:
            print("Your mentioned " + db_format + " format conversion is not available. Please write to us"
                                                  "if you want a conversion format")
            exit()


    else:
        for motif in db_path:
            if os.path.exists(motif):
                if 'abc4pwm' in db_format:
                    matrix2, matrix_string2, maximum_feq2, total_maximum2, info2 = read_energy_matrix(motif)
                elif 'transfac' in db_format:
                    matrix2, matrix_string2, maximum_feq2, total_maximum2, info2 = read_energy_matrix_transfac(motif)
                elif 'jaspar' in db_format:
                    matrix2, matrix_string2, maximum_feq2, total_maximum2, info2 = read_energy_matrix_jaspar(motif)
                else:
                    print("Your mentioned " + db_format + " format conversion is not available. Please write to us"
                                                          "if you want a conversion format")
                    exit()
                matrix2 = motif_weight2p(matrix2)
                score = compute_similarity_score4alignment(matrix1, matrix2)
                score_list[motif] = score


    rank = {k: v for k, v in sorted(score_list.items(), key=lambda item: item[1], reverse=True)}
    n_items = take(n, rank.items())
    matched_items = []
    score_list_top = []
    for item in n_items:
        matched_items.append(str(item[0]).split('/')[-1])
        score_list_top.append(item[1])
        shutil.copy(item[0], out)



    if 'abc4pwm' in db_format:
        plot_search_result(pwm, out, out, n_items, db_file_type, input_file_type)
    elif 'transfac' in db_format:
        Conversion(out, 'transfac2abc4pwm' , os.path.join(out,'converted'))
        plot_search_result(pwm, os.path.join(out, 'converted'), out, n_items, db_file_type, input_file_type)
    elif 'jaspar' in db_format:
        Conversion(out, 'jaspar2abc4pwm', os.path.join(out,'converted'))
        plot_search_result(pwm, os.path.join(out, 'converted'), out, n_items, db_file_type, input_file_type)
    else:
        print("Specified format ", db_format, " is not supported. ")
        exit()





    print("Searching Task finished. Please see results in ", out)
    return matched_items, score_list_top

def plot_search_result(pwm, resulted_matches_path, out, n_items, db_file_type, input_file_type):
    pdf = FPDF()
    pdf.add_page()
    plotdirectory(resulted_matches_path, db_file_type)

    src = resulted_matches_path
    dst = src + '/pngs'
    if not os.path.exists(dst):
        os.makedirs(dst)

    files = [i for i in os.listdir(src) if i.endswith(".png") and os.path.isfile(os.path.join(src, i))]
    for f in files:
        shutil.move(os.path.join(src, f), os.path.join(dst, f))

    path_to_pngs_of_cluster = os.path.join(resulted_matches_path, 'pngs')

    if not os.path.exists(os.path.join(resulted_matches_path,"in_pwm")):
        os.makedirs(os.path.join(resulted_matches_path,"in_pwm"), exist_ok=True)
    empty_dir(os.path.join(resulted_matches_path,"in_pwm"))
    shutil.copy(pwm, os.path.join(resulted_matches_path,"in_pwm"))

    path_to_rep_of_cluster = os.path.join(resulted_matches_path, 'in_pwm')
    plotdirectory(path_to_rep_of_cluster, input_file_type)

    scores_for_printing = []
    motif_images_clusterf = []
    for item in n_items:
        motif_images_clusterf.append(str(item[0]+'.png').split('/')[-1])
        scores_for_printing.append(item[1])
    with open(os.path.join(out,'search_results.txt'), 'w') as f:
        f.writelines("Matched File name\tSimilairty Score\n")
        for i,file in enumerate(motif_images_clusterf):
            f.writelines(str(file).split('.png')[0]+'\t'+str(scores_for_printing[i])+'\n')

    # motif_images_clusterf = sorted(os.listdir(path_to_pngs_of_cluster))
    # print(motif_images_clusterf)
    rep_motif = path_to_rep_of_cluster
    image_path_rep = os.path.join(path_to_rep_of_cluster, str(pwm).split('/')[-1] +'.png')
    # self.croping(image_path_rep)
    pdf.set_font('Arial', '', 8)
    # croping(image_path_rep, 'rep')
    pdf.image(image_path_rep)

    # Move to 10 cm to the right
    pdf.cell(100)
    # Centered text in a framed 20*10 mm cell and line break
    v = rep_motif
    caption_of_motif = 'Input Motif '
    pdf.cell((len(caption_of_motif) + 24), 6, caption_of_motif, 1, 1, 'C')

    # print(rep_motif, file=f)

    caption_of_motif = 'Matching Results: '
    pdf.cell((len(caption_of_motif) + 24), 6, caption_of_motif, 1, 1, 'C')

    for index, motif in enumerate(motif_images_clusterf):
        image_path = os.path.join(path_to_pngs_of_cluster, motif)
        croping(image_path)
        pdf.image(image_path)
        pdf.set_font('Arial', '', 8)
        # Move to 10 cm to the right
        pdf.cell(100)
        # Centered text in a framed 20*10 mm cell and line break
        v = motif
        tmp1 = v.split('-')[-1]
        # tmp2 = v.split('_')[0:2]

        stri = ''
        stri = "{:0.3f}".format(scores_for_printing[index])
        stri = str(stri) + '->' + str(tmp1)
        pdf.cell((len(stri) + 24), 6, str(index) + ': ' + stri, 1, 1, 'C')
        # print(index, motif, file=f)
    path_to_pdfs = out
    if not os.path.exists(path_to_pdfs):
        os.makedirs(path_to_pdfs, exist_ok=True)
    pdfname = os.path.join(out + 'search_result.pdf')
    pdf.output(pdfname, "F")

def croping(img, type = 'common'):           #cropping image
    # Opens a image in RGB mode
    im = Image.open(img)
    temp = img
    # Size of the image in pixels (size of orginal image)
    # (This is not mandatory)
    width, height = im.size

    # Setting the points for cropped image
    left = 0
    top = height / 2
    right = width
    bottom = height

    # Cropped image of above dimension
    # (It will not change orginal image)
    im1 = im.crop((left, top, right, bottom))
#    im1.save(temp)
    # Shows the image in image viewer
#    im1.show()

    if 'common' in type:
    # Resizing
        basewidth = 200
        wpercent = (basewidth/float(im1.size[0]))
        hsize = int((float(im1.size[1])*float(wpercent)))
        im1 = im1.resize((basewidth,hsize), Image.ANTIALIAS)
    im1.save(temp)
#    im1.show()



def plotdirectory(directoryname, db_file_type):                   #plotting motif images

    clustercheckspwms = [i for i in sorted(os.listdir(directoryname)) if i.endswith(db_file_type) ]

    for i in clustercheckspwms:
        print(i)
        if not i.endswith('.png'):
            make_logo(os.path.join(directoryname, i))

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



if __name__ == '__main__':
    pwms = glob.glob("../data/in/in_pwms/*.mlp")
    out = "../data/out/search_out/"

    # db_path = "../data/out/quality_assessed/"
    # db_path_list = []
    # dbds = [i for i in os.listdir(db_path) if not i.endswith('.DS_Store')]
    # for dbd in dbds:
    #     if os.path.exists(os.path.join(db_path, dbd, 'out/')):
    #         clusters = [i for i in sorted(os.listdir(os.path.join(db_path, dbd, 'out/'))) if
    #                     not i.endswith('.DS_Store')]
    #
    #         for cluster in clusters:
    #             rep = os.path.join(db_path, dbd, 'out/', cluster, 'repres/_' + str(cluster) + '_rep.mlp')
    #             db_path_list.append(rep)
    # print(db_path_list)
    # searching(pwms[4], db_path_list, 'list')

    # #
    # db_path = "../../../../../Downloads/JASPAR2020_CORE_vertebrates_non-redundant_pfms_transfac/"
    # searching(pwms[6], db_path, out, 'folder', 'transfac', 5)

    #
    # db_path = "../data/out/convert/to_jaspar/"
    # searching(pwms[6], db_path, out, 'folder', 'jaspar', 5)

    db_path = "../data/in/in_pwms/"
    searching(pwms[16], db_path, out, db_type='folder', db_format= 'abc4pwm', n=5)

    #SGD
    # db_path = "../data/sdg/SGD_yeast_p/"
    # searching('../data/sdg/demo_out/swi4_2k/demo_swi4_2Kseq_cacgaaaarepeatMasked.txt_mlpout_L11_1.mlp', db_path, out, 'folder',
    #           'abc4pwm', 5, db_file_type='.txt', db_count=False, input_prop=False, input_count=True, db_prob=True  )


    # db_path = "../data/out/quality_assessed_out/"
    # searching(pwms[4], db_path,out)

