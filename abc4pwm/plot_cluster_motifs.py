#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 17:41:09 2019

@author: omerali
"""


from fpdf import FPDF
import os, shutil
from abc4pwm.make_pwm_logo import make_logo
from PIL import Image

class plottingCluster():
    def __init__(self, path_to_clusters, output_path, clusters, dbd_plot='selected'):
        """

        :param path_to_clusters: folder path should contain clustered folders of a DBD.
        If --dbd_plot parameter is 'all' then this parameter should be path to folder which contain clustered dbds
        ;param output_path: folder where output file should go
        :param clusters: cluster numbers which you want to plot
        :param dbd_plot: dbd for which motifs to be printed. If 'all' then 'path_to_clusters' should be path
        to folder which contain clustered dbds
        """

        path_to_clusters = os.path.join(path_to_clusters,'')
        output_path = os.path.join(output_path,'')

        if 'all' in dbd_plot or '*' in dbd_plot:

            dbds = [i for i in sorted(os.listdir(path_to_clusters)) if not i.endswith('.DS_Store')]

            for dbd in dbds:
                if not os.path.exists(os.path.join(path_to_clusters,dbd,'out/')):
                    continue
                clusters = [i for i in sorted(os.listdir(os.path.join(path_to_clusters,dbd,'out/'))) if not i.endswith('.DS_Store')]
                print('Preparing motifs...')

                for c in clusters:
                    print('Cluster Number: ' + c)
                    self.plotclustersummary(self, os.path.join(path_to_clusters,dbd,'out/'), output_path, c)

            print("Task Completed. \nPlease see ", output_path, "to see plotted clusters in pdf files")

            exit()

        if 'all' in clusters:
            clusters = [i for i in sorted(os.listdir(path_to_clusters)) if not i.endswith('.DS_Store')]

        else:
            clusters = clusters.split(',')
        if not os.path.exists(output_path):
            os.makedirs(output_path, exist_ok=True)


        print('Preparing motifs...')

        for c in clusters:
            print('Cluster Number: ' + c)
            self.plotclustersummary(self, path_to_clusters, output_path, c)

        print("Task Completed. \nPlease see ",output_path,"to see plotted clusters in pdf files")

    def croping(self, img):
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
        # Shows the image in image viewer

        # Resizing
        basewidth = 200
        wpercent = (basewidth/float(im1.size[0]))
        hsize = int((float(im1.size[1])*float(wpercent)))
        im1 = im1.resize((basewidth,hsize), Image.ANTIALIAS)
        im1.save(temp)


    @staticmethod
    def plotdirectory(self, directoryname):                   #plotting motif images

        clustercheckspwms = [i for i in sorted(os.listdir(directoryname)) if i.endswith('.mlp') ]
        for i in clustercheckspwms:
            print(i)
            if not i.endswith('.png'):
                make_logo(os.path.join(directoryname, i))



    @staticmethod
    def folderizeclusters(folder_path):

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
    def plotclustersummary(self, indir, output_path, cluster):


        pdf = FPDF()
        pdf.add_page()

        path_to_clusters = indir
        path_to_mlps_of_cluster = os.path.join(path_to_clusters,str(cluster))
        self.plotdirectory(self, path_to_mlps_of_cluster)



        src = path_to_mlps_of_cluster
        dst = src + '/pngs'
        if not os.path.exists(dst):
            os.makedirs(dst)

        files = [i for i in os.listdir(src) if i.endswith(".png") and os.path.isfile(os.path.join(src, i))]
        for f in files:
            shutil.move(os.path.join(src, f),os.path.join(dst, f))



        path_to_pngs_of_cluster = os.path.join(path_to_mlps_of_cluster, 'pngs')

        if os.path.exists(os.path.join(src, 'repres')):
            path_to_rep_of_cluster = os.path.join(path_to_mlps_of_cluster, 'repres')
            self.plotdirectory(self, path_to_rep_of_cluster)
        else:
            print("Program Terminated. \n "
                  "Representative motif for this cluster does not exist. Please make representative first by "
                  "using make_representative function of the package.")
            exit()

        motif_images_clusterf = sorted(os.listdir(path_to_pngs_of_cluster))
        path_to_txts = output_path + '/txts/'
        if not os.path.exists(path_to_txts):
            os.makedirs(path_to_txts, exist_ok=True)
        with open(os.path.join(path_to_txts, str(indir).split('/')[-3] + '_'+ str(cluster)+'.txt'), 'w') as f:
            if len(files)>=2:
                rep_motif = path_to_rep_of_cluster
                image_path_rep = os.path.join(path_to_rep_of_cluster, str(cluster)+'_rep.mlp.png')
                pdf.set_font('Arial', '', 14)
                pdf.cell(40, 5, 'DBD \t\t\t Cluster: ', 0, 1)
                pdf.set_font('Arial', 'B', 14)
                pdf.cell(40, 5, str(indir).split('/')[-3] + ' \t\t\t '+ cluster, 0, 1)

                pdf.set_font('Arial', '', 8)
                pdf.image(image_path_rep)

                # Move to 10 cm to the right
                pdf.cell(100)
                # Centered text in a framed 20*10 mm cell and line break
                v = rep_motif
                caption_of_motif = 'Representative / Consencus Motif '
                pdf.cell((len(caption_of_motif) + 24), 6, caption_of_motif, 1, 1, 'C')
                print(rep_motif, file=f)



            for index, motif in enumerate(sorted(motif_images_clusterf)):
                image_path = os.path.join(path_to_pngs_of_cluster,motif)
                self.croping(image_path)
                pdf.image(image_path)
                pdf.set_font('Arial', '', 8)
                # Move to 10 cm to the right
                pdf.cell(100)
                # Centered text in a framed 20*10 mm cell and line break
                v = motif
                tmp1 = v.split('-')[-1]
                tmp2 = v.split('_')[0:2]

                stri = ''
                stri = tmp2[0]+tmp2[1]
                stri = stri+'_'+str(tmp1)
                pdf.cell((len(stri)+24), 6, str(index) + ': '+ stri, 1, 1, 'C')
                print(index, motif, file=f)
        path_to_pdfs = output_path
        if not os.path.exists(path_to_pdfs):
            os.makedirs(path_to_pdfs,exist_ok=True)
        pdfname = os.path.join(path_to_pdfs, str(indir).split('/')[-3] + '_' + str(cluster) + '.pdf')
        pdf.output(pdfname, "F")


if __name__ == "__main__":

    clusterPlotObj = plottingCluster('../data/out/quality_assessed_out/bZIP/out/','../data/out/plots/pdfs/', 'all')


