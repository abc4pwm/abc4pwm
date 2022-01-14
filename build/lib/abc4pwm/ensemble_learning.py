import numpy as np
import random
import multiprocessing as mp
import os
import shutil
import sys
from functools import partial
import time
from multiprocessing import Pool
from joblib import Parallel, delayed
import subprocess


class EnsembleLearning():
    def __init__(self, opt_dependence = 0,
                 numP = 15,
                 opt_numOfWeakReads = 0,
                 select_random_n= 200,
                 opt_strong_expFile = 'demo_in/demo_swi4_0.5Kseq_cacgaaaa_500.txt',
                 opt_weak_expFile = '',
                 opt_seqFile='demo_in/demo_swi4_0.5Kseq_cacgaaaa_500.fa',
                 opt_out = 'demo_out_ensemble/swi4/',
                 opt_loops=3,
                 opt_min_L = 9,
                 opt_max_L = 9,
                 opt_iteration=500,
                 opt_p_value=0.0001,
                 opt_strand=0,
                 opt_normalization = 2,
                 max_processors=0):

        if not os.path.exists(opt_out):
            os.makedirs(opt_out, exist_ok=True)
        if not os.path.exists(opt_strong_expFile):
            print(opt_strong_expFile, " does not exsist.")
            exit()
        if not os.path.exists(opt_seqFile):
            print(opt_seqFile, " does not exsist.")
            exit()

        self.empty_dir(opt_out)

        strong_exp, strong_genes, strong_head = self.read_exp_file(opt_strong_expFile)

        new_weak_exp = {}
        new_weak_genes = []
        if opt_numOfWeakReads and opt_numOfWeakReads > 0:
            weak_exp, weak_genes, weak_head = self.read_exp_file(opt_weak_expFile)

            print("weak exp: "+ weak_exp + "\n weak genes:  "+weak_genes+  "\n: weak head: "+ weak_head)

            select_weakreads  = opt_numOfWeakReads;
            select_strongreads = select_random_n;
            print("Randomly select "+ select_weakreads + " reads from weak input files\n")
            print("Randomly select "+ select_strongreads + "reads from strong input files\n")

             # exclude strong reads from weak reads

            new_weak_exp = set(strong_exp).intersection(set(weak_exp))
            new_weak_genes = set(strong_genes).intersection(set(strong_genes))

            print("There are "+ len(new_weak_exp) +" reads overlapping between two files \n")
            print("Within "+ len(weak_genes) + " reads, there are " + len(new_weak_genes) +" reads in the new weak exp files \n")



        pids = []
        pid = 0
        job_loops = 0

        if max_processors != 0:
            processors= max_processors
        else:
            processors= mp.cpu_count()

        #added wang parallel in number of random selections
        all_processes=[]
        for i in range(1,numP+1):
            all_processes.extend(self.random_run_bayesPI(i, select_random_n, opt_numOfWeakReads, opt_seqFile, opt_strong_expFile,
                                    opt_iteration, opt_loops, opt_normalization,
                                    opt_min_L, opt_max_L, opt_p_value, opt_out,
                                    new_weak_genes, strong_genes, new_weak_exp,
                                    strong_exp, strong_head, opt_strand, opt_dependence, max_processors=max_processors))
        #added wang
        #for p in all_processes:
        #    p.communicate()
        num_of_processes= max_processors
        run_processes=[]
        while len(all_processes) >0:
            if len(all_processes) < num_of_processes:
               num_of_processes= len(all_processes)
            for i in range(0, num_of_processes):
                run_processes.append(subprocess.Popen(all_processes.pop(), shell=True)) 
            for p in run_processes:
                p.communicate()
            run_processes=[]            

        print("Done with splitting jobs \n")

        #end added

        # remove temporary files
        for job in range(1,numP+1):
            temp_name = "_".join([str(opt_strong_expFile), str(select_random_n), "plus_weak", str(opt_numOfWeakReads), str(job)])
            del temp_name

    def random_run_bayesPI(self, job_loops, opt_numOfStrongReads, opt_numOfWeakReads, opt_seqFile,opt_strong_expFile,
                                opt_iteration, opt_loops, opt_normalization,
                                opt_min_L, opt_max_L, opt_p_value, opt_out,
                                new_weak_genes, strong_genes,  new_weak_exp,
                                strong_exp, strong_head, opt_strand, opt_dependence, max_processors):


        # random selection in weak and strong peak calls and export data
        weak_keys = []
        new_weak_keys = []
        strong_keys = list(strong_exp.keys())
        random.shuffle(strong_keys)
        new_strong_keys = strong_keys[:opt_numOfStrongReads]
        if opt_numOfWeakReads and opt_numOfWeakReads > 0:
            weak_keys = self.randomize(list(new_weak_exp.keys()))
            new_weak_keys = weak_keys[:-1]
            all_new_keys = new_weak_keys + new_strong_keys

        all_new_keys = new_strong_keys

        outExp = "_".join([str(opt_strong_expFile), str(opt_numOfStrongReads), "plus_weak", str(opt_numOfWeakReads), str(job_loops)])


        self.export_tag_counts(outExp, strong_head, all_new_keys, strong_exp, new_weak_exp)


        opt_numP = opt_max_L - opt_min_L + 1
        #added wang
        return self.run_bayesPI_parallel_by_motifLength(opt_dependence, opt_numP, outExp, opt_seqFile, opt_loops, opt_iteration, opt_normalization,opt_out, opt_min_L, opt_max_L, opt_p_value, opt_strand, max_processors=max_processors)


    def run_bayesPI_parallel_by_motifLength(self, opt_dependence, opt_numP, opt_expFile , opt_seqFile , opt_loops , opt_iteration , opt_normalization , opt_out , opt_min_L , opt_max_L , opt_strand, opt_p_value = 0.0001 , flank_sequence = 0, max_evidence ='', max_processors=0):


        num_of_motif_length =opt_max_L -opt_min_L + 1;
        num_of_motifs_in_each_process = int(num_of_motif_length / opt_numP)
        num_of_motifs_in_process ={}
        print("min_L= "+ str(opt_min_L)+ ", max_L= "+str(opt_max_L)+", numP="+str(opt_numP)+", in "+str(num_of_motifs_in_each_process))
        
        #added wang
        for i in range(1,opt_numP+1):
            temp_min = opt_min_L+num_of_motifs_in_each_process * (i-1)
            temp_max = temp_min +num_of_motifs_in_each_process - 1
            if i == opt_numP and temp_max > opt_max_L:
                temp_max = opt_max_L
            if i == opt_numP and temp_max < opt_max_L:
                temp_max =opt_max_L

            num_of_motifs_in_process[i] = [temp_min, temp_max]
            print(num_of_motifs_in_process[i] , temp_min, temp_max , " \n")

        pids = []
        pid = 0
        opt_flank_sequence = 0
        opt_max_evidence = 3


        if max_processors != 0:
            processors= max_processors
        else:
            processors= mp.cpu_count()

        arg_list = []

        #working======
        #added wang
        processes =[]
        for job_loops in range(opt_numP):

            temp_min_L = num_of_motifs_in_process[job_loops + 1][0]
            temp_max_L = num_of_motifs_in_process[job_loops+1][1]
            #pool.apply_async(self.run_bayesPI_in_motif_parallel, args= (job_loops, opt_seqFile, opt_expFile, opt_iteration, opt_loops, opt_normalization, opt_dependence, temp_min_L,  temp_max_L,  opt_p_value,
            #                              opt_out, opt_strand,  opt_max_evidence, opt_flank_sequence))
            tmp_cmd=self.run_bayesPI_in_motif_parallel(job_loops, opt_seqFile, opt_expFile, opt_iteration, opt_loops, opt_normalization,
                                 opt_dependence, temp_min_L,  temp_max_L,  opt_p_value, opt_out, opt_strand,  opt_max_evidence, opt_flank_sequence)
            print(tmp_cmd)
            #added wang
            #processes.append(subprocess.Popen(tmp_cmd,shell=True ))  
            processes.append(tmp_cmd)
        #============
    
        print("Done with splitting motifs")
        return processes


    def run_bayesPI_in_motif_parallel(self, job_loops, opt_seqFile, opt_expFile, opt_iteration, opt_loops, opt_normalization, opt_dependence, opt_min_L,  opt_max_L, opt_p_value,
                                          opt_out, opt_strand,  opt_max_evidence, opt_flank_sequence ):
        order = len(opt_expFile)
        ln = int(order*0.6)
        random = opt_expFile[order-ln+1:(order-ln)+ln]

        # if 'darwin' in sys.platform:
        #     bayes_binary = './bin/Mac/bayesPI'
        # elif 'linux' in sys.platform:
        #     bayes_binary= './bin/Linux/bayesPI'
        # else:
        #     print("Unknown platform. Can be run on Linux and Mac only")
        #     exit(1)

        bin_folder = os.path.abspath(
            os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "abc4pwm", "bin"))

        platform_subfolders = {"darwin": "Mac", "linux": "Linux"}

        platform = sys.platform
        bin_subfolder = None
        for platform_name, subfolder in platform_subfolders.items():
            if platform.startswith(platform_name):
                bin_subfolder = os.path.join(bin_folder, subfolder)
                break

        if bin_subfolder is None:
            print(
                "Unknown platform:" + platform + ". This package has binary components for the following platforms:" +
                ", ".join(platform_subfolders.keys()))
            exit(1)

        bayes_binary = os.path.join(bin_subfolder, "bayesPI")

        p1 = '-flank_sequence='+str(opt_flank_sequence)
        p2='-max_evidence='+str(opt_max_evidence)
        p3='-dependence='+str(opt_dependence)
        p4='-strand='+str(opt_strand)
        p5='-max_loop='+str(opt_loops)
        p6='-max_iteration='+str(opt_iteration)
        p7='-normalize='+str(opt_normalization)
        p8='-exp='+str(opt_expFile)
        p9='-seq='+str(opt_seqFile)
        p10='-min_L='+str(opt_min_L)
        p11='-max_L='+str(opt_max_L)
        p12='-p_value='+str(opt_p_value)
        p13='-out='+str(opt_out)
        p14='2>'+str(opt_out)+"log_"+str(opt_min_L)+"_"+str(job_loops)+"_"+str(random)

        cmd = bayes_binary + ' ' +p1+ ' ' + p2+ ' ' + p3+ ' ' + p4+ ' ' + p5+ ' ' + p6+ ' ' + p7+ ' ' + p8+ ' ' + p9+ ' ' + p10+ ' ' + p11+ ' ' + p12+ ' ' + p13 + ' ' + p14
        return cmd
        #os.system(cmd)

    def export_tag_counts(self, outfile, column_names_weak, ids, strong, weak):
        with open(outfile, 'w') as f:
            new_head = '\t'.join(["ID",str(column_names_weak)])

            f.writelines(new_head)
            f.writelines("\n")
            for id in ids:
                if id in strong:
                    f.writelines('\t'.join([id, strong[id]]))
                    f.writelines("\n")
                elif id in weak:
                    f.writelines('\t'.join([str(id), str(weak[id])]))
                    f.writelines("\n")
                else:
                    print("Not found")


    def read_exp_file(self, file):
        print("\n\t Getting gene expression data ", file, "\n")
        genes = []
        expressions = []
        hashes = {}
        head = ''
        with open(file) as f:
            lines = f.readlines()
            head = str(lines[0]).split('\t')[1].strip()
            for line in lines[1:]:
                gene = str(line).split('\t')[0].strip()
                ex = str(line).split('\t')[1].split('\n')[0].strip()
                genes.append(gene)
                expressions.append(ex)
                hashes[str(gene)]= ex

        unique_genes = np.unique(genes)
        # print(expressions)
        # print(hashes)
        # print(head)
        print("\t*", len(genes) ," total genes in 1 experiments (missing 0 of ",len(expressions)," expression values; ", len(genes) - len(unique_genes) ," duplicate entries")
        return hashes, unique_genes, head
    def empty_dir(self, folder):
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

if __name__ == '__main__':
    expfile = '../data/in/expression/demo_swi4_0.5Kseq_cacgaaaa_500.txt'
    seqfile = '../data/in/seq/demo_swi4_0.5Kseq_cacgaaaa_500.fa'
    random_n = 200
    out = '../data/out/ensemble/swi4/'
    EnsembleLearning(opt_strong_expFile=expfile, opt_seqFile=seqfile, select_random_n=random_n, opt_out=out, max_processors=5)
