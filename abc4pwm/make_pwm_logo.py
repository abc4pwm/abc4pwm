import sys
import subprocess
import os
import numpy as np
from abc4pwm.energy_to_p import read_energy_matrix



def make_logo(pwm_file):
    pwm, matrix_string, maximum_feq, total_maximum, info = read_energy_matrix(pwm_file)
    fields = os.path.basename(pwm_file).split("_")
    if len(fields) > 1:
        tf = fields[1]
        num = fields[-1].split(".")[0]
        comment = "TF " + tf + "/" + num
    else:
        comment = ""

    write_logo_file(pwm, pwm_file + ".png", comment)


def write_logo_file(pwm, file_name, comment):

    FNULL = open(os.devnull, 'w')
    p = subprocess.Popen(["weblogo", "-o", file_name, "-F", "png", "-D", "transfac", #"-U", "kT",
                          "--scale-width", "no", "-c", "classic", "-P", "", #"entropy: " + str(round(entropy, 1)),
                          "-t", comment, "--composition", "equiprobable", "--errorbars", "no", "--show-yaxis", "no",
                          "--resolution", "100"],
                         stdin=subprocess.PIPE, stderr=FNULL)

    transfac_matrix = 'PO\tA\tC\tG\tT\n'

    for i in range(len(pwm)):
        transfac_matrix += str(i + 1) + "\t" + "\t".join(map(str, pwm[i])) + "\n"

    p.communicate(transfac_matrix.encode('utf-8'))

def read_mlp(pwm_file, return_biases=False, return_r2=False):
    with open(pwm_file) as f:
        lines = [l.replace("\n", "") for l in f]

    b1, w2, b2, T, r2 = list(map(float, lines[3].split("\t")[:5]))
    data = [l.split("\t")[1:5] for l in lines[5:]]
    data = np.array(data, float)
    if w2 < 0:
        data = 1.0 / data

    if return_biases:
        if return_r2:
            return data, b1, w2, b2, r2
        else:
            return data, b1, w2, b2
    else:
        if return_r2:
            return data, r2
        else:
            return data

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python make_pwm_logo.py <mlp files>")
        exit()

    pwm_files = sys.argv[1:]
    for f in pwm_files:
        print(f)
        make_logo(f)

