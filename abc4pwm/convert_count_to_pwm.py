import os
import numpy as np


def read_bayesianReduce_mlp(file):
    info = list()
    motifstr = list()
    motif_energy = list()
    with open(file) as in_file:
        for line in in_file.readlines():
            if '#' in line:
                info.append(line.rstrip())
            else:
                splited = line.rstrip().split('\t')
                motifstr.append(splited[0])
                motif_energy.append(np.asarray([float(x) for x in splited[1:]]))
    return info, motifstr, np.asarray(motif_energy)


def motif_weight2p(motif_energy):
    w = np.transpose(motif_energy)
    new_w = (w + 0.25) / (w.sum(axis=0) + 1)
    psam = new_w / np.sum(new_w, axis=0)
    return psam


def convert_count_to_pwm(input_path, output_path):
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    out_file_path = os.path.join(output_path, input_path.split('/')[-1])
    with open(out_file_path, 'w') as out_file:
        info, motif_str, energy = read_bayesianReduce_mlp(input_path)
        psam = motif_weight2p(energy).transpose()
        out_file.write(('\n').join(info) + '\n')
        for line in range(psam.shape[0]):
            str_line = ('\t').join([str(x) for x in psam[line, :]])
            out_file.write(('\t').join([motif_str[line], str_line]) + '\n')
    return out_file_path
