import numpy as np
import os

def convert_energy_matrix_to_pwm(path, output_path):
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    matrix, matrix_string, maximum_feq, total_maximum, info = read_energy_matrix(path)
    normalized_matrix = (matrix.transpose() / maximum_feq).transpose()
    normalized_matrix = (normalized_matrix.transpose() + 0.25 / (np.sum(normalized_matrix, axis=1) + 1)).transpose()
    normalized_matrix = (normalized_matrix.transpose() / np.sum(normalized_matrix, axis=1)).transpose()
    out_file_path = os.path.join(output_path, path.split('/')[-1].split('.mlp')[0] + '_p.mlp')
    with open(out_file_path, 'w') as out_file:
        out_file.writelines(info)
        for i, line in enumerate(matrix_string):
            str_line = ('\t').join(
                [matrix_string[i], ('\t').join([str(x) for x in normalized_matrix[i, :]])]) + '\n'
            out_file.write((str_line))
    return str(out_file_path)


def read_energy_matrix(path):
    with open(path) as in_file:
        lines = in_file.readlines()
    start_index = 0
    for id, line in enumerate(lines):
        if '#' in line:
            start_index = id
    start_index += 1
    info = lines[0:start_index]
    matrix_string = []
    matrix_value = lines[start_index:]
    matrix = []
    for line in matrix_value:
        splited = line.rstrip().split('\t')
        matrix_string.append(splited[0])
        ls = [float(x) for x in splited[1:]]
        matrix.append(ls)
    matrix = np.asarray(matrix)
    matrix = np.clip(matrix, np.exp(-10), np.exp(10))

    return matrix, matrix_string, np.sum(matrix, axis=1), np.max(matrix), info

def read_energy_matrix_transfac(path):
    with open(path) as in_file:
        lines = in_file.readlines()
    start_index = 0
    ei = 0
    for id, line in enumerate(lines):
        if 'P0' in line:
            start_index = id

    start_index += 1
    info = lines[0:start_index]
    matrix_string = []
    for i, line in enumerate(lines[start_index:]):
        if 'XX' in line:
            ei = i
    end_index = start_index + (ei)
    matrix_value = lines[start_index:end_index]

    matrix = []
    for line in matrix_value:
        splited = line.rstrip().split('\t')
        matrix_string.append(splited[-1])
        ls = [float(x) for x in splited[1:-1]]
        matrix.append(ls)
    matrix = np.asarray(matrix)
    matrix = np.clip(matrix, np.exp(-10), np.exp(10))

    return matrix, matrix_string, np.sum(matrix, axis=1), np.max(matrix), info

def read_energy_matrix_jaspar(path):
    with open(path) as in_file:
        lines = in_file.readlines()
    start_index = 0
    ei = 0
    info = lines[0]
    matrix_string = []
    matrix_value = lines[1:]
    matrix = []
    for line in matrix_value:
        ls = line.split('[')[1].split(']')[0].split('\t')
        lss = [float(x) for x in ls]
        matrix.append(lss)
    matrix = np.asarray(matrix)
    matrix = np.clip(matrix, np.exp(-10), np.exp(10))

    return matrix.transpose(), matrix_string, np.sum(matrix, axis=1), np.max(matrix), info

if __name__ == "__main__":
    # convert_energy_matrix_to_pwm('test_energy2p/demo_swi4_0.05Kseq_cacgaaaa.txt_mlpout_L10_1.mlp','test_energy2p/out/')
    # read_energy_matrix('../data/in/in_pwms/AIRE_2_from_AIRE_transfac_M01000.mlp')
    # read_energy_matrix_transfac('../../../../../Downloads/JASPAR2020_CORE_vertebrates_non-redundant_pfms_transfac/MA0006.1.transfac')
    read_energy_matrix_jaspar('../data/out/convert/to_jaspar/GATA_known12_from_GATA2_2_from_GATA2_jaspar_MA0036.1.mlp')