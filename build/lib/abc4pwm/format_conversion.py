import os
from abc4pwm.energy_to_p import read_energy_matrix, \
                                read_energy_matrix_transfac, \
                                read_energy_matrix_jaspar
import glob
import shutil



class Conversion():
    def __init__(self, input_path, in2out, out):
        """

        :param input_path: path to the input files that you want to convert
        :param in2out: speicify format from and to
        :param out: out path where you want to store the output files.
        """

        if not os.path.exists(out):
            os.makedirs(out, exist_ok=True)

        self.empty_dir(out)

        pwms = glob.glob(input_path+"*.mlp")

        if in2out == 'abc4pwm2transfac':
            print("\nConverting from abc4pwm to TRANSFAC.\n")
            for pwm in pwms:
                self.abc4pwm2transfac(pwm, out)

        elif in2out == 'transfac2abc4pwm':
            print("\nConverting from TRANSFAC to abc4pwm.\n")
            for pwm in pwms:
                self.transfac2abc4pwm(pwm, out)

        elif in2out == 'abc4pwm2jaspar':
            print("\nConverting from abc4pwm to JASPAR.\n")
            for pwm in pwms:
                self.abc4pwm2jaspar(pwm, out)

        elif in2out == 'jaspar2abc4pwm':
            print("\nConverting from JASPAR to abc4pwm.\n")
            for pwm in pwms:
                self.jaspar2abc4pwm(pwm, out)
        else:
            "Wrong format specified, please assign value from the following. \n" \
            "--in2out 'abc4pwm2transfac' \n " \
            "--in2out 'transfac2abc4pwm'"
            exit()
        print("\nConversion task finished.\n")


    def abc4pwm2transfac(self, pwm, out):
        """

        :param pwm: position weight matrix
        :param out: path where you want to save converted to Transfac
        :return: output path of the file
        """

        if not os.path.exists(out):
            os.makedirs(out, exist_ok=True)

        output_file_path = os.path.join(out, str(pwm).split("/")[-1])
        with open(output_file_path, 'w') as f:

            matrix, matrix_string, sum, max, info  = read_energy_matrix(pwm)
            f.writelines("ID \t" + str(info[0]).split(' ')[0])
            f.writelines("BF \t " + str(info[0]).split(' ')[0])
            f.writelines("P0\tA\tC\tG\tT\n")
            for i, line in enumerate(matrix_string):
                str_line = ('\t').join(
                    ["{:02d}".format(i) , ('\t').join([str("{:.3f}".format(x)) for x in matrix[i, :]]), matrix_string[i]]) + '\n'
                f.write((str_line))
            f.writelines('XX\n//\n')
        return output_file_path

    def transfac2abc4pwm(self, pwm, out):
        """

        :param pwm: position weight matrix
        :param out: path where you want to save converted to abc4pwm format
        :return: output path of the file
        """

        if not os.path.exists(out):
            os.makedirs(out, exist_ok=True)

        output_file_path = os.path.join(out, str(pwm).split("/")[-1])
        with open(output_file_path, 'w') as f:

            matrix, matrix_string, sum, max, info  = read_energy_matrix_transfac(pwm)

            f.writelines(info[0].split(' ')[1:])
            f.writelines("dummy \t # p-value #File converted from TRANSFAC\n")
            f.writelines("dummy \t # benferroni\n")
            f.writelines("dummy 0 1 0 0 0 1 1 1 1 \t #derived from leading strand\n")

            f.writelines("#\tA\tC\tG\tT\n")
            for i, line in enumerate(matrix_string):
                str_line = ('\t').join(
                    [matrix_string[i], ('\t').join([str("{:.3f}".format(x)) for x in matrix[i, :]])]) + '\n'
                f.write((str_line))
        return output_file_path

    def abc4pwm2jaspar(self, pwm, out):
        """

        :param pwm: position weight matrix
        :param out: path where you want to save converted to Transfac
        :return: output path of the file
        """

        if not os.path.exists(out):
            os.makedirs(out, exist_ok=True)

        output_file_path = os.path.join(out, str(pwm).split("/")[-1])
        with open(output_file_path, 'w') as f:

            matrix1, matrix_string, sum, max, info  = read_energy_matrix(pwm)
            f.writelines(str(info[0]).split(' ')[0])
            matrix = matrix1.transpose()
            bases = ['A', 'C', 'T', 'G']
            for i, line in enumerate(matrix):
                str_line = ('\t').join(
                    [bases[i] , '[' + ('\t').join([str("{:.3f}".format(x)) for x in matrix[i, :]]) + ']']) + '\n'
                f.write((str_line))
        return output_file_path

    def jaspar2abc4pwm(self, pwm, out):
        """

        :param pwm: position weight matrix
        :param out: path where you want to save converted to abc4pwm format
        :return: output path of the file
        """

        if not os.path.exists(out):
            os.makedirs(out, exist_ok=True)

        output_file_path = os.path.join(out, str(pwm).split("/")[-1])
        with open(output_file_path, 'w') as f:

            matrix, matrix_string, sum, max, info  = read_energy_matrix_jaspar(pwm)

            f.writelines(info[0].split(' ')[1:])
            f.writelines("dummy \t # p-value #File converted from jaspar to abc4pwm\n")
            f.writelines("dummy \t # benferroni\n")
            f.writelines("dummy 0 1 0 0 0 1 1 1 1 \t #derived from leading strand\n")

            f.writelines("#\tA\tC\tG\tT\n")
            for i, line in enumerate(matrix):
                str_line = ('\t').join(
                    ["X", ('\t').join([str("{:.3f}".format(x)) for x in matrix[i, :]])]) + '\n'
                f.write((str_line))
        return output_file_path

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


if __name__ == '__main__':

    # input_path = "../data/in/in_pwms/"
    # Conversion(input_path, 'abc4pwm2transfac' ,"../data/out/convert/to_transfac")

    input_path = "../data/out/convert/to_transfac/"
    Conversion(input_path, 'transfac2abc4pwm' ,"../data/out/convert/to_abc")

    # input_path = "../data/in/in_pwms/"
    # Conversion(input_path, 'abc4pwm2jaspar' ,"../data/out/convert/to_jaspar")

    # input_path = "../data/out/convert/to_jaspar/"
    # Conversion(input_path, 'jaspar2abc4pwm' ,"../data/out/convert/to_abc")