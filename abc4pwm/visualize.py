import matplotlib.pyplot as plt
import os
import pandas as pd



class Visualize():
    def __init__(self, path_to_folder_of_assessment_file, path_to_folder_of_DBDs, output_folder, dbd_for_plot, task='boxplot'):
        """

        :param path_to_folder_of_assessment_file: folder path from where assesment file should be taken
        ;param path_to_folder_of_DBDs: this folder should contain clustered DBD folders
        :param output_folder: folder where visualization output should go
        :param dbd_for_plot: path of dbd which is needed to be visualized
        :param task: which visualization task needed. e.g., boxplot etc
        """

        if os.path.exists(os.path.join(path_to_folder_of_assessment_file, 'quality_df.json')):
            quality_df = pd.read_json(os.path.join(path_to_folder_of_assessment_file, 'quality_df.json'))
        else:
            print("Quality assesment file is not available. Please generate assessment by quality_cluster")
            exit()

        if 'boxplot' in task:
            if 'all' in dbd_for_plot:
                for x in sorted(os.listdir(path_to_folder_of_DBDs)):
                    if x.startswith("."):
                        continue
                    dbd_for_plot = os.path.join(path_to_folder_of_DBDs, x)
                    dbd_for_plot = dbd_for_plot + '/'
                    self.print_boxplot_dbd(quality_df, dbd_for_plot, output_folder)
            else:
                self.print_boxplot_dbd(quality_df, dbd_for_plot, output_folder)
        print("You can see plots in ", output_folder)

    @staticmethod
    def print_boxplot_dbd(quality_df, plot_dbd, output_path):
        if len(os.listdir(plot_dbd)) == 0:
            return 1
        dbd_to_plot = str(plot_dbd).split('/')[-2]

        selected_dbd = quality_df.loc[quality_df['dbd'] == dbd_to_plot]
        similarityMatrices = []
        cluster_for_plotting = []

        for row in range(len(selected_dbd)):
            similarityMatrices.append(selected_dbd.iloc[row, 9])
            cluster_for_plotting.append(str(row) + '(' + str(selected_dbd.iloc[row, 7]) + ')')

        labels = cluster_for_plotting
        if not len(labels):
            return 1


        fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(24, 12))

        # rectangular box plot
        ax1.boxplot(similarityMatrices,
                             vert=True,  # vertical box alignment
                             patch_artist=True,  # fill with color
                             labels=labels)  # will be used to label x-ticks
        ax1.set_title('Rectangular box plot for ' + str(dbd_to_plot))

        # adding horizontal grid lines
        for ax in [ax1]:
            ax.yaxis.grid(True)
            ax.set_xlabel('Clusters')
            ax.set_ylabel('Observed Similarities')

        plt.xticks(rotation=45)
        path_to_boxplots = output_path
        if not os.path.exists(path_to_boxplots):
            os.makedirs(path_to_boxplots, exist_ok=True)
        plt.savefig(os.path.join(path_to_boxplots, dbd_to_plot + '_boxplot.png'))
        plt.close()

if __name__ == "__main__":
    Visualize('../data/in/','../data/out/clustering_out/','../data/out/plots/boxplots/', 'all')
