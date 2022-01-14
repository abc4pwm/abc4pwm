
============================Affinity Based Clustering for Position Weight Matrices============================


A new tool for the clustering of Transcription Factors is designed that work on RNAseq, ChIPseq, Human Transfaction factors,
and other biological datasets to give easy interpretition for biologists. Details are as follows

usage: abc4pwm [-h] {cleandatabase_for_classification, classification, clustering, representative_motif, quality_assessment, visualize, plot_cluster_motifs, text_tfdb, searching,conversion,ensemble_learning,ensemble_investigate}


Tasks available for using:
    cleandatabase_for_classification
    classification
    clustering
    representative_motif
    quality_assessment
    visualize
    plot_cluster_motifs
    text_tfdb
    conversion
    ensemble_learning
    ensemble_investigate
    conversion
    searching


optional arguments:
  -h, --help            show this help message and exit

You should run the pipeline in the following order to generate all necessary files:
1- 'cleandatabase_for_classification' :In the first step you should generate clean database to classify inputs files
2- 'classification' :In this step classify files according to their DBD by adding labels
3- 'clustering ' :In the step apply clustering on output folder produced by step 2 
4- 'quality_assessment' :This step calculates quality of clusters and remove bad quality pwms from clusters 
5- 'representative_Motif': This step prepare a representative motif of the cluster
6- 'plot_cluster_motifs': This will plot the motifs inside cluster in along with the names in a pdf file. Also generate a text report.

Example of tasks:

------------------------------------cleandatabase_for_classification----------------------------------

usage: abc4pwm cleandatabase_for_classification [-h]
                                                      [--pwm_files_directory FOLDER]
                                                      [--read_new NUMBER]


Arguments:
  --pwm_files_directory FOLDER       This folder should contain the input pwm files in . mlp formate


optional arguments:
  --read_new NUMBER     (Optional) Select 1 if you want a new updated read from sources(internet Connection Required). Default 0
  -h, --help            show this help message and exit



----------------------------------classification----------------------------------

usage: abc4pwm classification [-h] [--pwm_files_directory FOLDER]
                                        [--output_directory FOLDER]
                                        [--original_pwm_files_directory FOLDER]
                                        [--load_new_db NUMBER]


Arguments:
  --pwm_files_directory             FOLDER      This folder should contain the input pwm files in . mlp formate
  --output_directory                FOLDER      This folder should point to folder where files should go after getting labels of DBD
  --original_pwm_files_directory    FOLDER      This folder should contain a copy of input pwm files in . mlp formate

optional arguments:
  --load_new_db                     NUMBER      (Optional) This should be 1 if you want to download a new update db from sources 0 by default
  -h, --help            show this help message and exit


----------------------------------clustering----------------------------------


usage: abc4pwm clustering [-h] [--dbd_folders_directory FOLDER]
                                    [--output_directory FOLDER]
                                    [--in_dbd NUMBER]
                                    [--minimum_pwms_in_dbd NUMBER]
                                    [--max_processors NUMBER]
                                    [--seed Number]
                                    [--damp Number]
                                    [--max_iter Number]
                                    [--convergence_iter Number] [--preference Number]


Arguments:
--dbd_folders_directory 	FOLDER	This folder should contain DBD folders. Output of classifcation_pwm should be this task input
  --output_directory 		FOLDER	This folder should to output folder. Output of classifcation_pwm should be input to this task


optional arguments:
  -h, --help                show this help message and exit
  --in_dbd 			        NUMBER      This should be 0 if you want to cluster all together (non DBD) 1 by default
  --minimum_pwms_in_dbd 	NUMBER		minimum number of pwms in a dbd to be clusteredDefault value is 5
  --max_processors 		    NUMBER		maximum number of processors for paralle processingDefault value is 5
  --seed                    Number      Seed for random selection of cluster center. Input 1 to fix. Default Seed is 0.
  --damp                    Number      Damping factor (between 0.5 and 1) is the extent to which the current value is maintained relative to incoming values (weighted 1 - damping).
                                        This in order to avoid numerical oscillations when updating these values (messages).
  --max_iter                Number      Maximum number of iterations.
  --convergence_iter        Number      Number of iterations with no change in the number of estimated clusters that stops the convergence.
  --preference              Number      Preferences for each point - points with larger values of preferences are more likely to be chosen as exemplars.
                                        The number of exemplars, ie of clusters, is influenced by the input preferences value. If the preferences are not passed as arguments,
                                        they will be set to the median of the input similarities.


----------------------------------quality_assessment----------------------------------

usage: abc4pwm quality_assessment [-h] [--dbd_folders_directory FOLDER]
                                     [--output_directory FOLDER]
                                     [--dbd_for_plotting FOLDER]
                                     [--load_new_assesment <class 'bool'>]
                                     [--mean_threshold NUMBER]
                                     [--z_score_threshold NUMBER]
                                     [--top_occurrences NUMBER]
                                     [--occurrences_threshold NUMBER]


Arguments:
  --dbd_folders_directory 			FOLDER	This folder should contain clustered DBD folders. Output of clusterings should be this task input
  --out_path_for_qa_clusters 			FOLDER	This folder should point to output folderwhere quality assessed clusters will be stored.
  --output_folder_for_text_report 		FOLDER	Specify a folder where report in txt file will be stored.

optional arguments:
  -h, --help            show this help message and exit
  --output_path_for_quality_assessment_file 	FOLDER	(Optional) folder where quality assessment .json file will be stored.Default, data/in/
  --load_new_assesment 				NUMBER	 1 if you want to do new assesment, 0 if you want load existing assesment
  --mean_threshold 				    NUMBER	 mean threshold for uncertain clustersDefault value is 0.80
  --z_score_threshold 				NUMBER	 max negative threshold of zscore for similarity values of pwmDefault value is -1.0
  --top_occurrences 				NUMBER	 This value corresponds to occurrence of a pwm less than a threshold z-scoreValue is between 0 to 1. Default value is 0.15
  --occurrences_threshold 			NUMBER	 This value corresponds to threshold of occurrence from top occurrencesValue is between 0 to 1. Default value is 0.05


----------------------------------representative_motif------------------------------------
usage: abc4pwm representative_motif [-h] [--path_to_clusters FOLDER]
                                          [--dbd string] 
					                      [--clusters string]
                                          [--ic NUMBER]
                                          [--best_match_initial_motif NUMBER]
                                          [--mean_threshold NUMBER]
                                          [--z_score_threshold NUMBER]
                                          [--top_occurrences NUMBER]
                                          [--occurrences_threshold NUMBER]

Arguments:
  --path_to_clusters 		FOLDER		This folder should contain the clusters folders
  --clusters 			string     	This argument should be cluster numbers as string For example, 0,1,2,3,4 if you want to plot these clusterswrite 'all' if you want to make 						representative for all clusters 
optional arguments:
  -h, --help            show this help message and exit

      --dbd 			            string      Default value is 'selected'. Representative of clusters pathmentioned in --path_to_clusters parameter will be calculated.  						Write 'all' if representative calculation of all dbd and all clusters is required.
      --best_match_initial_motif 	NUMBER		This should be 0 if you want initial motif to be random 1 by default
      --mean_threshold	 	        NUMBER		mean threshold for uncertain clustersDefault value is 0.80
      --z_score_threshold	 	    NUMBER		max negative threshold of zscore for similarity values of pwmDefault value is -1.0
      --top_occurrences 		    NUMBER		This value corresponds to occurrence of a pwm less than a threshold z-scoreValue is between 0 to 1. Default value is 0.15
      --occurrences_threshold 	    NUMBER		This value corresponds to threshold of occurrence from top occurrencesValue is between 0 to 1. Default value is 0.05
      --ic			                NUMBER      Information Content for trimming edges.Default value is 0.4

----------------------------------plot_cluster_motifs----------------------------------
usage: abc4pwm plot_cluster_motifs [-h] 
				[--path_to_clusters FOLDER]
                [--output_folder FOLDER]
				[--clusters string]
				[--dbd string]

Arguments:
  --path_to_clusters    FOLDER      This folder should contain clusters of a DBD, Write all if you want to plot all dbds.
  --output_folder       FOLDER      This folder should contain clusters of a DBD
  --clusters            string      This arguement should be cluster numbers as string For example, 0,1,2,3,4 if you want to plot these clusters. write 'all' if you want to plot all clusters
optional arguments:
  -h, --help            show this help message and exit

  --dbd 		string       Default value is 'selected'. Clusters pathmentioned in --path_to_clusters parameter will be printed  Write 'all' if plotting of all dbd and all clusters is required.



---------------------------------visualize----------------------------------
usage: abc4pwm visualize [-h]
                               [--path_to_folder_of_assessment_file FOLDER]
                               [--path_to_folder_of_DBDs FOLDER]
                               [--output_folder FOLDER]
                               [--dbd_for_plot FOLDER] [--task string]

Arguments:
  --path_to_folder_of_assessment_file 	FOLDER		folder path from where quality assessment file should be taken
  --path_to_folder_of_DBDs 		        FOLDER		this folder should contain clustered DBD folders
  --output_folder 			            FOLDER		folder where visualization output should be saved
  --dbd_for_plot 			            FOLDER		path of dbd which is needed to be visualized. Write 'all' if you want to plot for all dbds.


optional arguments:
  -h, --help            show this help message and exit
  --task 				string         Specify visualization taskFor example, boxplot, pichart, etc Default is boxplot.





----------------------------------text_tfdb---------------------------------
usage: abc4pwm text_tfdb [-h]   [--pwm_files_directory FOLDER]
                                [--output_directory FOLDER]

optional arguments:
  -h, --help            show this help message and exit
Arguments:
  --pwm_files_directory FOLDER      This folder should contain the input pwm files in . mlp format
  --output_directory      FOLDER      This folder should to output folder. Boxplot will go to this folder

----------------------------------searching----------------------------------
usage: abc4pwm searching [-h] 
			[--pwm file]
			[--db_path path or list]
            [--output_directory FOLDER]
			[--db_type string]
            [--db_format string]
			[--top_n NUMBER]
		    [--tf_name string]
            [--input_count NUMBER]
            [--db_count NUMBER]
            [--db_file_type String]
            [--input_file_type String]
            [--input_prob Number]
            [--db_prob Number]


Arguments:
  --pwm 		        file            position weight matrix file (motif) which you want to search. .mlp format
  --db_path 		    path or list	path to clustered dbds according to the hierarchy of abc4pwm.if db_type=list then then this paramter should be list of pwms, against whcih you are searching the pwm
  --output_directory 	FOLDER		    This folder should to output folder. search_result.pdf output file will be stored here.

optional arguments:
  -h, --help            show this help message and exit
  --tf_name             string      If you want to search specific tf in a folder then use this parameter
  --db_type 		    string      If database for comparison is folder hierarchy like abs4pwm then this willbe db_type=path. Write list if providing list of pwms for comparison
  --db_format 		    string    	If database for comparison have format, please mention.Supported formats are abc4pwm, Tranfac, Jaspar. default is abc4pwm
  --top_n 		        NUMBER      Number of top matches from the database. Default 5
  --input_count         NUMBER      1 if input file contains values in counts
  --db_count            NUMBER      1 if database file contains values in counts
  --db_file_type        String      mention the extension of file type. e.g., .mlp, .txt
  --input_file_type     String      mention the extension of file type. e.g., .mlp, .txt
  --input_prob          Number      1 if input file contains values in probabilities
  --db_prob             Number      1 if databas file contains values in probabilities
----------------------------------conversion----------------------------------
usage: abc4pwm conversion [-h]
                    [--pwm_files_directory FOLDER]
                    [--in2out string]
                    [--output_folder FOLDER]

optional arguments:
  -h, --help                            show this help message and exit
  --pwm_files_directory     FOLDER      This folder should contain the input pwm files in . mlp formate  which you want to convert
  --in2out                  string      Specify conversion like the following.
                                        --in2out 'abc4pwm2transfac'
                                        --in2out 'transfac2abc4pwm'
  					--in2out 'abc4pwm2jaspar'
                                        --in2out 'jaspar2abc4pwm'
  --output_folder           FOLDER      folder where converted files should be saved

----------------------------------ensemble learning----------------------------
usage: abc4pwm ensemble_learning [-h] [--opt_dependence Number]
                                 [--numP Number] [--opt_numOfWeakReads Number]
                                 [--number_of_genes Number]
                                 [--expFile File path]
                                 [--opt_weak_expFile File path]
                                 [--opt_seqFile File path]
                                 [--opt_out File path.] [--opt_loops Number]
                                 [--opt_min_L Number] [--opt_max_L Number]
                                 [--opt_iteration Number]
                                 [--opt_p_value float] [--opt_strand Number]
                                 [--opt_normalization Number]
                                 [--max_processors Number]


Arguments:
  --expFile File        path        Strong Expression file path
  --opt_seqFile         path        Strong sequence file path.



optional arguments:
  -h, --help                        show this help message and exit
  --opt_dependence      Number      Define dependence, Default is 0.
  --numP                Number      number of times random selections will be done. Default is 15
  --opt_numOfWeakReads  Number      number of weak read if any. Default 0.
  --number_of_genes     Number      number of genes to be selected from input. Default is 200
  --opt_weak_expFile    File path   Weak expression file path.
  --opt_out             File path.  output folder path.
  --opt_loops           Number      Number of loops to repeat calculations. Default is 3.
  --opt_min_L           Number      Minimum length for predicted pwm (motif), Default is 9.
  --opt_max_L           Number      Maximum length for predicted pwm (motif). Default is 9.
  --opt_iteration       Number      Number of iterations. Default is 500.
  --opt_p_value         float       p value. Default is 0.0001
  --opt_strand          Number      Strand.     Default is 0
  --opt_normalization   Number      Normalization value. Default is 2.
  --max_processors      Number      Define maximum number of processors for parallel computing. Default is mp.cpu_count()
  --seed Number         Number      Seed for random selection. Default seed is 0(random).


----------------------------------ensemble investigate--------------------------
usage: abc4pwm ensemble_investigate [-h]
                                    [--path_to_predicted_files . mlp file]
                                    [--db_folder path or list]
                                    [--output_folder FOLDER]
                                    [--tf_name string]
                                    [--db_type string]
                                    [--top_n Number]
                                    [--dst_for_bad_pwms path]
                                    [--mean_threshold NUMBER]
                                    [--z_score_threshold NUMBER]
                                    [--top_occurrences NUMBER]
                                    [--occurrences_threshold NUMBER]
                                    [--ic_for_rep NUMBER]
                                    [--min_pwms_in_cluster Number]
                                    [--db_format string]
                                    [--input_count NUMBER]
                                    [--db_count NUMBER]
                                    [--db_file_type String]
                                    [--input_file_type String]
                                    [--input_prob Number]
                                    [--db_prob Number]
                                    [--qa Number]
                                    [--seed Number]
                                    [--damp Number]
                                    [--max_iter Number]
                                    [--convergence_iter Number]
                                    [--preference Number]

Arguments:
  --path_to_predicted_files     .mlp files      Folder which contain predicted files.
  --db_folder                   path or list    path to clustered dbds according to the hierarchy of abc4pwm.if db_type=list then then this paramter should be list of pwms, against whcih you are searching the pwm
  --output_folder               FOLDER          This folder should to output folder. search_result.pdf output file will be stored here.


optional arguments:
  -h, --help            show this help message and exit
  --tf_name                     String          If you want to search specific tf in a folder then use this parameter
  --db_type                     String          If database for comparison is folder hierarchy like abs4pwm then this will be db_type=path. Write list if providing list of pwms for comparison
  --top_n                       Number          Specify how many top matches from database is required. Default 2
  --min_pwms_in_cluster         Number          Number of minimum acceptable pwms in a cluster made from predicted pwms. Default 3
  --db_format                   string          If database for comparison have format, please mention.Supported formats are abc4pwm, Tranfac, Jaspar. default is abc4pwm
  --input_count                 NUMBER          1 if input file contains values in counts
  --db_count                    NUMBER          1 if database file contains values in counts
  --db_file_type                String          mention the extension of file type. e.g., .mlp, .txt
  --input_file_type             String          mention the extension of file type. e.g., .mlp, .txt
  --input_prob                  Number          1 if input file contains values in probabilities
  --db_prob                     Number          1 if database file contains values in probabilities
  --ic			                NUMBER          Information Content for trimming edges.Default value is 0.4
  --qa                          Bool            1 if quality assessment of predicted files (clusters) is need, 0 is by default.
  --seed                        Number          Seed for random selection of cluster center. Input 1 to fix. Default Seed is 0.
  --damp                        Number          Damping factor (between 0.5 and 1) is the extent to which the current value is maintained relative to incoming values (weighted 1 - damping).
                                                This in order to avoid numerical oscillations when updating these values (messages).
  --max_iter                    Number          Maximum number of iterations.
  --convergence_iter            Number          Number of iterations with no change in the number of estimated clusters that stops the convergence.
  --preference                  Number          Preferences for each point - points with larger values of preferences are more likely to be chosen as exemplars.
                                                The number of exemplars, ie of clusters, is influenced by the input preferences value. If the preferences are not passed as arguments,
                                                they will be set to the median of the input similarities.


----------------------------Example-------------------
Examples of all modules are given in demo folder
