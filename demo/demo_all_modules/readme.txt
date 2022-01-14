
============================demos of Affinity Based Clustering for Position weight Matrices============================


A new tool for the clustering of Transcription Factors is designed. Details are as follows

usage: abc4pwm [-h] {cleandatabase_for_classification, classification, clustering, representative_motif, quality_assessment, visualize, plot_cluster_motifs, text_tfdb}

positional arguments:
{cleandatabase_for_classification,classification,clustering,representative_motif,quality_assessment,visualize,plot_cluster_motifs,text_tfdb}                        


Demo of tasks:

------------------------------------cleandatabase_for_classification----------------------------------
./1demo_clean_database.sh


This demo creates a clean database for input pwms given in ../data/in/in_pwms/

----------------------------------classification----------------------------------
./2demo_classification.sh

This demo classify the input pwms in their respective DNA binding domains



----------------------------------clustering----------------------------------

./3demo_clustering.sh

This demo cluster the pwms that were classify into their respective DBDs and produce results into ../data/out/clustering_out/


----------------------------------quality_assessment----------------------------------
./4demo_quality_assesment.sh

This demo evaluate the clusters produce in the clustering step. Stores evaluated results in ../data/out/quality_assesment_out/

----------------------------------representative_motif------------------------------------
./5demo_representative.sh

This demo produce representatives of the all clusters in HMG DNA binding domain.

----------------------------------plot_cluster_motifs----------------------------------
./6demo_plot_cluster.sh
 
This demo plots the motifs of all the clusters in HMG DNA binding domain.


---------------------------------Conversions----------------------------------
./7demo_conversion.sh

This demo convert file formats of abc4pwm to transfac and jasper and vice versa.

----------------------------------searching----------------------------------
./8demo_searching.sh

This demo search for a pwm ../data/in/in_pwms/AP1_known1_from_AP1_2_from_AP-1_transfac_M00172.mlp in database ../data/in/in_pwms/ and extract the top 5 matches with the input pwm.

----------------------------------text_tfdb---------------------------------
./9demo_text_tfdb.sh

This demo produces a database for input pwms in a text file which can be used for different purposes.


----------------------------------ensemble learning----------------------------
./10demo_ensemble.sh

This demo find motifs from expression and sequence files in ../data/in/expression and ../data/in/seq by randomly select 200 expressions and repeat this process 15 times to predict pwms and use BayesPI program.


----------------------------------ensemble investigate--------------------------
./11demo_investigate.sh

This demo takes predicted pwms as inputs. It run clustering and eliminates small cluster (which are identified as a parameter) , then make representative of each cluster and search a database against each of the pwm to give results.

