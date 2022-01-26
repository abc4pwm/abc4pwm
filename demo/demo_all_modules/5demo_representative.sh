#!/usr/bin/env bash


abc4pwm representative_motif --path_to_clusters ../../data/out/quality_assessed_out/HMG/out/ --clusters all --ic 0.45




#abc4pwm ensemble_investigate --qa 1 --dst_for_bad_pwms ../../data/out/ensemble_new/bad_pwms/ \
#        --path_to_predicted_files ../../data/out/ensemble_new/er_er/ \
#        --db_folder ../../data/in/in_pwms/ --output_folder ../../data/out/ensemble_new/search_out/ \
#        --db_type folder --input_count 1 --db_count 1 --db_file_type '.mlp' \
#        --input_file_type '.mlp' --input_prob 0 --db_prob 0 --min_pwms_in_cluster 10 --top_n 2 \
#        --seed 1 --ic_for_rep 0 --mean_threshold 0.5 --z_score_threshold 0.3 --damp 0.7 --max_iter 1000 --convergence_iter 30 \
#        --top_occurrences 0.5 --occurrences_threshold 0.5