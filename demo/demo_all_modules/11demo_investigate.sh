#!/usr/bin/env bash

abc4pwm ensemble_investigate --path_to_predicted_files ../../data/out/ensemble/swi4/ --db_folder ../../data/in/SGD/SGD_yeast_pwm/ --dst_for_bad_pwms ../../data/out/ensemble/swi4/uncertain_pwms \
    --output_folder ../../data/out/ensemble/swi4/search_out/ --db_type 'folder' --input_count 1 --db_count 1 --db_file_type '.txt' \
    --input_file_type '.mlp' --input_prob 0 --db_prob 1 --min_pwms_in_cluster 3 --top_n 2 --qa 1 --ic_for_rep 0 \
    --mean_threshold 0.7 --z_score_threshold -1.5 --top_occurrences 0.15 --occurrences_threshold 0.15 \
    --damp 0.8 --max_iter 500 --seed 1