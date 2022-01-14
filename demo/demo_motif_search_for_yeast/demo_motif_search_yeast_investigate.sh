#!/usr/bin/env bash

abc4pwm ensemble_investigate --path_to_predicted_files ../../data/out/ensemble/swi4/ --db_folder ../../data/in/SGD/SGD_yeast_pwm/ --output_folder ../../data/out/ensemble/search_out/ --db_type folder --input_count 1 --db_count 0 --db_file_type '.txt' --input_file_type '.mlp' --input_prob 0 --db_prob 1 --tf_name SWI4