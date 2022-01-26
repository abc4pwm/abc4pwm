#!/usr/bin/env bash


abc4pwm quality_assessment --dbd_folders_directory ../../data/out/Comparison/in_abc4pwm_motif_clustering \
                            --out_path_for_qa_clusters ../../data/out/Comparison/quality_assessed_out_abc4pwm_new/ \
                            --output_folder_for_text_report ../../data/out/Comparison/abc4pwm/reports_in_text/ \
                             --output_path_for_quality_assessment_file ../../data/in/ --load_new_assesment 1


abc4pwm visualize   --path_to_folder_of_assessment_file ../../data/in/ \
                    --path_to_folder_of_DBDs ../../data/out/Comparison/quality_assessed_out_abc4pwm_new/ \
                    --output_folder ../../data/out/Comparison/plots/abc4pwm/boxplots/ \
                    --dbd_for_plot all

abc4pwm quality_assessment --dbd_folders_directory ../../data/out/Comparison/stamp_in_abc4pwm_motif \
                            --out_path_for_qa_clusters ../../data/out/Comparison/quality_assessed_out_stamp/ \
                            --output_folder_for_text_report ../../data/out/Comparison/stamp/reports_in_text/ \
                             --output_path_for_quality_assessment_file ../../data/in/ --load_new_assesment 1



abc4pwm visualize   --path_to_folder_of_assessment_file ../../data/in/ \
                    --path_to_folder_of_DBDs ../../data/out/Comparison/quality_assessed_out_stamp/ \
                    --output_folder ../../data/out/Comparison/plots/stamp/boxplots/ \
                    --dbd_for_plot all

