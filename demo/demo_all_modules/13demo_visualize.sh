#!/usr/bin/env bash

#this creates a uniformly named database of pwms from three different sources.

abc4pwm visualize   --path_to_folder_of_assessment_file ../../data/in/ \
                    --path_to_folder_of_DBDs ../../data/out/Comparison/quality_assessed_out_abc4pwm_old/ \
                    --output_folder ../../data/out/Comparison/plots/abc4pwm/boxplots/ \
                    --dbd_for_plot all