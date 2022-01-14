#!/usr/bin/env bash

#for format conversion
#abc4pwm to transfac
abc4pwm conversion --pwm_files_directory ../../data/in/in_pwms/ --in2out abc4pwm2transfac --output_folder ../../data/out/convert/to_transfac/


#transfac to abc4pwm
abc4pwm conversion --pwm_files_directory ../../data/out/convert/to_transfac/ --in2out transfac2abc4pwm --output_folder ../../data/out/convert/to_abc/

#abc4pwm to jaspar
abc4pwm conversion --pwm_files_directory ../../data/in/in_pwms/ --in2out abc4pwm2jaspar --output_folder ../../data/out/convert/to_jaspar/