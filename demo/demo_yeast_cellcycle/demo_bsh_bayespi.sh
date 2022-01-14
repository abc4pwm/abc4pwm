#!/usr/bin/env bash

#For strand1
#see log in out.log.alpha files.

if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    platform='Linux'

elif [[ "$OSTYPE" == "darwin"* ]]; then
    # Mac OSX
    platform='Mac'

else
    echo 'This Operating System is not supported.'
fi


echo ${platform}

(./../../abc4pwm/bin/${platform}/bayesPI -max_loop=6 -seed=1 -information_content_threshold=0.2 -max_evidence=3 \
-dependence=0 -strand=1 -max_iteration=1200 -normalize=2 -min_L=6 -max_L=14 -p_value=0.05 \
-out=yeast_out/strand1/cycle_alpha7/ -exp=yeast_in/yeast_GE_7.txt -seq=yeast_in/yeast_GE_800bp_fasta.fa
 ) 2>out.log.alpha7

(../../abc4pwm/bin/${platform}/bayesPI -max_loop=6 -seed=1 -information_content_threshold=0.2 -max_evidence=3 \
-dependence=0 -strand=1 -max_iteration=1200 -normalize=2 -min_L=6 -max_L=14 -p_value=0.05 \
-out=yeast_out/strand1/cycle_alpha42/ -exp=/yeast_in/yeast_GE_42.txt -seq=yeast_in/yeast_GE_800bp_fasta.fa
 ) 2>out.log.alpha42


(../../abc4pwm/bin/${platform}/bayesPI -max_loop=6 -seed=1 -information_content_threshold=0.2 -max_evidence=3 \
-dependence=0 -strand=1 -max_iteration=1200 -normalize=2 -min_L=6 -max_L=14 -p_value=0.05 \
-out=yeast_out/strand1/cycle_alpha49/ -exp=/yeast_in/yeast_GE_49.txt -seq=yeast_in/yeast_GE_800bp_fasta.fa
 ) 2>out.log.alpha49


#for strand0

(../../abc4pwm/bin/${platform}/bayesPI -max_loop=6 -seed=1 -information_content_threshold=0.2 -max_evidence=3 \
-dependence=0 -strand=0 -max_iteration=1200 -normalize=2 -min_L=6 -max_L=14 -p_value=0.05 \
-out=yeast_out/strand1/cycle_alpha7/ -exp=/yeast_in/yeast_GE_7.txt -seq=yeast_in/yeast_GE_800bp_fasta.fa
 ) 2>out.log.alpha7

(../../abc4pwm/bin/${platform}/bayesPI -max_loop=6 -seed=1 -information_content_threshold=0.2 -max_evidence=3 \
-dependence=0 -strand=0 -max_iteration=1200 -normalize=2 -min_L=6 -max_L=14 -p_value=0.05 \
-out=yeast_out/strand1/cycle_alpha42/ -exp=/yeast_in/yeast_GE_42.txt -seq=yeast_in/yeast_GE_800bp_fasta.fa
 ) 2>out.log.alpha42


(../../abc4pwm/bin/${platform}/bayesPI -max_loop=6 -seed=1 -information_content_threshold=0.2 -max_evidence=3 \
-dependence=0 -strand=0 -max_iteration=1200 -normalize=2 -min_L=6 -max_L=14 -p_value=0.05 \
-out=yeast_out/strand1/cycle_alpha49/ -exp=/yeast_in/yeast_GE_49.txt -seq=yeast_in/yeast_GE_800bp_fasta.fa
 ) 2>out.log.alpha49


