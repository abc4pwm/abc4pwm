#!/bin/sh
# optimized parameters: -strand=0 -min_L=20, -max_L=20 -max_iteration= 1500 -initial_beta=10 -normalize=4

echo 'run motif find ...'

# Give a name for the out fold:
out='out_motif_Find'

echo ${out}
mkdir ${out}
#../bayesPI -information_content_threshold=0.2 -seed=1 -split_sequences=1 -initial_beta=10 -dependence=0 -strand=0 -max_iteration=1500 \
#        -max_loop=3 -normalize=4 \
#	-out=${out} \
#        -exp=data_in/TIME0_er1_rmDup_narrowPeaks_EncodeBlackListFiltered.bed_bsites_Tag \
#        -seq=data_in/TIME0_er1_rmDup_narrowPeaks_EncodeBlackListFiltered.bed_bsites_Seq_600_bp \
#        -min_L=20 -max_L=20 -p_value=0.001
#


if [[ "$OSTYPE" == "linux-gnu"* ]]; then

        ../../../abc4pwm/bin/Linux/bayesPI -information_content_threshold=0.2 -seed=1 -split_sequences=1 -initial_beta=10 -dependence=0 -strand=0 -max_iteration=1500 \
        -max_loop=3 -normalize=4 \
	-out=${out} \
        -exp=data_in/TIME0_er1_rmDup_narrowPeaks_EncodeBlackListFiltered.bed_bsites_Tag \
        -seq=data_in/TIME0_er1_rmDup_narrowPeaks_EncodeBlackListFiltered.bed_bsites_Seq_600_bp \
        -min_L=20 -max_L=20 -p_value=0.001

elif [[ "$OSTYPE" == "darwin"* ]]; then
        # Mac OSX
        ../../../abc4pwm/bin/Mac/bayesPI -information_content_threshold=0.2 -seed=1 -split_sequences=1 -initial_beta=10 -dependence=0 -strand=0 -max_iteration=1500 \
        -max_loop=3 -normalize=4 \
	-out=${out} \
        -exp=data_in/TIME0_er1_rmDup_narrowPeaks_EncodeBlackListFiltered.bed_bsites_Tag \
        -seq=data_in/TIME0_er1_rmDup_narrowPeaks_EncodeBlackListFiltered.bed_bsites_Seq_600_bp \
        -min_L=20 -max_L=20 -p_value=0.001

else
    echo 'This Operating System is not supported.'
fi


