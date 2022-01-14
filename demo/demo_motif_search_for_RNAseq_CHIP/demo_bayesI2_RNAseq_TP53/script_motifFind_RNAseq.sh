#!/bin/sh
# Optimized parameters: -strand=0 -min_L=10 -max_L=20 -max_iteration=1000  -initial_beta=1

echo 'run motif find...'

out='out_HCT116_motif_Find'
echo ${out}
mkdir ${out}

#../bayesPI -information_content_threshold=0.2 -seed=1 -split_sequences=1 \
#	-initial_beta=1 -dependence=0 \
#	-strand=0 -max_iteration=1000 \
#	-min_L=10 -max_L=20 -p_value=0.001 \
#        -max_loop=3 -normalize=0 \
#	-out=${out} \
#        -exp=data_in/HCT116_nutlin_S2_RNAseq_data_log2.txt \
#        -seq=data_in/HCT116_nutlin_S2_RNAseq_data_log2_seq.fasta



if [[ "$OSTYPE" == "linux-gnu"* ]]; then

        ../../../abc4pwm/bin/Linux/bayesPI -information_content_threshold=0.2 -seed=1 -split_sequences=1 \
	-initial_beta=1 -dependence=0 \
	-strand=0 -max_iteration=1000 \
	-min_L=10 -max_L=20 -p_value=0.001 \
        -max_loop=3 -normalize=0 \
	-out=${out} \
        -exp=data_in/HCT116_nutlin_S2_RNAseq_data_log2.txt \
        -seq=data_in/HCT116_nutlin_S2_RNAseq_data_log2_seq.fasta

elif [[ "$OSTYPE" == "darwin"* ]]; then
        # Mac OSX
        ../../../abc4pwm/bin/Mac/bayesPI -information_content_threshold=0.2 -seed=1 -split_sequences=1 \
	-initial_beta=1 -dependence=0 \
	-strand=0 -max_iteration=1000 \
	-min_L=10 -max_L=20 -p_value=0.001 \
        -max_loop=3 -normalize=0 \
	-out=${out} \
        -exp=data_in/HCT116_nutlin_S2_RNAseq_data_log2.txt \
        -seq=data_in/HCT116_nutlin_S2_RNAseq_data_log2_seq.fasta

else
    echo 'This Operating System is not supported.'
fi