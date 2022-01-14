#!/bin/sh
for i in L20_1 L20_2 L20_3

do
echo $i
abc4pwm searching --pwm in/in_pwms/*$i.mlp \
	--db_path ../data/in/in_pwms/ \
	--output_directory out/search_out_$i/ \
	--db_type folder \
	--top_n 10 \
	--input_count 1 --input_prob 0 \
	--db_count 1 --db_prob 0
done
