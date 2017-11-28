#!/usr/bin/env bash

for i in {4194304..67108864..1048576}
do
	for j in 1 2 3 4 5 6 7 8 9
	do
	echo "n = %i"
	./precise_parallel_fp -n $i -o experiments1.csv
	done
done
