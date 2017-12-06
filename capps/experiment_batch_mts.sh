#!/usr/bin/env bash

for i in 20 21 22 23 24 25 26 27 28 29 30
do
echo "Experiment for $i"
./precise_parallel_fp $i
done