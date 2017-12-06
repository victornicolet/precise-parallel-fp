#!/usr/bin/env bash

for i in 20 21 22 23 24 25 26 27 28 29 30
do
echo "MTS experiment for $i ------------------------"
./precise_parallel_fp 0 $i
echo "POLY experiment for $i ------------------------"
done