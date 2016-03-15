#!/bin/sh
cd runs
#../bin/rpower ../data/russell_1000_cov.txt -s 0 -q 1 -w 1 -r 1 # only 1 unperturbed job
../bin/rpower ../data/russell_1000_cov.txt -s 10 -q 5 -w 2 -r 2 -t 1e-6
