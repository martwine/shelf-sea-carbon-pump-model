#!/bin/bash

for run in experiments/$1/*.R; do
	Rscript experiment_run.R $run $1
done
