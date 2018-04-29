#!/bin/bash

# Run one of the two loops below, but not both

############
# Using qsub to run simulation code
# Arguments:
#	m, n, iteration_start, iteration_end
############
for iter in $(eval echo {$3..$4}) do
	qsub -cwd simu_finalized.R $1 $2 $iter 
done

############
# Run simulation code sequentially, may want to parallelize simu_finalized.R instead for efficiency
# Arguments:
#	m, n, iteration_start, iteration_end
############
for iter in $(eval echo {$3..$4}) do
	simu_finalized.R $1 $2 $iter
done