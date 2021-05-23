#! /usr/bin/env python

import sys
import os
import glob

def compile_script(epsilon, number_runs):

	script = []
	j = 0

	while j < number_runs:
		scr_1 = "mpirun -np 10 ./b_m_XII_run_mpi.out 1000000 >> results/numbers__nu__5__ncg__6__e__" +str(epsilon) +"__raw_" +str(j) +".txt"
		scr_2 = "awk -f state_pull_split.awk results/numbers__nu__5__ncg__6__e__" +str(epsilon) +"__raw_" +str(j) +".txt"
		scr_3 = "awk -f awk_ave.awk 0 1 2 3 4 5 6 7 8 9 >> results/numbers__nu__5__ncg__6__e__" +str(epsilon) +"__state_ave_" +str(j) +".txt"
		scr_4 = "rm 0 1 2 3 4 5 6 7 8 9"

		script.append(scr_1)
		script.append(scr_2)
		script.append(scr_3)
		script.append(scr_4)
		j = j +1;

	scr_5 = "awk -f awk_ave.awk results/numbers__nu__5__ncg__6__e__" +str(epsilon) +"__state_ave_*.txt >> results/numbers__nu__5__ncg__6__e__" +str(epsilon) +"__state_averaged_all.txt"
	script.append(scr_5)

	return script

def dump_script(script, filename):
	output = open(filename,"w")
	
	output.write("#!/bin/bash\n")
	
	for i in script:
		print i
		output.write(i+"\n")
		
	return 1
