#! /usr/bin/env python

import sys
import os
import glob

def compile_script(epsilon, number_runs):

	script = []
	j = 0

	while j < number_runs:
		scr_1 = "awk -f state_pull_split.awk results/numbers__nu__5__ncg__6__e__" +str(epsilon) +"__raw_" +str(j) +".txt"
		scr_2 = "awk -f awk_ave_length_sq.awk 0 1 2 3 4 5 6 7 8 9 >> results/numbers__nu__5__ncg__6__e__" +str(epsilon) +"__sq_length_ave_" +str(j) +".txt"
		scr_3 = "rm 0 1 2 3 4 5 6 7 8 9"

		script.append(scr_1)
		script.append(scr_2)
		script.append(scr_3)
		j = j +1;

	scr_4 = "awk -f awk_sum.awk results/numbers__nu__5__ncg__6__e__" +str(epsilon) +"__sq_length_ave_*.txt >> results/numbers__nu__5__ncg__6__e__" +str(epsilon) +"__sq_length_ave_all.txt"
	script.append(scr_4)

	return script

def dump_script(script, filename):
	output = open(filename,"w")
	
	output.write("#!/bin/bash\n")
	
	for i in script:
		print i
		output.write(i+"\n")
		
	return 1
