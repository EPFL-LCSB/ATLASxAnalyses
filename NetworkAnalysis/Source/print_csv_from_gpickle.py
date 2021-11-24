# LCSB EPFL 2020, Jasmin Hafner
# This script translates a gpickle graph file into a csv file suitable for Gephi import

import networkx as nx
import pickle
import timeit
import sys
import os

start_time = timeit.default_timer()
################## PARAMETERS ###################

CAR_cutoff = 0.34 # Best predictor for KEGG RPAIR of type "main"

################# OPTIONS ###################
assert 1 <= len(sys.argv) <= 2, "ERROR: only 0 or 1 arguments allowed."
if len(sys.argv) > 1 and sys.argv[1] == 'data_sources':
	scopes = ['KEGG','MetaCyc','SEED', 'BiGG', 'BKMS', 'Brenda', 'MetaNetX', 'Reactome', 'Rhea', 'bioDB', 'bioATLAS', 'chemATLAS']
else:	
	scopes = ['bioDB', 'bioATLAS', 'chemATLAS']


################## FILE NAME ###################

input_dir_path = "../Data/"
output_folder = "../../../ATLASxAnalyses_output"
project_folder = "/CSV_files/"
if not os.path.isdir(output_folder): 
	os.mkdir(output_folder)
if not os.path.isdir(output_folder+project_folder): 
	os.mkdir(output_folder+project_folder)
output_dir_path = output_folder + project_folder

################## MAInn ###################
def main(input_dir, output_dir, scope):
	rpair_file_path = input_dir + "rpairs_" +  scope + "_hp.gpickle"
	G = nx.Graph()
	print("Load network:", scope)
	G = nx.read_gpickle(rpair_file_path)
	outfile = open(output_dir + 'EdgeList_' + scope + '_' + str(CAR_cutoff) + '.csv', 'w')
	outfile.write('Source,Target,CAR\n')
	#if CAR >= threshold -> write to table
	for (u, v, c) in G.edges.data('car'):
		if c > CAR_cutoff:
			car = round(c, 3)
			outfile.write('%d,%d,%s\n' %(u,v,str(car)))
	outfile.close()




for scope in scopes:
	main(input_dir_path, output_dir_path, scope)

print("Runtime:", str(timeit.default_timer() - start_time) + ' s')

