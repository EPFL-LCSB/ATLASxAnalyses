# LCSB EPFL 2019, Jasmin Hafner
# This script plots the network analys

import networkx as nx
import os
import timeit
import re
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
import sys

################## PARAMETERS ###################

CAR_cutoff = 0.34 # Best predictor for KEGG RPAIR of type "main"

################## FILE NAME ###################

rpair_dir_path = "../Data/"
output_folder = "../../../ATLASxAnalyses_output"
project_folder = "/ComponentDistribution/"
if not os.path.isdir(output_folder): 
	os.mkdir(output_folder)
if not os.path.isdir(output_folder+project_folder): 
	os.mkdir(output_folder+project_folder)
	
################# OPTIONS ###################
assert 1 <= len(sys.argv) <= 2, "ERROR: only 0 or 1 arguments allowed."
if len(sys.argv) > 1 and sys.argv[1] == 'data_sources':
	scopes = ['KEGG','MetaCyc','SEED', 'BiGG', 'BKMS', 'Brenda', 'MetaNetX', 'Reactome', 'Rhea']
	colors_pal = ['#E86685', '#C58C22', '#97A424', '#57B363', '#58AFA4', '#58A8D2', '#A087F8', '#E550DE', '#B2B2B2'] # for data sources
	output_file_path = output_folder + project_folder + "component_distribution_ATLASx_datasources_CAR"+ str(CAR_cutoff) +".png"
else:	
	scopes = ['bioDB', 'bioATLAS', 'chemATLAS']
	colors_pal = ['#48a921','#105a37','#339098'] # for bioATLAS, chemATLAS
	output_file_path = output_folder + project_folder + "component_distribution_ATLASx_scopes_CAR"+ str(CAR_cutoff) +".png"
	
################## FUNCTIONS GRAPH ########################
def getComponentStats(indir, outfile):
	data_dict = {}
	for index, db in enumerate(scopes): 
		G = nx.Graph()
		print("Load network:", db)
		G = nx.read_gpickle(indir + 'rpairs_' + db + '_hp.gpickle')
		data_dict[db] = calculateStatistics(G, db)
	with open(outfile, 'wb') as handle:
		pickle.dump(data_dict, handle)

	
def plotData(infile, outdir):
	with open(infile, 'rb') as handle:
		data_dict = pickle.load(handle)

	plt.figure(figsize=(7,5), dpi=300)
	ax = plt.axes(yscale='log', xscale = 'log')
	linestyles = [':', '--', '-', ':', '--', '-', ':', '--', '-', ':', '--', '-']
	for index, db in enumerate(scopes): 
		color = colors_pal[index]
		y = range(1,len(data_dict[db])+1)
		print(data_dict[db][0:100])
		plt.plot(y, data_dict[db], linewidth=1, color = color, linestyle = linestyles[index], marker='o', markersize=2, label = db)
	plt.ylabel('Component size (number of compounds)')	
	plt.xlabel('Components, ordered by size')	
	plt.legend()
	plt.savefig(output_file_path, bbox_inches = 'tight', dpi=300)


def calculateStatistics( G, db_name):
	#create new unweighted graph with CAR >= threshold
	simple_G = nx.Graph()
	for (u, v, c) in G.edges.data('car'):
		if c > CAR_cutoff:
			simple_G.add_edge(u, v)
	component_sizes = []
	[component_sizes.append(len(c)) for c in sorted(nx.connected_components(simple_G), key=len, reverse=True)]
	print(db_name, 'component_sizes: ', len(component_sizes))
	print(db_name, 'connected_components: ', nx.number_connected_components(simple_G))
	plot_degree_dist(simple_G, output_folder + project_folder +'degree_dist_'+db_name+'_'+str(CAR_cutoff)+'.png', db_name)
	return component_sizes
	

def getDate(outfile):
	today = date.today()
	d = today.strftime("%B %d, %Y")
	outfile.write('#date_download\t%s\n'%d)
		
def group(s):
	if not s.isdigit():
		return s
	groups = []
	while s and s[-1].isdigit():
		groups.append(s[-3:])
		s = s[:-3]
	return s + "'".join(reversed(groups))

def plot_degree_dist(G, outpath, scope):
    degrees = [G.degree(n) for n in G.nodes()]
    x = []
    y = []
    probabilities = []
    for d in sorted(list(set(degrees))):
    	x.append(d) 
    	y.append(degrees.count(d)/len(degrees))
    	# print(d, degrees.count(d)/len(degrees))
    # plt.figure(figsize=(7,5), dpi=300)
    plt.yscale('log')
    plt.xscale('log')
    plt.scatter(x,y, s=10, label=scope)
    plt.ylabel('P(k)')
    plt.xlabel('k')
    plt.legend()
    plt.savefig(outpath, bbox_inches = 'tight', dpi=300)


    
################## MAIN ########################

def main(rpair_dir_path, output_file_path): 
	start_time = timeit.default_timer()
	getComponentStats(rpair_dir_path, output_folder + project_folder +'component_data.pickle')
	plotData(output_folder + project_folder +'component_data.pickle', output_folder + project_folder )
	print("Runtime:", str(timeit.default_timer() - start_time) + ' s')

main(rpair_dir_path, output_file_path)






