# LCSB EPFL 2019, Jasmin Hafner

import networkx as nx
import time
import pandas as pd
import os
import multiprocessing # pathway search is parallelized

# Files
input_precursor_target_csv = 'data/MetacycExtractedLinearPathways.csv'
db_scopes = ['MetaCyc', 'MetaCyc_hp', 'chemATLAS', 'chemATLAS_hp']

distance_type = 'dist' # simple distance calculate as 1/CAR
CAR_cutoff = 0.34 # standard CAR cutoff as described in NICEpath

max_num_pathways = 10000 # increasing this threshold increases the time for pathway search

#################### Defaults #####################
header = ['id', 'length', 'intermediates', 'pairs_id', 'pairs_weight', 'if_native']
NetworkFilesFolder = '../NetworkAnalysis/Data/'

##### Graph truncation ####
def getSimpleG(G):
    # Removing the edges that are below the defined CAR
    simple_G = nx.Graph()
    for (u, v, c) in G.edges.data('car'):
        if c >= CAR_cutoff:
            simple_G.add_edge(u, v)
            simple_G[u][v]['car'] = G[u][v]['car']
            simple_G[u][v]['dist'] = G[u][v]['dist']
    return simple_G

########################################################
#      Check the right pathway as edges of Network     #
########################################################
def addEdgesCoverage(goldenDF):
    # Check for each edge is it is included into the network
    for NT in db_scopes:
        graph_name = 'rpairs_'+NT
        G = nx.read_gpickle(NetworkFilesFolder + graph_name + ".gpickle")
        Gs = getSimpleG(G)
        # getting simpler version of edges for faster comparison
        G_edges = []
        for edge in Gs.edges:
            e = list(edge)
            e.sort()
            G_edges.append(str(e[0])+'_'+str(e[1]))
        G_edges = set(G_edges)
        print('edges for', NT, 'calculated')
        goldenDF[NT+'_coveredByNetwork']= goldenDF.apply(check_all_edges_in_network, args = (G_edges), axis=1)
        # Check if any pathway exists between the two compounds in the network
        goldenDF[NT+'_hasPath']=goldenDF.apply(check_if_has_path, args=(G,), axis=1)

    goldenDF.to_csv('output/networkEdgesCoverage.csv', index = False)

def check_all_edges_in_network(row, G_edges):
    # loading the native pathways
    native_pathway = load_pathway(row)
    # checking if all the edges of the native pathway are in the network
    edges_native = []
    for edge in native_pathway:
        e = list(edge)
        e.sort()
        edges_native.append(str(e[0])+'_'+str(e[1]))
    edges_covered = G_edges.intersection(set(edges_native))
    if len(edges_covered)==len(native_pathway):
        return 1
    else:
        return 0

def load_pathway(row):
    comps = [int(i) for i in row['Pathway_original_in_refV'].split('|')]
    edges = [(comps[i-1], comps[i]) for i in range(1, len(comps))]
    return edges

def check_if_has_path(row, G):
    if row['Precursor LCSB ID'] in G.nodes and row['Target LCSB ID'] in G.nodes:
        return nx.has_path(G, row['Precursor LCSB ID'], row['Target LCSB ID'])
    return False

#################### PATHWAY SEARCH FUNCTIONS ######################

def shortest_paths(Graph, source, target, weight):
    try:
        return nx.shortest_simple_paths(Graph, source, target, weight=weight)
    except Exception as e:
        print(source, target, e)
        return list()

def find_pathway_rank(param_list):
    #param_list: 0:prec; 1:target; 2:pwid; 3:ref_seq_of_intermediates
    entry = 1

    reference_path = [int(i) for i in param_list[3].split('|')]

    for path in shortest_paths(G, param_list[0], param_list[1], distance_type):
        if entry > max_num_pathways:
            return 0
        pathway_match = [i for i, j in zip(reference_path, path) if i == j]
        if len(pathway_match) == len(reference_path):
            return entry
        entry += 1

##################### MAIN ######################

if not os.path.exists('output/networkEdgesCoverage.csv'):
    metaCycPwDataset = pd.read_csv(input_precursor_target_csv)
    addEdgesCoverage(metaCycPwDataset)

metaCycPwDataset = pd.read_csv('output/networkEdgesCoverage.csv')

# Pathway search
local_time = time.ctime(time.time())
print("Time execution starts:", local_time)

for db_scope in db_scopes:
    column_network = db_scope+"_coveredByNetwork"
    column_path = db_scope+"_hasPath"
    df_scope = metaCycPwDataset.query('{}==1&{}==1'.format(column_network, column_path))
    print(db_scope)
    G = nx.read_gpickle(NetworkFilesFolder + "rpairs_" + db_scope + ".gpickle")
    G = getSimpleG(G)

    print("Find paths...", db_scope)

    param_list = list(zip(df_scope['Precursor LCSB ID'].to_list(), df_scope['Target LCSB ID'].to_list(), df_scope['MetaCyc Pathway ID'].to_list(), df_scope['Pathway_original_in_refV'].to_list()))

    # pathway search is parallelized
    pool = multiprocessing.Pool()
    pw_ranks = pool.map(find_pathway_rank, param_list)
    # pool has to be closed to avoid OSError: [Errno 24] Too many open files
    pool.close()

    df_scope['pw_rank'] = pw_ranks

    df_scope.to_csv('output/pw_ranking_'+db_scope+'.csv', index = False)

local_time = time.ctime(time.time())
print("Time execution ends:", local_time)