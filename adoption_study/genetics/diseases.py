
#
# libraries
#
import os
import json
import pprint as pp
import networkx as nx
import math

import matplotlib.pyplot as plt

#
# user settings
#
compute_ego_network = False
output_directory = 'output'
diseases_list_file = 'data/diseases_associated_with_violence.txt'
skip_words_file = 'data/words_and_phrases_to_skip.txt'


#
# load files
#
f = open(diseases_list_file)
diseases = list(set([x.strip().lower() for x in f.readlines() if x.strip() != '' and x[0] != '#']))
f.close()
f = open(skip_words_file)
skip_words = list(set([x.strip().lower() for x in f.readlines() if x.strip() != '' and x[0] != '#']))
f.close()

#
# extract from database
#
os.system('cat sql.sql | sqlite3 data/disgenet_2017.db > ' + output_directory + '/db_output.txt')

#
# remove duplicates in DB output
#
f = open('output/db_output.txt')
de_dup_dict = {}
for line in f:
    line = line.strip()
    de_dup_dict[line] = None
f.close()

#
# process lines
#
disease_2_gene = {}
gene_2_disease = {}
gene_id_2_symbol = {}
for line in de_dup_dict.keys():
    line = line.split('|')
    if line[0] == 'geneId':
        continue

    gene_id = int(line[0])
    gene_symbol = line[1]
    cui = line[2]
    name = line[3].lower()

    gene_id_2_symbol[gene_id] = gene_symbol

    for d in diseases:
        if name.find(d) >= 0:
            skip_word = False
            for s in skip_words:
                if name.find(s) >= 0:
                    skip_word = True
            if not skip_word:
                
                if not name in disease_2_gene:
                    disease_2_gene[name] = {}
                disease_2_gene[name][gene_id] = None

                if not name in gene_2_disease:
                    gene_2_disease[gene_id] = {}
                gene_2_disease[gene_id][name] = None

#
# make graph
#
graph_as_dict = {}
nodes = {}
for i in gene_2_disease.keys():
    for disease in gene_2_disease[i].keys():
        for j in disease_2_gene[disease]:
            if i != j:
                ii = min(i, j)
                jj = max(i, j)

                nodes[ii] = {'symbol' : gene_id_2_symbol[ii]}
                nodes[jj] = {'symbol' : gene_id_2_symbol[jj]}
                
                key = str(ii) + '_' + str(jj)
                if not key in graph_as_dict:
                    graph_as_dict[key] = 0
                graph_as_dict[key] += 1

#
# build graph
#
G = nx.Graph()
for node in nodes.keys():
    G.add_node(node)
for key in graph_as_dict.keys():
    ii = int(key.split('_')[0])
    jj = int(key.split('_')[1])
    weight = graph_as_dict[key]
    G.add_edge(ii, jj, weight=weight)

#
# ego network
#
if compute_ego_network:
    for n in nodes.keys():
        eG = nx.ego_graph(G, n)
        eG_size = eG.size()
        nodes[n]['ego_network_size'] = eG_size
        
    with open(output_directory + '/nodes.json', 'w') as f:
        json.dump(nodes, f, indent=4)




#
# node clique number
#
#for node in nx.node_clique_number(G):
#    print(node)

#
# plot
#
n = len(nodes.keys())
k = 1. / math.sqrt(n)
pos = nx.spring_layout(G)
nx.draw_networkx(G, pos=pos, with_label=False, node_size=2, node_color='blue', edge_color='grey')

nx.draw_spring(G)
plt.savefig(output_directory + '/plot.png')
plt.close()




# nx.cliques_containing_node(G, nodes=None, cliques=None)
# 





            

