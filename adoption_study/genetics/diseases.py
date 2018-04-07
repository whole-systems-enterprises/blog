
#
# libraries
#
import os
import json
import pprint as pp
import networkx as nx
import math
import numpy as np

import matplotlib.pyplot as plt

#
# user settings
#
compute_ego_network = False
make_plot = True
make_table = False
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
os.system('cat diseases_sql.sql | sqlite3 data/disgenet_2017.db > ' + output_directory + '/db_output.txt')

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


with open(output_directory + '/nodes.json') as f:
    nodes = json.load(f)

#
# make table
#
if make_table:
    content = {}
    for gene_id in nodes.keys():
        size = nodes[gene_id]['ego_network_size']
        symbol = nodes[gene_id]['symbol']
        if not size in content:
            content[size] = []
        content[size].append({
                'gene_id' : gene_id,
                'gene_symbol' : symbol,
                })

    score_list = sorted(content.keys())
    score_list.reverse()

    html = ''
    html += '<table class="emily-table">' + '\n'
    html += '<tr class="emily-tr"><th class="emily-th">' + '</th><th class="emily-th">'.join(['Network Size', 'NCBI Gene Symbol', 'NCBI Gene ID']) + '</th></tr>' + '\n'

    for size in score_list:
        for item in content[size]:
            gene_id = item['gene_id']
            gene_symbol = item['gene_symbol']
            html += '<tr class="emily-tr"><td class="emily-td">' + '</td><td class="emily-td">'.join([str(size), gene_symbol, gene_id]) + '</td></tr>' + '\n'
    html += '</table>\n'

    f = open(output_directory + '/table.html', 'w')
    f.write(html)
    f.close()



#
# plot
#
if make_plot:

    size_list = []
    to_remove = []
    for gene_id in nodes.keys():
        size_list.append(nodes[gene_id]['ego_network_size'])
    p = np.percentile(size_list, [20])[0]
    for gene_id in nodes.keys():
        if nodes[gene_id]['ego_network_size'] > p:
            to_remove.append(int(gene_id))


    GG = G.copy()
    GG.remove_nodes_from(to_remove)


    n = len(nodes.keys())
    k = 1. / math.sqrt(n)

    plt.figure(figsize=[15, 12])

    pos = nx.spring_layout(GG, k=k)
    nx.draw_networkx(GG, pos=pos, with_labels=False, node_size=2, node_color='red', edge_color='darkgray')


    plt.savefig(output_directory + '/plot.png')
    plt.close()



            

