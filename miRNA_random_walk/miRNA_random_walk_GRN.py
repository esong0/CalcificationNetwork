# Author: Jun-Seop Song (= Euijun Song)


import networkx as nx
import numpy as np
from random_walk_restart import *
from randomization import *
import warnings
warnings.filterwarnings("ignore")


disease_set = ['Coronary artery disease',
               'Myocardial Infarction',
               'Aneurysm',
               'Cardiomyopathies',
               'Atherosclerosis',
               'Hypertension',
               'Atrial Fibrillation']
disease = disease_set[0]
print(disease)


# GRN import
G = nx.DiGraph()

for line in open('../GRN_human/GRN_human_202006.txt', 'r'):
    line_data = line.strip().split('\t')
    gene1, gene2 = line_data[0], line_data[1]
    G.add_edge(gene1, gene2)

print("[Interactome]\n%s\n" % nx.info(G))



# Seed import
seed = []
for line in open('disease_miRNA_gene_dataset/disease_gene_lcc.txt', 'r'):
    if line[0] == '#':
        continue
    line_data = line.strip().split('\t')
    if line_data[0] == disease:
        seed = line_data[2].split(',')
print("%s original seed genes: %d\n" % (disease, len(seed)))
##seed = list(max(nx.weakly_connected_components(G.subgraph(seed)), key=len))
seed = list(G.subgraph(seed))
print("%s seed genes: %d\n" % (disease, len(seed)))


# miRNA-target network
DT = nx.Graph()
miRNA_set = set()

for line in open('../miRNA_database/miRTarBase_human_v7.0_strongevidence.txt', 'r'):
    line_data = line.strip().split('\t')
    miRNA_id = line_data[0]
    gene = line_data[1]
    if gene not in G.nodes:
        continue
    DT.add_edge(gene, miRNA_id)
    miRNA_set.add(miRNA_id)

print("[miRNA-target network]\n%s\n" % nx.info(DT))
print("miRNA: %d\n" % len(miRNA_set))


### RWR analysis
lengths = None
np.random.seed(0) # for reproducibility


# Generate random samples
print('Generating random samples...')
num_random = 500
bins = get_degree_binning(G, 100, lengths)
rnd_nodes_S_set = get_random_nodes(seed, G, bins = bins,
                                   n_random = num_random, min_bin_size = 100, seed = 666)


# Calculate RWR score
print('Calculating RWR...\n')
f = open('random_walk_restart_GRN_' + disease + '.txt', 'w')
f.write('#miRNA\tRWR_score\tZ-score\tP\tm\tsd\tTarget_size\tTargets\n')

for idx, miRNA in enumerate(sorted(miRNA_set)):
    print(idx)
    try:
        target_size = len(list(DT.neighbors(miRNA)))
        rnd_nodes_T_set = get_random_nodes(list(DT.neighbors(miRNA)), G, bins = bins,
                               n_random = num_random, min_bin_size = 100, seed = 666)
        s, z, p, m, sd = calc_stat_rwr_score(G, seed, list(DT.neighbors(miRNA)),
                                             num_random, rnd_nodes_S_set, rnd_nodes_T_set)
        
        target_print = ",".join(sorted(list(DT.neighbors(miRNA))))
        
        print('%s\t%e\t%e\t%e\t%e\t%e\t%d\t%s'
              % (miRNA, s, z, p, m, sd, target_size, target_print))
        f.write('%s\t%e\t%e\t%e\t%e\t%e\t%d\t%s\n'
                % (miRNA, s, z, p, m, sd, target_size, target_print))
        f.flush()
    except:
        continue

f.close()
