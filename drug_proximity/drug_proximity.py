# Author: Jun-Seop Song (= Euijun Song)


import networkx as nx
import statistics as stat
from scipy.stats import norm
import warnings
warnings.filterwarnings("ignore")

INFINITY = 999999

def calc_proximity(network, nodes_S, nodes_T, num_random,
                   rnd_nodes_S_set, rnd_nodes_T_set):
    if len(nodes_S) == 0 or len(nodes_T) == 0:
        return INFINITY, INFINITY

    dist = calc_closest_distance(network, nodes_S, nodes_T)

    dist_distribution = []
    for i in range(num_random):
        rnd_nodes_S = rnd_nodes_S_set[i]
        rnd_nodes_T = rnd_nodes_T_set[i]
        d = calc_closest_distance(network, rnd_nodes_S, rnd_nodes_T)
        if d != INFINITY:
            dist_distribution.append(d)
    m, sd = stat.mean(dist_distribution), stat.stdev(dist_distribution)

    z = 0.0 if sd == 0 else (dist - m) / sd
    p = norm.cdf(z)

    return dist, z, p # proximity, z-score, p-value


def calc_closest_distance(network, nodes_S, nodes_T):
    min_dist = []
    for node_T in nodes_T:
        min_d = INFINITY
        for node_S in nodes_S:
            try:
                d = nx.shortest_path_length(network, node_S, node_T)
            except:
                continue  # alternative: d = INFINITY
            if d < min_d:
                min_d = d
                if d == 0:
                    break
        min_dist.append(min_d)

    return stat.mean(min_dist)
