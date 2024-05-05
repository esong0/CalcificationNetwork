# Author: Jun-Seop Song (= Euijun Song)


import networkx as nx
import numpy as np
import math
from scipy.stats import norm
import warnings
warnings.filterwarnings("ignore")

INFINITY = 999999


def differential_network_proximity(d1, m1, sd1, d2, m2, sd2):
    d_diff = d1 - d2
    m_diff = m1 - m2
    sd_diff = math.sqrt(sd1**2 + sd2**2)
    z_diff = (d_diff - m_diff) / sd_diff
    p_diff = norm.cdf(z_diff)
    return d_diff, z_diff, p_diff, m_diff, sd_diff



def calc_closeness(network, nodes_S, nodes_T, num_random, rnd_nodes_S_set, rnd_nodes_T_set):
    dist = calc_average_shortest_distance(network, nodes_S, nodes_T)
    
    dist_distribution = []
    for i in range(num_random):
        rnd_nodes_S = rnd_nodes_S_set[i]
        rnd_nodes_T = rnd_nodes_T_set[i]
        try:
            dist_distribution.append(calc_average_shortest_distance(network, rnd_nodes_S, rnd_nodes_T))
        except:
            continue
    m, sd = np.nanmean(dist_distribution), np.nanstd(dist_distribution)

    z = 0.0 if sd == 0 else (dist - m) / sd
    p = norm.cdf(z)

    return dist, z, p, m, sd



def calc_average_shortest_distance(network, nodes_S, nodes_T):
    dist = []
    for node_T in nodes_T:
        dist_ = []
        for node_S in nodes_S:
            try:
                d = nx.shortest_path_length(network, node_S, node_T)
            except:
                continue # d = INFINITY
            dist_.append(d)
        dist.append(np.nanmean(dist_))
    return np.nanmean(dist)
