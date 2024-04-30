# Author: Jun-Seop Song (= Euijun Song)


import networkx as nx
import scipy
from scipy.stats import norm
import statistics as stat
import warnings
warnings.filterwarnings("ignore")


def calc_stat_rwr_score(network, nodes_S, nodes_T, num_random,
                        rnd_nodes_S_set, rnd_nodes_T_set):
    '''
    Calculate z-score of RWR score

    nodes_S: disease module genes
    nodes_T: perturbation genes (miRNA targets)
    num_random: number of randomization (e.g., 1000)
    rnd_nodes_S_set: set of randomized nodes for nodes_S
    rnd_nodes_T_set: set of randomized nodes for nodes_T
    '''
    W = normalized_adjacency_matrix(network)
    score = calc_rwr_score(network, W, nodes_S, nodes_T)
    
    score_distribution = []
    for i in range(num_random):
        rnd_nodes_S = rnd_nodes_S_set[i]
        rnd_nodes_T = rnd_nodes_T_set[i]
        s = calc_rwr_score(network, W, rnd_nodes_S, rnd_nodes_T)
        score_distribution.append(s)
    m, sd = stat.mean(score_distribution), stat.stdev(score_distribution)
    
    z = 0.0 if sd == 0 else (score - m) / sd
    p = norm.sf(z)
    
    return score, z, p, m, sd # RWR score, z-score, p-value, mean, sd



def calc_rwr_score(G, Wp, nodes_S, nodes_T):
    '''
    Calculate RWR score for disease module

    nodes_S: disease module genes
    nodes_T: perturbation genes (miRNA targets)
    Wp: normalized directed adjacency matrix of G
    '''
    rwr_result = rwr_process(G, Wp, nodes_T, alpha=0.5, max_iter=100, tol=1.0e-6)
    rwr_score = scipy.mean([rwr_result[v] for v in nodes_S])

    return rwr_score



def rwr_process(G, Wp, initial_seed, alpha=0.5, max_iter=100, tol=1.0e-6):
    '''
    Random walk with restart (RWR) process
    x_update = alpha*W*x + (1-alpha)*x0

    Wp: normalized directed adjacency matrix
    initial_seed: perturbation genes (miRNA targets)
    alpha: restart probability is (1-alpha)
    '''
    N = len(G)
    nodelist = list(G)
    
    if Wp is None:
        W = normalized_adjacency_matrix(G)
    else:
        W = Wp
    
    # initial vector
    x0 = scipy.array([int(v in initial_seed) for v in nodelist], dtype=float)
    x0 = x0 / x0.sum()
    x = x0
    
    # iteration
    for _ in range(max_iter):
        x_update = x
        x = alpha * W.dot(x) + (1 - alpha) * x0
        err = scipy.absolute(x - x_update).sum()
        if err < N * tol:
            break
    
    return dict(zip(nodelist, map(float, x)))



def normalized_adjacency_matrix(G):
    '''
    Calculate normalized directed adjacency matrix of G (directed graph)
    '''
    nodelist = list(G)
    A = nx.to_scipy_sparse_matrix(G, nodelist=nodelist, dtype=float)
    M1 = scipy.array(A.sum(axis=1)).flatten()
    M2 = scipy.array(A.sum(axis=0)).flatten()
    M1[M1 != 0] = 1.0 / scipy.sqrt(M1[M1 != 0])
    M2[M2 != 0] = 1.0 / scipy.sqrt(M2[M2 != 0])
    D1 = scipy.sparse.spdiags(M1.T, 0, *A.shape, format='csr')
    D2 = scipy.sparse.spdiags(M2, 0, *A.shape, format='csr')
    W = D1 * A * D2

    return W
