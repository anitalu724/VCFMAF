import nimfa
from collections import defaultdict, Counter
import urllib
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from sklearn import preprocessing
import scipy.cluster.hierarchy as sch
import pandas as pd

df = (pd.read_csv("sample_data/out_20201013/ms_input.tsv", sep="\t")).T
data = df.to_numpy()

rank_cands = range(2, 8, 1)
snmf = nimfa.Snmf(data, seed='random_vcol', max_iter=100)
summary = snmf.estimate_rank(rank_range=rank_cands, n_run=10, what='all')

rss = [summary[rank]['rss'] for rank in rank_cands]
coph = [summary[rank]['cophenetic'] for rank in rank_cands]
disp = [summary[rank]['dispersion'] for rank in rank_cands]
spar = [summary[rank]['sparseness'] for rank in rank_cands]
spar_w, spar_h = zip(*spar)
evar = [summary[rank]['evar'] for rank in rank_cands]

plt.plot(rank_cands, rss, 'o-', label='RSS', linewidth=2)
plt.plot(rank_cands, coph, 'o-', label='Cophenetic correlation', linewidth=2)
plt.plot(rank_cands, disp,'o-', label='Dispersion', linewidth=2)
plt.plot(rank_cands, spar_w, 'o-', label='Sparsity (Basis)', linewidth=2)
plt.plot(rank_cands, spar_h, 'o-', label='Sparsity (Mixture)', linewidth=2)
plt.plot(rank_cands, evar, 'o-', label='Explained variance', linewidth=2)
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=3, numpoints=1);
