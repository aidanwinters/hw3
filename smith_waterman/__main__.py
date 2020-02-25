import sys

from .algs import sw, plot_roc_line, plot_roc_final, auc, tp_fp, optimize_score_matr
from .io import read_score_matrix, read_fa, get_pos_pairs, get_neg_pairs
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


#function to perform a grid search to optimize SW params (gap open penalty, gap extension penalty)
def grid_search(open_range, extend_range, score_matrix='./BLOSUM50', file_out="grid_search.csv"):

    pos = get_pos_pairs()
    neg = get_neg_pairs()

    sm = read_score_matrix(score_matrix)

    temp = np.zeros((len(open_range),len(extend_range)))

    for op in range(len(open_range)):
        for ex in range(len(extend_range)):
            print('Testing:', open_range[op], extend_range[ex])
            pos_scores = [sw(p[0], p[1], sm, gap_open = open_range[op], gap_extend = extend_range[ex]) for p in pos]
            neg_scores = [sw(n[0], n[1], sm, gap_open = open_range[op], gap_extend = extend_range[ex]) for n in neg]
            aa = auc(*tp_fp(pos_scores, neg_scores))
            print(aa)
            temp[op,ex] = aa

    df = pd.DataFrame(temp)
    df.index = open_range
    df.columns = extend_range
    df.to_csv(file_out)

def test_score_matrices():

    score_matrices = ['BLOSUM50', 'BLOSUM62', 'MATIO', 'PAM100', 'PAM250']
    pos = get_pos_pairs()
    neg = get_neg_pairs()

    for s in score_matrices:
        print('Reading: %s' % ('./' + s))

        sm = read_score_matrix('./' + s)
        pos_scores = [sw(p[0], p[1], sm) for p in pos]
        neg_scores = [sw(n[0], n[1], sm) for n in neg]
        res = tp_fp(pos_scores, neg_scores)
        aa = auc(*res)
        plot_roc_line(*res, s + '- AUC:' + str(round(aa,2)))

    plot_roc_final()

def test_normalize_matr():

    norm = [True, False]
    pos = get_pos_pairs()
    neg = get_neg_pairs()

    for s in norm:
        print('Norm?: %s' % ( s))

        sm = read_score_matrix('./BLOSUM50')
        pos_scores = [sw(p[0], p[1], sm, normalize=s) for p in pos]
        neg_scores = [sw(n[0], n[1], sm, normalize=s) for n in neg]
        res = tp_fp(pos_scores, neg_scores)
        aa = auc(*res)
        plot_roc_line(*res, 'Normalized:' + str(s))

    plot_roc_final()

# optimize_score_matr('./BLOSUM50', orig_auc=0.84)
# optimize_score_matr('./MATIO', orig_auc=0.7)

# test_score_matrices()
# test_normalize_matr()
# grid_search(range(1,21), range(1,6), file_out='gridsearch_BLOSUM50_1to20_1to5.csv')
