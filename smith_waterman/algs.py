import numpy as np
from sklearn import metrics
import random
from .io import read_score_matrix, get_pos_pairs, get_neg_pairs

def sw(seqA, seqB, score_matrix, gap_open = 6, gap_extend = 5, score_only=True, normalize=False):

    H = np.zeros((len(seqA) + 1, len(seqB) + 1), np.int)
    H_gaps = np.zeros(H.shape, np.int) # -1 = verticla gap, 1 = horizontal gap

    for i in range(1, len(seqA) + 1):
        for j in range(1, len(seqB) + 1):

            #get the score for a match, which is the score of diagonal plus the score for matching these two chars
            match = H[i-1, j-1] + score_matrix.loc[seqA[i-1], seqB[j-1]]
            vertical = H[i - 1, j] - (gap_extend if H_gaps[i-1, j] == -1 else gap_open)
            horizontal = H[i, j-1] - (gap_extend if H_gaps[i, j-1] ==  1 else gap_open)

            #using a trick from Hasan (instead of a couple If statements), simultaneously update the Score matrix and gap matrix depending on which is max
            H[i,j], H_gaps[i,j] = max((match, 0), (vertical, -1), (horizontal, 1), (0,0))

    if(normalize):
        H = H/min(len(seqA), len(seqB))

    if(score_only):
        return max(H.flatten())
    else:
        return H

#calculate true positive rate and false positve rate given
#a list of positive scores and negative scores
def tp_fp(pos_vals, neg_vals):

    #currently threshold wherever there is a score
    thresholds = (pos_vals + neg_vals)
    thresholds.sort()
    tpr, fpr = [], []

    for thresh in thresholds:
        total_pos = 0
        total_neg = 0

        for pos in pos_vals:
            if pos >= thresh:
                total_pos += 1

        for neg in neg_vals:
            if neg >= thresh:
                total_neg += 1

        tpr.append(total_pos/len(pos_vals))
        fpr.append(total_neg/len(neg_vals))

    return tpr, fpr

def auc(tpr, fpr):
    return metrics.auc(fpr, tpr)


#mutate the matrix based on rules in the params
#mutate_num is a fixed number of elements to randomly select for mutation
#mutate_change is the amount the score position can change
def mutate_score_matrix(sm, mutate_num=50, mutate_change=3):

    sm_copy = sm.copy()

    x = np.random.randint(0, high=sm_copy.shape[0], size=mutate_num)
    y = np.random.randint(0, high=sm_copy.shape[0], size=mutate_num)
    sign = [1 if random.random() < 0.5 else -1 for x in range(mutate_num)]

    for i in range(mutate_num):
        sm_copy.iloc[x[i], y[i]] += mutate_change*sign[i]
        sm_copy.iloc[y[i], x[i]] += mutate_change*sign[i]

    return sm_copy

#Optimization function
#This function is essentially a genetic aaproach.
#The score matrix is mutated num_mutants number of ways
#Then positives and negatives are evaluated on mutants
#
#Allows for optional original auc argument in order to skip recomputing
def optimize_score_matr(sm_path, auc_desired=0.05, max_iter=15, num_mutants=3, orig_auc=None):

    orig_sm = read_score_matrix(sm_path)

    pos = get_pos_pairs()
    neg = get_neg_pairs()

    if orig_auc is None:
        pos_scores = [sw(p[0], p[1], orig_sm) for p in pos]
        neg_scores = [sw(n[0], n[1], orig_sm) for n in neg]
        orig_auc = auc(*tp_fp(pos_scores, neg_scores))
        print('Done with original Score Matrix')
    else:
        print('Original AUC provied: ', orig_auc)

    curr_sm = orig_sm
    curr_auc = orig_auc
    print(curr_sm)

    #only continue for certain num iterations
    for m in range(max_iter):
        print('Working on Iter:', m)

        #generate mutant scoring matrices
        mutants = [mutate_score_matrix(curr_sm) for s in range(num_mutants)]
        for mut in mutants:
            print('\tWorking on Mutant')
            print(mut)
            pos_scores = [sw(p[0], p[1], mut) for p in pos]
            neg_scores = [sw(n[0], n[1], mut) for n in neg]
            aa = auc(*tp_fp(pos_scores, neg_scores))
            print('\tAUC:', aa)
            #check auc
            if aa > curr_auc:
                if aa > (orig_auc + auc_desired):
                    print('Reached desired AUC:', aa)
                    mut.write_csv(sm_path + '_optimized')
                print('Old Auc: %s - Mutant AUC: %s' % (str(round(curr_auc)), str(round(aa))))
                curr_sm = mut
                curr_auc = aa

    print('final auc:', curr_auc)
    curr_sm.to_csv(sm_path + '_optimized2')
    return curr_sm
