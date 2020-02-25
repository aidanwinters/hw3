#!/usr/bin/env python3

import numpy as np
from smith_waterman.algs import sw
from smith_waterman.io import read_score_matrix

def test_smithwaterman():
    seq1 = 'SLEAAQ'
    seq2 = 'SLEAAQ'

    sm = read_score_matrix('./BLOSUM50')
    score_matr = sw(seq1, seq2, sm, score_only=False)

    #check that the max position in SW matrix is the bottom diagonal, this would be the case for an exact match
    assert np.argmax(score_matr) == ((len(seq1) + 1) * (len(seq2) + 1) - 1)
