import pandas as pd
from Bio import SeqIO


#file: a string file path
#
def read_score_matrix(file, comment_char='#'):

    #read in file with at least 1 space as delimiter
    df = pd.read_csv(file, sep='\s+', comment=comment_char, engine='python')
    df.index = df.columns

    return df

#returns a SeqIO FASTA object with two fields:
#id: the id of the fasta seq

def read_fa(file):
    return SeqIO.read(file, "fasta")


def get_pos_pairs(pos_pair_file='./Pospairs.txt'):

    fas = []

    with open(pos_pair_file) as file:
        for line in file:
            spl = line.strip().split(' ')
            fas.append([str(read_fa(spl[0]).seq).upper(), str(read_fa(spl[1]).seq).upper()])
    return fas

def get_neg_pairs(neg_pair_file='./Negpairs.txt'):

    fas = []

    with open(neg_pair_file) as file:
        for line in file:
            spl = line.strip().split(' ')
            fas.append([str(read_fa(spl[0]).seq).upper(), str(read_fa(spl[1]).seq).upper()])
    return fas

#add a line to the current ROC plot
def plot_roc_line(tpr, fpr, lbl):

    plt.plot(fpr, tpr, '-', label=lbl)
    plt.xlabel('FPR')
    plt.ylabel('TPR')
    plt.xlim(0,1)
    plt.ylim(0,1)

#plot the final roc plot (assuming all lines have already been plotted)
def plot_roc_final():
    plt.plot([0,1], [0,1], '--', c='red')
    plt.legend()
    plt.show()
