import sys
sys.path.append("/Users/joazofeifa/Lab/motif_displacement_analysis/motif_displacement_analysis_package/")
import matplotlib.pyplot as plt
import numpy as np,math
from scipy.interpolate import spline
from scipy.special import erf
import scipy.stats as ss
import time
from matplotlib import cm
import matplotlib as mpl
from statsmodels.stats import multitest
import pandas as pd,os
from scipy import stats
import normal_laplace as nl
import itertools, display,experiment
import collect_stats 
class mds_frame:
    def __init__(self):
        self.EXPS           = {} #inserted experiment IDS
    def load_MD_score_file(self, *args, **kwargs):
        if len(args) < 2:
            print "please specify at least one motif displacement file and one experiment ID"
            return -1
        ID  = args[-1]
        for i in range(len(args)-1):
            exp  = experiment.dataset( args[i] )

            if ID not in self.EXPS:
                self.EXPS[ID]   = list()
            self.EXPS[ID].append(exp)
    def re_index(self,a,b):
        X,Y             = list(),list()
        '''
            look at tss or look at eRNAs?
        '''
        motifs  = list(set(self.EXPS[a][0].df.ID))
        for i in range(len(self.EXPS[a])):
            df  = self.EXPS[a][i].df.reindex(motifs)
            X.append(np.array(df.as_matrix()[:,1:],dtype=float))
        for i in range(len(self.EXPS[b])):
            df  = self.EXPS[b][i].df.reindex(motifs)
            Y.append(np.array(df.as_matrix()[:,1:],dtype=float))
        '''
            reshuffle the rows
        '''
        return motifs,X,Y
    def differential(self, A,B,h=150):

        if A not in self.EXPS or B not in self.EXPS:
            print A, "or", B,"has not been added as an experiment to this mds_frame"
            return -1
        n,m             = len(self.EXPS[A]),len(self.EXPS[B])
        if n <= m:
            comparisons     =  [(i,j) for i in range(n) for j in range(i,m)  ]
        else:
            comparisons     =  [(i,j) for j in range(m) for i in range(j,n)  ]
        '''
            re-index
        '''
        motifs, X,Y         = self.re_index(A,B)
        SB                  = np.zeros((len(motifs), 7*len(comparisons)))
        for k,(i,j) in enumerate(comparisons):
            st                  = collect_stats.cs(X[i], Y[j])
            stats_bundle        = st.stats()
            SB[:,k*7:(k+1)*7] = stats_bundle
        FINAL               = np.zeros((len(motifs), 7))
        '''
            compute average of stats
        '''
        FINAL[:,:4]=np.array([np.mean([ SB[:,k+i*7] for i in range(len(comparisons))],axis=0) for k in range(4)]).T
        '''
            combine pvalues fisher's method
        '''
        FINAL[:,4:]=1.0-ss.chi2(len(comparisons)*2).cdf(np.array([np.sum([ -2*np.log(SB[:,k+i*7]) for i in range(len(comparisons))],axis=0) for k in range(4,7)]).T)
        

        df  = pd.DataFrame(columns=("motif","counts", "distance", "mds", "ks_stat", "pval_distance", "pval_mds", "pval_ks") )
        for i,m in enumerate(motifs):
            df.loc[i]   = [m]+list(FINAL[i,:])
        df  = df.sort_values("pval_ks")
        return df
    





