import numpy as np
import matplotlib.pyplot as plt,math
import seaborn as sns,scipy.stats as ss
from scipy.special import erf
def check_pval(pv):
	if math.isnan(pv) or math.isinf(pv):
		return 0.5
	return pv

class cs:
	def __init__(self, X,Y,a=-1500,b=1500,percent=0.1):
		'''
			x and y come in as histograms
			and they should have the same dimension;

			they are really just two experiments
		'''
		assert X.shape == Y.shape
		self.X 		= X
		self.Y 		= Y
		self.xs 		= np.abs(np.linspace(a,b,self.X.shape[1]))
		self.xs2 	= np.linspace(a,b,self.X.shape[1])**2
		self.bools 	= np.array([int(i < (b-a)*percent*0.5) for i in self.xs ])
	def ks(self, x,y):
		xu 	= [d for i,d in enumerate(self.xs) for j in range(int(x[i]))]
		yu 	= [d for i,d in enumerate(self.xs) for j in range(int(y[i]))]
		score,pv 	= ss.ks_2samp(xu,yu)
		return score,pv
	def compute_pvalues(self,xy):
		'''
			should get 3 pvalues and 7 stats

			0) N (X and Y)
			1) MD score (X and Y)
			2) MD mean  (X and Y)
			3) KS_stat			
		'''

		x,y 	= xy[:len(xy)/2],xy[len(xy)/2:]
		nX 	= sum(x)+1
		nY 	= sum(y)+1

		mX 	= self.xs.dot(x)/nX 
		vX 	= self.xs2.dot(x)/nX 
		pX 	= self.bools.dot(x)

		mY 	= self.xs.dot(y)/nY 
		vY 	= self.xs2.dot(y)/nY 
		pY 	= self.bools.dot(y)


		'''
			KS -test
		'''
		try:
			ks_stat, ks_pval 	= self.ks(x,y)
		except:
			ks_stat, ks_pval 	= 0.5,0.5

		'''
			z-test on proportions
		'''
		PX      = pX/nX
		PY      = pY/nY
		p       = (PX*nX  + PY*nY)  / ( nX+nY )
		SE      = math.sqrt(p*(1-p)*( (1.0/ nX) + (1.0/nY)  ))
		Z       = ( PX - PY   ) / SE
		if Z > 0:
			mds_pval    = 2*(1.0-(0.5+0.5*erf(Z / math.sqrt(2.0) )))
		else:
			mds_pval    = 2*(0.5+0.5*erf(Z / math.sqrt(2.0) ))

		'''
			t-test on means
		'''
		SE 		= math.sqrt( (nX*vX + nY*vY)/(nX+nY+1) )*math.sqrt((1.0/nX) + (1.0/nY))
		DF 		= nX + nY -2
		if mX > mY:
			mean_pval 		= 2.0-2*ss.t.cdf((mX-mY)/SE,DF)
		else:
			mean_pval 		= 2*ss.t.cdf((mX-mY)/SE,DF)
		return nX+nY, mY-mX, PY-PX, ks_stat, check_pval(mean_pval), check_pval(mds_pval), check_pval(ks_pval)



	def stats(self,mean=True, MDS=True, KS=True):
		'''
			pv_stats is 641 x 10 column

			rows are motif models
			columns in this order
			1) N (control) + N (treatment) 
			2) mean (treatment) - mean (control)
			3) mds  (treatment) - mds  (control)
			4) ks difference (difference in empiracle cdfs)
			5) pvalue mean
			6) pvalue mds
			7) pvalue ks
			
		'''
		AP 		= np.apply_along_axis
		XY 		= np.concatenate([self.X.T,self.Y.T]).T
		pvs_stats= AP(self.compute_pvalues, 1, XY)
		return pvs_stats
