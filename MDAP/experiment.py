import pandas as pd,seaborn as sns
import matplotlib.pyplot as plt
import numpy as np, seaborn as sns
import matplotlib as mpl,matplotlib.cm as cm
from scipy.spatial.distance import pdist,squareform
import scipy.cluster.hierarchy as sch
from collections import OrderedDict
import normal_laplace as nl
class dataset:
	def __init__(self, File):
		self.f 	= File
		self.df 	= None

		self.norm= False
		self._load_file()
	def _load_file(self):
		try:
			open(self.f,"r")
		except:
			raise TypeError, "couldn't open:"+self.f
		try:
			self.df 		= pd.read_csv(self.f)
			self.df 		= self.df.set_index(self.df.ID)
			self.motifs = self.df.ID
		except:
			raise TypeError, "couldn't parse:"+self.f
	def collapse(self, motif=None):
		if motif is not None:
			X 		= self.df[(self.df.ID==motif)].as_matrix()
		else:
			X 		= self.df.as_matrix()
		motifs= X[:,0]

		return np.array(X[:,1:],dtype=float),motifs
	def _compute_mds(self,x, percent=0.1):
		c 		= len(x)/2
		w 		= len(x)*percent/2
		return sum(x[c-w:c+w])/(sum(x)+1)
	def get_mds(self, motif=None, percent=0.1):
		X,motifs 		= self.collapse(motif=motif)
		mds 	= OrderedDict([(motifs[i], self._compute_mds(X[i,:],percent=percent)) for i in range(X.shape[0])])
		return mds
	def get_distribution(self, motif, unroll=False):
		ds 	= self.df[self.df.ID==motif]
		if not unroll:
			return ds
		return [d for i,d in enumerate(np.linspace(-1500,1500,len(ds))) for j in range(ds[i])  ]


	def get_barcode(self,motif, ax=None,tss=False,
									comb = False,
									cmap = sns.cubehelix_palette(as_cmap=True,dark=0, light=1),
									norm = None,
									bins = 150):
		X,motifs 	= self.collapse(motif=motif)

		sns.set(font_scale=2)
		sns.set_style("white")
		if ax is None:
			F 		= plt.figure(tight_layout=True,figsize=(10,2))
			ax 	= plt.gca()			

		counts,edges 	= np.histogram(range(X.shape[1]),weights=X[0,:],bins=bins)
		counts 			= np.array([c/float(max(counts)) for c in counts])
		if norm is None:
			norm 				= mpl.colors.Normalize(vmin=counts.min(), vmax=counts.max())
		m 					= cm.ScalarMappable(norm=norm, cmap=cmap)
		w 					= 3000.0/bins
		ax.bar(np.linspace(-1500,1500,bins),np.ones((bins,)), width=w,color=map(m.to_rgba, counts),edgecolor="white" )
		sns.despine(ax=ax,left=True,bottom=True)
		ax.set_xlabel("bidirectional origin")

		ax.set_yticks([]);ax.set_xticks([-1500,-750,0,750,1500])
		ax.set_xlim([-1500,1500])
		
	def get_histogram(self,motif, ax=None,bins=100,normed=0):
		X,motifs 	= self.collapse(motif=motif)
		sns.set(font_scale=2)
		sns.set_style("white")
		if ax is None:
			F 		= plt.figure(tight_layout=True,figsize=(10,10))
			ax 	= plt.gca()
		counts,edges 	= np.histogram(np.linspace(-1500,1500,X.shape[1]),weights=X[0,:],bins=bins,normed=int(normed))
		if normed:
			counts 		   = [float(c)/float(max(counts)) for c in counts]

		w 					= 3000.0/bins
		ax.bar((edges[:-1] + edges[1:])/2., counts,edgecolor="white",width=w,alpha=1)
		sns.despine(ax=ax)
		ax.set_xlim([-1500,1500])
		ax.set_ylabel("Frequency of Motif")
		ax.set_xlabel("bidirectional origin")


		
	def random_motif(self,ax=None, as_histogram = False ,
									cmap = sns.cubehelix_palette(as_cmap=True,dark=0, light=1),
									norm = None,
									bins = 150):
		motif 	= self.df.motif[np.random.randint(0,len(self.df))]
		if not as_histogram:
			self.get_barcode(motif,ax=ax,cmap=cmap, norm=norm,bins=bins)
		else:
			self.get_histogram(motif, ax=ax,bins=bins)
	def get_all(self,ax=None,
									cmap = sns.cubehelix_palette(as_cmap=True,dark=0, light=1),
									norm = None,
									bins = 150,clustering=False,sort_md=False):
		if ax is None:
			sns.set(font_scale=2)
			sns.set_style('white')
			F 	= plt.figure(tight_layout=True,figsize=(4,10))
			ax = plt.gca()	
		
		'''
			collapse data
		'''
		X,motifs 		= self.collapse()
		X 		= np.array([np.histogram(np.linspace(-1500,1500,len(x)),weights=x,bins=bins)[0] for x in X],dtype=float)
		X 		= np.array([x/max(x+1) for x in X])
		'''
			sort the rows?
		'''

		idx 	= range(X.shape[0])

		if clustering:
			Y 			= sch.linkage(X, method='ward')
			Z2 		= sch.dendrogram(Y, no_plot=True)
			idx 		= Z2['leaves']
		elif sort_md:
			mds 			= self.get_mds().values()
			idx 		= np.array(sorted(range(len(mds)), key=lambda k: mds[k],reverse=True))
		X 		= X[idx,:]
		ax.imshow(X,cmap=cmap,norm=norm,aspect="auto",interpolation="nearest")
		ax.set_yticks([])
		ax.set_xticks(np.linspace(0, X.shape[1], 5))
		ax.set_xticklabels(["-1500","-750","0","750","1500"])
		ax.set_xlabel("bidirectional origin")
		sns.despine(ax=ax,left=True,bottom=True)
		return motifs[idx]
	def show_motif_displacement(self,ax=None, All=False,
												Spec=False,comb=False, tss=False,
												as_histogram=False,
												Random=False,bins=150,
												cmap = sns.cubehelix_palette(as_cmap=True,dark=0, light=1),
												norm = None,clustering=False,sort_md=False,normed=0):
		stf = None
		if Spec:
			if not self.df.ID.isin([Spec]).any():
				print Spec,"not a valid motif identifier"
				print "try: the get_motifs() method to see the list..."
				return -1
			if as_histogram:
				stf 	= self.get_histogram(Spec, ax=ax,  bins=bins,normed=normed )
			else:
				stf 	= self.get_barcode(Spec, ax=ax,  bins=bins,cmap=cmap,norm=norm )

		elif Random:
			stf 	= self.random_motif(ax=ax,tss=tss,as_histogram=as_histogram,comb=comb,cmap=cmap,norm=norm,bins=bins)
		elif All:
			stf 	= self.get_all(ax=ax, bins=bins, cmap=cmap, norm=norm,clustering=clustering,sort_md=sort_md)
		return stf

	def uni_bi(self,comb=False, tss=False,motif=False,bins=150,norm=False):
		X,motifs 	= self.collapse(motif=motif, tss=tss, comb=comb)
		em 			= np.vectorize(lambda x: nl.EM().get_best(x,norm=norm) )
		self.norm 	= norm
		'''
			bin for speed?
		'''
		X 		= np.array([np.histogram(np.linspace(-1500,1500,X.shape[1]),weights=X[i,:],bins=bins)[0] for i in range(X.shape[0])],dtype=float)
		R 		= np.array(np.apply_along_axis(nl.EM().get_best,1,X))
		return R
	def draw_density_function(self,mu, weight, b, A=-1500,B=1500,ax=None):
		if ax is None:
			sns.set(font_scale=2)
			sns.set_style("white")

			ax 	= plt.gca()
		l 		= abs(B-A)
		xs 	= np.linspace(A,B,300)
		ys 	= [nl.EM().density(x,mu,weight, b,l,norm=self.norm ) for x in xs]
		ys 	= [y/max(ys) for y in ys]
		sns.despine(ax=ax)
		ax.plot(xs,ys)

	def get_motifs(self):
		return list(set(self.df.motif))





