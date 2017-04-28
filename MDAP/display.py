import os,pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import numpy as np
def despine(ax,right=True, top=True,left=False, bottom=False):
	if right:
	  ax.yaxis.set_ticks_position('left')
	  ax.spines['right'].set_visible(False)
	if left:
	  ax.yaxis.set_ticks_position('left')
	  ax.spines['left'].set_visible(False)
	if bottom:
	  ax.xaxis.set_ticks_position('bottom')
	  ax.spines['bottom'].set_visible(False)
	if top:
	  ax.spines['top'].set_visible(False)
	  ax.xaxis.set_ticks_position('bottom')
class show:
	def __init__(self, df=None, path=None,ax=None, FDR=pow(10,-1), N_by_delta=True,MD_by_MD=False,delta_HIST=False):
		if df is None and path is None:
			print "Please specific either pandas dataframe or path"
			return None

		self.ax 			= ax
		self.N_by_delta= N_by_delta
		self.MD_by_MD	= MD_by_MD
		self.delta_HIST= delta_HIST
		self.FDR 		= FDR
		if df is not None:
			self.df 		= df
		else:
			try:
				self.df 	= pd.read_csv(path)
			except:
				print "could not open: ", path
				return None
		self.display()
	def display(self):
		self.display_single()


	def display_histogram(self, df,ax,ylim=None):
		MAX 				= max([abs(x) for x in df.mds])
		counts,edges 	= np.histogram(df.mds, bins=40, range=(-MAX,MAX))
		ax.barh((edges[1:] + edges[:-1])*0.5,counts,height=MAX*2/40.0,edgecolor="white"	)
		ax.set_xlabel("Frequency")
		if ylim:
			ax.set_ylim(ylim)
	def display_MA(self, df, ax):
		current_palette 	= sns.color_palette()
		df_non 				= df[df.pval_mds> self.FDR]
		df_sig 				= df[df.pval_mds< self.FDR]

		ax.scatter(df_non.counts+1, df_non.mds,edgecolor="white",color=current_palette[0])	
		ax.scatter(df_sig.counts+1, df_sig.mds,edgecolor="white",color=current_palette[2])	
		ax.set_xlim(10,max(df_non.counts+1))	
		ax.set_xlabel("Total Motifs between A and B")
		ax.set_ylabel("MD Score(B)-MD Score(A)")
		ax.set_xscale("log")
		despine(ax)


	def display_single(self):
		df 				= self.df
		if self.ax is not None:
			if self.N_by_delta:
				print "Warning: User specified an axis will show N by delta MD scatter"
				self.display_MA(df,self.ax)
			elif self.delta_HIST:
				print "Warning: User specified an axis will show histogram of delta MD scores"
				self.display_histogram(df,self.ax)
		else:
			sns.set(font_scale=2)
			sns.set_style("white")
			F 				= plt.figure(facecolor='white', figsize=(15,5),tight_layout=True)
			ax1,ax2 		= F.add_subplot(1,2,1),F.add_subplot(1,2,2)
			self.display_MA(df, ax1)
			self.display_histogram(df, ax2,ylim=ax1.get_ylim())
			for ax in (ax1,ax2):
				despine(ax)
			ax2.set_yticks([])
			plt.show()

if __name__ == "__main__":
	show("test.csv")