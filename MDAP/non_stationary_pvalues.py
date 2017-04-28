import numpy as np,pandas as pd
import matplotlib.pyplot as plt,seaborn as sns
import experiment,scipy.stats as ss

class simulation:
	def __init__(self, FILE):
		self.FILE 	= FILE
		self.G 		= {}
		self.GC 		= list()
		self._load_DB()
	def _load_DB(self):
		FH,collect 	= open(self.FILE,"r"),False
		for line in FH:
			if line[0]==">":
				ID 		  = line[1:].strip('\n')
			elif line[0]=="~":
				self.G[ID] = map(float,line.strip("\n").split('|')[1].split(','))
			elif "Est" in line[:4]:
				collect=True
			elif collect:
				self.GC.append(map(float,line.strip("\n").split("\t")[-1].split(",")))
		FH.close()
		self.GC 	= np.array(self.GC)
	def draw_GC(self):
		F 		= plt.figure(tight_layout=True)
		sns.set(font_scale=2);sns.set_style("ticks")
		ax 	= plt.gca()
		char 	= "ACGT"
		for j in range(self.GC.shape[1]):
			ax.plot(np.linspace(-1500,1500,self.GC.shape[0]), self.GC[:,j],label=char[j])
		ax.legend(loc="best")
		sns.despine()
		plt.show()
	def _compute_mds(self,x, percent=0.1):
		c 		= len(x)/2
		w 		= len(x)*percent/2
		return sum(x[c-w:c+w])/float(sum(x)+1)
	def make_draws(self,motif, N,BSN,percent=0.1):
		assert motif in self.G
		vals 		= np.array(self.G[motif])/sum(self.G[motif])
		draws 	= np.random.multinomial(N,vals,size=BSN)
		mds 		= [self._compute_mds(draws[i,:],percent=percent) 
						for i in range(draws.shape[0])]
		mean,std = np.mean(mds),np.std(mds)
		return mean,std
	def get_pvalue(self,obs,mean,std):
		if obs  > mean:
			return (1.0-ss.norm(mean,std).cdf(obs))*2
		return ss.norm(mean,std).cdf(obs)*2
def main():
	SRR 		= "SRR1105737"
	MDS_DIR 	= "../motif_displacement_files/"
	DB 		= MDS_DIR + "human_non.db"
	d 			= experiment.dataset(MDS_DIR + SRR + "_MDS.csv")
	S 			= simulation(DB)
	BSN 		= 1000
	percent 	= 0.1
	OUT  		= pd.DataFrame(columns=["ID", "MDS", "mean", "std", "pvalue"])
	'''
		compute MD scores
	'''
	OBS 		= d.get_mds(percent=percent)
	'''
		iterate across motifs in d
	'''
	for k,m in enumerate(d.motifs):
		print k+1,"/",len(d.motifs)
		N 				= sum(d.df.loc[m][1:])
		mean,std 	= S.make_draws(m,N,BSN,percent=percent)
		obs 			= OBS[m]
		pval 			= S.get_pvalue(obs,mean,std)
		OUT.loc[k] 	= [m,obs,mean,std,pval]
	OUT.to_csv(MDS_DIR+SRR + "_ns.csv",index=False)




if __name__ == "__main__":
	main()