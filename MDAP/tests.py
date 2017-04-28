import experiment,seaborn as sns, matplotlib.pyplot as plt
import numpy as np,os
from md_stats import mds_frame
import pandas as pd
from display import show
def test_single_cases():
	'''
		change this file accordingly
	'''
	SRR 		= "SRR1105737"


	MDS_DIR 	= "../motif_displacement_files/"

	d 			= experiment.dataset(MDS_DIR + SRR + "_MDS.csv")
	sns.set(font_scale=2)
	sns.set_style("ticks")
	F  		= plt.figure(figsize=(15,6))
	ax 		= F.add_axes([0.6,0.3,0.34,0.6])
	ax2 		= F.add_axes([0.6,0.035,0.34,0.06])
	d.show_motif_displacement(Spec="HO_NRF1_HUMAN.H10MO.A",ax=ax,
											as_histogram=True,norm=False)

	ax.set_xlabel("\ndistance to $\mu$ (base pairs)" )
	d.show_motif_displacement(Spec="HO_NRF1_HUMAN.H10MO.A",ax=ax2)
	ax2.set_xlabel("")
	ax2.set_xticks([])

	ax 		= F.add_axes([0.1,0.3,0.34,0.6])
	ax2 		= F.add_axes([0.1,0.035,0.34,0.06])

	sns.set_style("ticks")
	d.show_motif_displacement(Spec="HO_MEOX1_HUMAN.H10MO.D",ax=ax,as_histogram=True,norm=False,)
	ax.set_xlabel("\ndistance to $\mu$ (base pairs)")
	d.show_motif_displacement(Spec="HO_MEOX1_HUMAN.H10MO.D",ax=ax2)
	ax2.set_xlabel("")
	ax2.set_xticks([])

	plt.show()
def test_all():
	'''
		change this file accordingly
	'''
	SRR 		= "SRR1105737"


	MDS_DIR 	= "../motif_displacement_files/"

	d 			= experiment.dataset(MDS_DIR + SRR + "_MDS.csv")
	sns.set(font_scale=2)
	sns.set_style("ticks")

	F  		= plt.figure(figsize=(5,10))
	ax1 		= F.add_subplot(1,1,1)
	ax1.set_title("Ranked by MD-score")
	motif_order = d.show_motif_displacement(All=True,ax=ax1,sort_md=True,
		clustering=False,bins=150)
	plt.show()
def test_differential():
	RECOMP 	= True
	MDS_DIR 	= "../motif_displacement_files/"
	f1,f2 	= "SRR1015583", "SRR1015587"
	OUT 		= MDS_DIR+"differential_" + f1 + "_" + f2 + ".csv"
	if RECOMP:
		mds 		= 	mds_frame()
		mds.load_MD_score_file(MDS_DIR+f1+"_MDS.csv", "WT")
		mds.load_MD_score_file(MDS_DIR+f2+"_MDS.csv", "KO")
		df 		= mds.differential( "WT", "KO",h=150)
		df.to_csv(OUT,index=False)

	df 			= pd.read_csv(OUT)
	sns.set(font_scale=2)
	sns.set_style("ticks")
	F 				= plt.figure(tight_layout=True,figsize=(12,10))
	
	ax1,ax2 		= F.add_subplot(1,2,1),F.add_subplot(1,2,2)
	S 				= show(df=df,ax=ax2,FDR=pow(10,-5))
	ax1.hist(df.pval_mds,bins=60,color="steelblue",edgecolor="white",label="N="+str(df.shape[0]))
	ax1.set_xlabel("p-value");ax1.set_ylabel("frequency")
	ax1.legend(loc="best")
	sns.despine()
	plt.show()
if __name__ == "__main__":
	#test_single_cases()
	#test_all()
	test_differential()