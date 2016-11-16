import sys
sys.path.append("/Users/joazofeifa/Lab/motif_displacement_analysis/motif_displacement_analysis_package/")
import matplotlib.pyplot as plt
import numpy as np,math
from scipy.interpolate import spline
from scipy.special import erf
import time
from matplotlib import cm

import matplotlib as mpl
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
class PSSM:
    def __init__(self, name):
        self.name   = name
        self.P      = {}
        self.N      = {}
        self.D      = {}
    def insert_EXP(self, ID, p,n):
        if ID not in self.P:
            self.P[ID]    = list()
            self.N[ID]    = list()
        self.P[ID].append(p)
        self.N[ID].append(n)
    def insert_displacements(self, ID, x):
        if ID not in self.D:
            self.D[ID]  = list()
        self.D[ID].append(x)
def format_str(x,k=5):
    if len(x) > k:
        return x[:k]
    if len(x) < k:
        return x + " ".join(["" for i in range(k-len(x))])
    return x




class mds_frame:
    def __init__(self):
        self.motif_models   = {} #inserted motif models
        self.EXPS           = {} #inserted experiment IDS
        self.motifs         = list()
    def insert(self, motif, exp_ID, p, n):
        if motif not in self.motif_models:
            self.motif_models[motif]    = PSSM(motif)
        self.motif_models[motif].insert_EXP(exp_ID,p,n)
    def insert_displacements(self, motif, exp_ID, displacements):
        if motif not in self.motif_models:
            self.motif_models[motif]    = PSSM(motif)
        self.motif_models[motif].insert_displacements(exp_ID,displacements)
    def to_array(self):
        self.motifs=self.motif_models.keys()
        for motif in self.motif_models:
            for exp_ID in self.motif_models[motif].P:
                self.motif_models[motif].P[exp_ID] = np.array(self.motif_models[motif].P[exp_ID])
                self.motif_models[motif].N[exp_ID] = np.array(self.motif_models[motif].N[exp_ID])
    def show_barcodes(self, motif, EXPS,ax=None,BINS=100):
        P       = self.motif_models[motif]
        X       = np.zeros((len(EXPS), BINS))
        for i,exp in enumerate(EXPS):
            disps   = P.D[exp]
            for d in disps:
                counts,edges    = np.histogram(np.linspace(-1000,1000,len(d)), weights=d,bins=BINS)
                X[i,:]+=counts
            X[i,:]/=float(len(disps))
        if ax is None:
            ax  = plt.gca()

        ax.imshow(X, cmap=cm.GnBu,interpolation='nearest', aspect='auto',vmin=np.min(X), vmax=np.max(X))
        ax.set_xticks([])
        ax.set_yticks(range(len(EXPS)))
        ax.set_yticklabels(EXPS)

        despine(ax,left=True, bottom=True)


    def load_MD_score_file(self, FILES, exp_ID, P=2,DISPS=False):
        self.EXPS[exp_ID]   = 1
        for FILE in FILES:
            FH      = open(FILE,'r')
            C       = None
            for line in FH:
                if "Binned" in line and not DISPS:    
                    break
                elif "Binned" in line and DISPS:
                    C   = True
                elif C and "Empiracle" in line:
                    break
                elif C:
                    motif, scores  = line.strip("\n").split("\t")[:2]
                    scores          = map(float, scores.split(","))
                    self.insert_displacements(motif, exp_ID, scores)
                elif line[0]!="#":
                    line_array=line.strip("\n").split('\t')
                    motif, NS, scores = line_array[0], map(float, line_array[1].split(",")),map(float, line_array[2].split(","))    
                    self.insert(motif,exp_ID,scores[P], NS[P])
            FH.close()

    def compute_pvalue(self, motif,  b,a,mean=0):
        P1      = self.motif_models[motif].P[a]
        P2      = self.motif_models[motif].P[b]
        N1      = self.motif_models[motif].N[a]
        N2      = self.motif_models[motif].N[b]  
        if not sum(N1) or not sum(N2):
            return 0.5
        p       = (P1.dot(N1) + P2.dot(N2))  / (sum(N1)+sum(N2))
        SE      = math.sqrt(p*(1-p)*( (1.0/sum(N1)) + (1.0/sum(N2))  ))
        Z       = (np.mean(P2)-np.mean(P1)-mean) / SE
        pval    = 0.5*(1.0 + erf(Z / math.sqrt(2.0) ))
        return pval

    def dMDS(self, motif, b,a):

        return np.mean(self.motif_models[motif].P[b])- np.mean(self.motif_models[motif].P[a])
    def dN(self, motif, b,a ):
        return np.mean(self.motif_models[motif].N[a]) + np.mean(self.motif_models[motif].N[b])

    def differential_single(self, A,B, pval_threshold=pow(10,-4),ax=None,OUT="",
        xlabel=False,ylabel=False,legend=True,AX2=None,title=None):
        self.to_array()
        for i in (A,B):
            if len([1 for m in self.motifs if i not in self.motif_models[m].N])>0:
                print "Experiment Entry: \""+i+"\"is not in current mds_frame"
                print "loaded and available options are: " + ",".join(self.EXPS.keys())
                return False
        S       = False
        if ax is None:
            F   = plt.figure(figsize=(12,7),facecolor="white")
            ax  = F.add_axes([0.1,0.1,0.7,0.8])
            ax2  = F.add_axes([0.8,0.1,0.2,0.8])
            AX2     = True
            S      = True
        FHW     = None
        if OUT:
            try:
                FHW = open(OUT, "w")
            except:
                pass
        

        x     = [ self.dN(m, B,A) for m in self.motifs if "CPEB" not in m]
        y     = [ self.dMDS(m, B,A) for m in self.motifs if "CPEB" not in m]
        MEAN   = np.mean(y)
        z     = [ self.compute_pvalue(m,B,A,mean=MEAN) for m in self.motifs if "CPEB" not in m]


        ups   = [ (k,i,j,self.motifs[l]) for l,(i,j,k) in enumerate(zip(x,y,z)) if i > 30 and j > 0.05 and  k > (1.0 - pval_threshold)]
        downs = [ (k,i,j,self.motifs[l]) for l,(i,j,k) in enumerate(zip(x,y,z)) if i > 30 and j < -0.045 and  k < pval_threshold]
        nC    = [ (k,i,j,self.motifs[l]) for l,(i,j,k) in enumerate(zip(x,y,z)) if i > 30 and  pval_threshold < k < (1.0 - pval_threshold)]
        lbl   = r"$p-value<10^{" + str(int(math.log(pval_threshold,10))) + "}$"
        if len(ups):
            ax.scatter([i[1] for i in ups],[i[2] for i in ups],color="red",edgecolor="red",
                label=lbl)
        if len(downs):
            ax.scatter([i[1] for i in downs],[i[2] for i in downs],color="green",edgecolor="green",
                label=lbl)
        ax.scatter([i[1] for i in nC],[i[2] for i in nC],color="blue",edgecolor="blue",alpha=0.5)
        ax.set_xlim(10,max(x)+pow(10,5))
        if title is None:
            ax.set_title( A + r" vs " + B + r"",fontsize=20)
        else:
            ax.set_title(title, fontsize=20)
        if ylabel:
            ax.set_ylabel( r"$\Delta MDS$"  ,fontsize=30 )
        if xlabel:
            ax.set_xlabel( r"$N$" ,fontsize=30 )
        if legend and (len(ups) or len(downs)) :
            ax.legend(loc="best")
        ax.set_xscale("log")





        if FHW is not None:
            FHW.write("#Motif Model\tAverage Delta MD Score\tAverage N, base 10\tp-value\tComments\n")
            ups.sort()            
            for u in ups:
                FHW.write(format_str(u[-1],k=30)+"\t" + format_str(str(u[2]),k=10) 
                    + "\t" + format_str(str(math.log(u[1],10))) + "\t" + format_str(str(1.0-u[0])) + "\tUP\n" )

            downs.sort(reverse=True)
            for u in downs:
                FHW.write(format_str(u[-1],k=30)+"\t" + format_str(str(u[2]),k=10) 
                    + "\t" + format_str(str(math.log(u[1],10))) + "\t" + format_str(str(u[0])) + "\tDOWN\n" )
            for u in nC:
                FHW.write(format_str(u[-1],k=30)+"\t" + format_str(str(u[2]),k=10) 
                    + "\t" + format_str(str(math.log(u[1],10))) + "\t" + format_str(str(u[0])) + "\tNo Change\n" )
        NN=0
        if AX2:
            for i in np.arange(0,len(ups),2):
                lbl     = ",".join([x[-1].lstrip("HO_").split("_")[0] + "(" + str(x[2])[:5] + ")" for x in ups[i:i+2]])
                ax2.text(0,-NN, lbl,verticalalignment="top",color="red")
                NN+=1
            for i in np.arange(0,len(downs),2):
                lbl     = ",".join([x[-1].lstrip("HO_").split("_")[0] + "(" + str(x[2])[:5] + ")" for x in downs[i:i+2]])
                ax2.text(0,-NN, lbl,verticalalignment="top", color="green")
                NN+=1
            ax2.set_ylim(-NN,0)
            ax2.set_xlim(0,3)
            despine(ax2,bottom=True, left=True)        
            ax2.set_yticks([])
            ax2.set_xticks([])
            
        despine(ax)

        if S:
            plt.show()
    def differential_multiple(self, EXPS, base=None,
                                dt=False,smooth=False,
                                filter_static=False,ax=None,
                                pval_threshold=pow(10,-3),xlabel=True,ylabel=True,title="",
                                AX2=None,FIG=""):
        self.to_array()
        if base is not None and len(EXPS)!=len(base):
            print "You specified a base time series" 
            print "but this does not match in length your observed time series"
            return False
        for i in EXPS:
            if len([1 for m in self.motifs if i not in self.motif_models[m].N])>0:
                print "Experiment Entry: \""+i+"\"is not in current mds_frame"
                print "loaded and available options are: " + ",".join(self.EXPS.keys())
                return False
        S       = False
        if ax is None:
            F   = plt.figure(figsize=(15,6),facecolor="white")
            ax  = F.add_axes([0.1,0.4,0.8,0.4])
            ax2 = F.add_axes([0.1,0.1,0.8,0.2])
            S   = True
        ax.set_title(title)
        xs      = range(len(EXPS))
        if base:
            xs      = range(len(EXPS)+1)

        lines   = list()
        SIG     = {}
        for m in self.motifs:
            if base is not None:
                ys      = [0] + [ self.dMDS(m, EXPS[i], base[i]) for i in range(len(EXPS))]
                pvals   = [0.5]+ [ self.compute_pvalue(m, EXPS[i],base[i]) for i in range(len(EXPS))]
            elif dt:
                ys      = [0] + [ self.dMDS(m, EXPS[i+1],EXPS[i]) for i in range(len(EXPS)-1)]
                pvals   = [0.5] + [ self.compute_pvalue(m,EXPS[i+1],EXPS[i]) for i in range(len(EXPS)-1)]
            else:
                ys      = [ self.dMDS(m, EXPS[i],EXPS[0]) for i in range(len(EXPS))] 
                pvals   = [ self.compute_pvalue(m,EXPS[i],EXPS[0]) for i in range(len(EXPS)-1)]
            for i,p in enumerate(pvals):
                if i not in SIG:
                    SIG[i]  = list()
                if (p > (1.0 - pval_threshold) or p < pval_threshold) and abs(ys[i]) > 0.0025:
                    SIG[i].append((ys[i],m.lstrip("HO_").split('_')[0]))

            X,Y     = xs,ys
            if smooth:
                x_sm,y_sm= np.array(xs), np.array(ys)
                x_smooth = np.linspace(x_sm.min(), x_sm.max(), 100)
                y_smooth = spline(xs, ys, x_smooth,order=3)
                X,Y      = x_smooth,y_smooth
            if not filter_static or max(pvals) > (1.0 - pval_threshold) or min(pvals)<pval_threshold:
                if (max(pvals) > (1.0 - pval_threshold) or min(pvals)<pval_threshold) and abs(ys[i]) > 0.0025: 
                    l=ax.plot(X,Y,lw=1.0,color="red",alpha=1.0)
                else:
                    l=ax.plot(X,Y,lw=1.0,color="blue",alpha=0.07)

        ax.set_xticks(xs)
        if AX2:
            ax2.set_xticks(xs)
            ax2.set_yticks([])
            ax2.set_xticklabels([])
            ax2.set_yticklabels([])

            for i in SIG:
                SIG[i].sort()
                labels  = [y + "(" + str(x)[:5] + ")" for x,y in SIG[i]]    
                ax2.text(i-0.25,0, "\n".join([",".join( labels[i:i+2] )  for i in np.arange(0,len(labels),2)]),verticalalignment="top")
            ax2.set_ylim(-10,1)
            despine(ax2,left=True,bottom=True)
        if ylabel:
            ax.set_ylabel( r"$\Delta MDS(k-j)$"  ,fontsize=30 )

        if xlabel:
            if dt:
                ax.set_xticklabels(["Start"] +[ EXPS[i+1] + "-" + EXPS[i] for i in range(len(EXPS)-1)],rotation=12,fontsize=15)
            else:
                ax.set_xticklabels([ EXPS[i] + "-" + EXPS[0] for i in range(len(EXPS))],rotation=12,fontsize=15)            
        else:
            ax.set_xticks([])
        despine(ax)
        if FIG:
            plt.savefig(FIG)
        if S:
            plt.show()
        return SIG










if __name__ == "__main__":
    DIR     = "/Users/joazofeifa/Lab/new_motif_distances/motif_hits_mouse/"
    DIR     = "/Users/joazofeifa/Lab/new_motif_distances/motif_hits_human/"
    def add(b,a=DIR , c="_enrichment_stats.tsv"):
        return a+b+c
    

    MDS3    = mds_frame()

    MDS3.load_MD_score_file(map(add,["SRR1015583", "SRR1015584"]), "0(min)")
    MDS3.load_MD_score_file(map(add,["SRR1015585", "SRR1015586"]), "10(min)")
    MDS3.load_MD_score_file(map(add,["SRR1015587", "SRR1015588"]), "30(min)")
    MDS3.load_MD_score_file(map(add,["SRR1015589", "SRR1015590"]), "120(min)")
    
    MDS3.differential_single("0(min)", "30(min)")

    # MDS3.differential_multiple(["0(min)","10(min)", 
    #     "30(min)","120(min)"  ],
    #     smooth=True,filter_static=True,pval_threshold=pow(10,-5),dt=True,title="Time\n(TNF-alpha Treatment)\nCell Type: AC16\n(Luo,2016)")




