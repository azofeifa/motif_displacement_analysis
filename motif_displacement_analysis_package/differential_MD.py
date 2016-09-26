import matplotlib.pyplot as plt
import numpy as np,math
from scipy.interpolate import spline
from scipy.special import erf
import time
import matplotlib as mpl
from load import despine


class PSSM:
    def __init__(self, name):
        self.name   = name
        self.P      = {}
        self.N      = {}
    def insert_EXP(self, ID, p,n):
        if ID not in self.P:
            self.P[ID]    = list()
            self.N[ID]    = list()
        self.P[ID].append(p)
        self.N[ID].append(n)
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
    def to_array(self):
        self.motifs=self.motif_models.keys()
        for motif in self.motif_models:
            for exp_ID in self.motif_models[motif].P:
                self.motif_models[motif].P[exp_ID] = np.array(self.motif_models[motif].P[exp_ID])
                self.motif_models[motif].N[exp_ID] = np.array(self.motif_models[motif].N[exp_ID])
    def load_MD_score_file(self, FILES, exp_ID, P=2):
        self.EXPS[exp_ID]   = 1
        for FILE in FILES:
            FH      = open(FILE,'r')
            for line in FH:
                if "Binned" in line:    break
                elif line[0]!="#":
                    line_array=line.strip("\n").split('\t')
                    motif, NS, scores = line_array[0], map(float, line_array[1].split(",")),map(float, line_array[2].split(","))    
                    self.insert(motif,exp_ID,scores[P], NS[P])
            FH.close()

    def compute_pvalue(self, motif,  b,a):
        P1      = self.motif_models[motif].P[a]
        P2      = self.motif_models[motif].P[b]
        N1      = self.motif_models[motif].N[a]
        N2      = self.motif_models[motif].N[b]  
        if not sum(N1) or not sum(N2):
            return 0.5
        p       = (P1.dot(N1) + P2.dot(N2))  / (sum(N1)+sum(N2))
        SE      = math.sqrt(p*(1-p)*( (1.0/sum(N1)) + (1.0/sum(N2))  ))
        Z       = (np.mean(P2)-np.mean(P1)) / SE
        pval    = 0.5*(1.0 + erf(Z / math.sqrt(2.0) ))
        return pval

    def dMDS(self, motif, b,a):

        return np.mean(self.motif_models[motif].P[b])-np.mean(self.motif_models[motif].P[a])
    def dN(self, motif, b,a ):
        return np.mean(self.motif_models[motif].N[a]) + np.mean(self.motif_models[motif].N[b])

    def differential_single(self, A,B, pval_threshold=pow(10,-5),ax=None,OUT=""):
        self.to_array()
        for i in (A,B):
            if len([1 for m in self.motifs if i not in self.motif_models[m].N])>0:
                print "Experiment Entry: \""+i+"\"is not in current mds_frame"
                print "loaded and available options are: " + ",".join(self.EXPS.keys())
                return False
        
        if ax is None:
            F   = plt.figure(figsize=(12,7))
            ax  = plt.gca()
        FHW     = None
        if OUT:
            try:
                FHW = open(OUT, "w")
            except:
                pass
        

        x     = [ self.dN(m, B,A) for m in self.motifs]
        y     = [ self.dMDS(m, B,A) for m in self.motifs]
        z     = [ self.compute_pvalue(m,B,A) for m in self.motifs]


        ups   = [ (k,i,j,self.motifs[l]) for l,(i,j,k) in enumerate(zip(x,y,z)) if i > 30 and  k > (1.0 - pval_threshold)]
        downs = [ (k,i,j,self.motifs[l]) for l,(i,j,k) in enumerate(zip(x,y,z)) if i > 30 and  k < pval_threshold]
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
        ax.set_title( A + r"$(j)$ vs" + B + r"$(k)$",fontsize=20)
        ax.set_ylabel( r"$\Delta MDS(k-j)$"  ,fontsize=30 )
        ax.set_xlabel( r"$(n_k + n_j)$" ,fontsize=30 )
        ax.legend(loc="best")
        ax.set_xscale("log")

        print "-------------------------------------------------"
        print "UPs      (" + str(len(ups))+") -> " + ", ".join([ "(+"+ str(i[2])[:5] + ")" + i[-1]  for i in ups])
        print "Downs    (" + str(len(downs))+") -> " + ", ".join(["("+ str(i[2])[:6] + ")" + i[-1] for i in downs])
        print "-------------------------------------------------"

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
        despine(ax)
        plt.tight_layout()
        plt.show()
    def differential_multiple(self, EXPS,
                                dt=False,smooth=False,
                                filter_static=False,ax=None,pval_threshold=pow(10,-3)):
        self.to_array()
        for i in EXPS:
            if len([1 for m in self.motifs if i not in self.motif_models[m].N])>0:
                print "Experiment Entry: \""+i+"\"is not in current mds_frame"
                print "loaded and available options are: " + ",".join(self.EXPS.keys())
                return False
        if ax is None:
            F   = plt.figure(figsize=(15,6))
            ax  = plt.gca()
        xs      = range(len(EXPS))
        lines   = list()
        for m in self.motifs:
            if dt:
                ys      = [0] + [ self.dMDS(m, EXPS[i+1],EXPS[i]) for i in range(len(EXPS)-1)]
                pvals   = [ self.compute_pvalue(m,EXPS[i+1],EXPS[i]) for i in range(len(EXPS)-1)]
            else:
                ys  = [ self.dMDS(m, EXPS[i],EXPS[0]) for i in range(len(EXPS))] 
                pvals   = [ self.compute_pvalue(m,EXPS[i],EXPS[0]) for i in range(len(EXPS)-1)]
            X,Y     = xs,ys
            if smooth:
                x_sm,y_sm= np.array(xs), np.array(ys)
                x_smooth = np.linspace(x_sm.min(), x_sm.max(), 100)
                y_smooth = spline(xs, ys, x_smooth,order=2)
                X,Y      = x_smooth,y_smooth
            if not filter_static or max(pvals) > (1.0 - pval_threshold) or min(pvals)<pval_threshold:
                if max(pvals) > (1.0 - pval_threshold) or min(pvals)<pval_threshold: 
                    l=ax.plot(X,Y,lw=1.0,color="blue",alpha=1.0)
                else:
                    l=ax.plot(X,Y,lw=1.0,color="blue",alpha=0.2)
                
        plt.xticks(xs)
        plt.ylabel( r"$\Delta MDS(k-j)$"  ,fontsize=30 )
        if dt:
            ax.set_xticklabels(["Start"] +[ EXPS[i+1] + "-" + EXPS[i] for i in range(len(EXPS)-1)],rotation=10,fontsize=15)
        else:
            ax.set_xticklabels([ EXPS[i] + "-" + EXPS[0] for i in range(len(EXPS))],rotation=10,fontsize=15)            
        despine(ax)
        plt.tight_layout()
        plt.show()











if __name__ == "__main__":
    DIR     = "/Users/joazofeifa/Lab/new_motif_distances/motif_hits_mouse/"
    #DIR     = "/Users/joazofeifa/Lab/new_motif_distances/motif_hits_human_2/"
    def add(b,a=DIR , c="_enrichment_stats.tsv"):
        return a+b+c
    
    SRRS    = ["SRR1145801", "SRR1145808","SRR1145815","SRR1145822", "SRR1145829"]  
    


    # MDS     = mds_frame()
    
    # MDS.load_MD_score_file(map(add,["SRR1105737","SRR1105740"]), "HCT116 DMSO")
    # MDS.load_MD_score_file(map(add,["SRR1105739","SRR1105738"]), "HCT116 Nutlin")


    # MDS2     = mds_frame()
    
    # MDS2.load_MD_score_file(map(add,["SRR1015583","SRR1015584"]), "AC16 DMSO")
    # MDS2.load_MD_score_file(map(add,["SRR1015585","SRR1015586"]), "AC16 TNFalpha(10m)")
    # MDS2.load_MD_score_file(map(add,["SRR1015587","SRR1015588"]), "AC16 TNFalpha(30m)")
    # MDS2.load_MD_score_file(map(add,["SRR1015589","SRR1015590"]), "AC16 TNFalpha(1hr)")




    MDS3    = mds_frame()

    MDS3.load_MD_score_file(map(add,["SRR930691", "SRR930692", "SRR930693", "SRR930694","SRR930684"]), "KLA_DMSO")
    MDS3.load_MD_score_file(map(add,["SRR930670"]), "KLA_10min") 
    MDS3.load_MD_score_file(map(add,["SRR930671", "SRR930659","SRR930660","SRR930661","SRR930662"]), "KLA_1hr") 
    MDS3.load_MD_score_file(map(add,["SRR930672", "SRR930663","SRR930664"]), "KLA_6hr") 
    MDS3.load_MD_score_file(map(add,["SRR930673", "SRR930665","SRR930666"]), "KLA_12hr") 
    MDS3.load_MD_score_file(map(add,["SRR930674", "SRR930667","SRR930668"]), "KLA_24hr") 


    MDS3.differential_multiple(["KLA_DMSO","KLA_10min", 
        "KLA_1hr","KLA_6hr","KLA_12hr","KLA_24hr"  ],smooth=True,filter_static=True,pval_threshold=pow(10,-10),dt=False)




    # labels  = SRRS
    # dlabels = [SRRS[i] +"-" + SRRS[i+1] for i in range(len(labels)-1)]
    # SRRS    = map(load_MD_score_file, [DIR+x+"_enrichment_stats.tsv" for x in SRRS])


    # show_time(SRRS, deriv=True, labels=labels,smooth=True, P=0, FILTER_STATIC=True)


