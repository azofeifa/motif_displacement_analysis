from math import log
import numpy as np
import math
import matplotlib.pyplot as plt
def despine(ax,right=True, top=True):
    if right:
        ax.yaxis.set_ticks_position('left')
        ax.spines['right'].set_visible(False)
    if top:
        ax.spines['top'].set_visible(False)
        ax.xaxis.set_ticks_position('bottom')

def draw_A(width, height,x,y, ax, color="blue",lw=2.0,alpha=0.5):
    ax.plot([x, x+width*0.5 ], [ y, y+height ],lw=lw, color=color ,alpha=alpha)
    ax.plot([ x+width*0.5, x+width ], [ y+height, y ],lw=lw, color=color ,alpha=alpha)
    ax.plot([x+width*0.25,x+width*0.75], [height*0.5, height*0.5], color=color, lw=lw,alpha=alpha)
    return y+height
def draw_T(width, height,x,y, ax, color="red",lw=2.0,alpha=0.5):
    center  = width*0.5
    ax.plot([x+center, x+center], [y, y+height],lw=lw, color=color,alpha=alpha)
    ax.plot([x, x+width], [y+height, y+height],lw=lw, color=color,alpha=alpha)
    return y+height
def draw_C(width, height,x,y, ax, color="green",lw=2.0,alpha=0.5):
    width/=2.0
    x       += width
    y       -=math.sin(3*math.pi/2.)*height
    ts      = np.linspace(3.85*math.pi/2., math.pi/6., 100)
    xs      = [ math.cos(t)*width+x for t in ts]
    ys      = [ math.sin(t)*height+y for t in ts]
    ax.plot(xs,ys,lw=lw, color=color,alpha=alpha)

    return y + math.sin(math.pi/2.)*height 
def draw_G(width, height,x,y, ax, color="orange",lw=2.0,alpha=0.5):
    width/=2.0
    x       += width

    y       += math.sin(math.pi/2.)*height
    ts      = np.linspace(3.95*math.pi/2., math.pi/6., 100)
    xs      = [ math.cos(t)*width+x for t in ts]
    ys      = [ math.sin(t)*height+y for t in ts]

    ax.plot(xs,ys,lw=lw, color=color,alpha=alpha)
    ptl     = math.cos(3*math.pi/2.)*width+x,  math.sin(3*math.pi/2.)*height+y

    ptu     = xs[0],ys[0]

    ax.plot([ ptu[0],ptu[0] ], [ptu[1], ptl[1] ],lw=lw, color=color,alpha=alpha)
    ax.plot([ ptu[0]-width*0.75,ptu[0] ], [ptu[1], ptu[1] ],lw=lw, color=color,alpha=alpha)
    return y + math.sin(math.pi/2.)*height
def draw_seq_log(PSSM,ax,lw=3.0,alpha=1):
    A       = np.array(PSSM)
    for i in range(A.shape[0]):
        methods     = (draw_A, draw_C, draw_G, draw_T)
        x           = i
        width       = 0.9
        y           = 0.0
        for j in range(4):
            d       = methods[j]
            h       = A[i,j]
            if j == 1 or j ==2:
                h/=2.0
            y       = d(width, h , x,y,ax,lw=lw,alpha=alpha)
    ax.set_xlim(0,A.shape[0])
    ax.set_xticks([])
    ax.set_xlabel("Position",fontsize=15)
    ax.set_ylabel("Probability",fontsize=15)

            
    pass
class PSSM:
    def __init__(self, model,
                 background=np.array([0.25,0.25,0.25,0.25]),
                 bins=100,pval_threshold=pow(10,-6),name=""):
        self.pval_threshold   = pval_threshold
        self.model,self.bins  = model,bins
        self.bg,self.name     = background,name
        self.ePDF, self.eCDF  = None,None
        self.ll_threshold     = None
    def draw_logo(self,ax=None):
        if ax is None:
            ax= plt.gca()
        draw_seq_log(self.model,ax )
        despine(ax)
        plt.show()
    def compute_binned_llr_distribution(self):
        N  = self.model.shape[0]
        R  = int(N*self.bins)
        C  = list()
        for i in range(N): #DP Programming ALG
            vals   = [ 2*(log(self.model[i,j]) - log(self.bg[j])) for j in range(4)]
            vals.sort()
            r      = min(pow(4,i+1), R)
            if i == 0:
                C  = [[v,1] for v in vals]
            else:
                min_x, max_x = vals[0] + C[0][0] , vals[-1]+C[-1][0]
                nC           = [[c,0] for c in np.linspace(min_x,max_x,r) ]
                for v in vals:
                    k=0
                    for c in C:
                        while k < r and nC[k][0]< (v+c[0]) :
                            k+=1
                        if k < r:
                            nC[k][1]+=c[1]
                C=nC
        C=np.array(C)
        TOTAL=float(sum(C[:,1]))
        self.ePDF, self.eCDF= C,np.array([[C[i,0],sum(C[:i,1]/TOTAL ) ]  for i in range(C.shape[0]) ])
        self.get_ll_threshold()
    def get_ll_threshold(self):
        i=0
        while i < self.eCDF.shape[0] and 1.0-self.eCDF[i,1] > self.pval_threshold:
            i+=1
        if i < self.eCDF.shape[0]:
            self.ll_threshold = self.eCDF[i,0]
        else:
            self.ll_threshold = self.eCDF[-1,0]        
    def scan(self, seq):
        ALPHABET  = {"A":0, "C":1,"G":2, "T":3}
        ALPHABETr = {"A":3, "C":2,"G":1, "T":0}
        forward   = [ ALPHABET[s] for s in seq]
        reverse   = [ ALPHABETr[s] for s in seq]        
        N         = self.model.shape[0]
        positions = list()
        ll        = 0
        for i in range(0,len(seq)-N):
            j,k    = 0,N-1
            llf,llr= 0.0,0.0
            while k >= 0:
                llf+=2*(log(self.model[j,forward[i+j] ])   - log(self.bg[forward[i+j]]))
                llr+=2*(log(self.model[k,reverse[i+N-j] ]) - log(self.bg[forward[i+N-j]]))
                j+=1
                k-=1
            if llf> self.ll_threshold:
                positions.append(i)
            if llr>self.ll_threshold:
                positions.append(i)
        return positions
        