import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
import sys
sys.path.append("/Users/joazofeifa/Lab/motif_displacement_analysis/motif_displacement_analysis_package/")

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

def show_motif_displacement(G, motif_model,title=""):
    if motif_model not in G:
        print ("motif_model: " + motif_model + " not in load_motif_displacement_raw structure")
        return False
    F  = plt.figure(figsize=(15,7))
    axh= F.add_axes([0.1,0.2,0.8,0.6])
    axb= F.add_axes([0.1,0.1,0.8,0.05])

    axh.set_title(title+motif_model, fontsize=25)

    bins         = 100
    counts,edges = np.histogram(np.linspace(-1500,1500,len(G[motif_model])), weights=G[motif_model],bins=bins )
    edges        = (edges[1:] + edges[:-1])/2.
    axh.bar(edges, counts,width=3000.0/bins,edgecolor="white",color="green")

    axb.set_yticks([])
    
    
    despine(axh)

    for ax in (axh, axb):
        ax.set_xlim(-1500,1500)
    axb.set_xlabel("Displacement From BTE", fontsize=15)
    norm = mpl.colors.Normalize(vmin=0, vmax=max(counts))
    cmap = cm.GnBu
    M    = cm.ScalarMappable(norm=norm, cmap=cmap)
    
    colors = [M.to_rgba(xi) for xi in counts ]
    axb.bar(np.linspace(-1500,1500,len(counts)),
            np.ones((len(counts), )), color=colors, edgecolor=colors,width=3000.0/len(counts) )
    axb.set_ylim(0,1)
    axh.set_ylabel("Frequency", fontsize=15)
    axh.set_xticks([])
    
    plt.show()

    


def MD_score_file(FILE):
    FH,G=open(FILE, 'r'),{}
    for line in FH:
        if "Binned" in line:
            break
        elif line[0]!="#":
            line_array=line.strip("\n").split("\t")
            motif, NS, scores = line_array[0], map(float, line_array[1].split(",")),map(float, line_array[2].split(","))    
            G[motif]=(NS, scores)
    FH.close()
    return G

def motif_displacements_raw(FILE):
    FH,G=open(FILE, 'r'),{}
    collect=False
    for line in FH:
        if collect and ("Empiracle" in line or "Multi" in line):
            break
        elif not collect and "Binned" in line:
            collect=True
        elif collect:
            line_array=line.strip("\n").split("\t")
            motif, vals = line_array[0], map(float, line_array[1].split(","))
            G[motif]=vals
    FH.close()
    return G

def get_obs_and_rand(FILE):
    
    FH,G=open(FILE, 'r'),{}
    for line in FH:
        if "Binned" in line:
            break
        elif line[0]!="#":
            line_array=line.strip("\n").split("\t")
            motif     = line_array[0]
            obs, rand = map(float, line_array[2].split(",")), map(float, line_array[-1].split(","))      
            G[motif]  = obs,rand
    return G

def load_motif_displacements_background_MD(FILE):
    FH,G=open(FILE, 'r'),{}
    collect=False
    for line in FH:
        if not collect and "Motif Id" in line[:20]:
            collect=True
        elif collect:
            line_array=line.strip("\n").split("\t")
            motif, vals = line_array[0], map(float, line_array[1].split(","))
            G[motif]=vals
    FH.close()
    return G

def gc_content_from_simulated_DB(FILE):
    FH = open(FILE, 'r')
    collect, L = False, list()
    for line in FH:
        line_array = line.strip("\n").split("\t")
        if collect and "Estim" in line:
            break
        elif "Estim" in line:
            collect = True
        elif collect and len(line_array)==2:
            x = [float(i) for i in line_array[1].split(",")]
            L.append(x)
    return L
def PSSM_and_simulated_draws(FILE):
    motif,FH  = None,open(FILE, 'r')
    P,D       = {},{}
    for line in FH:
        if line[0]==">":
            motif       = line[1:].strip("\n")
            P[motif]    = list()
        elif motif is not None and line[0]=="~":
            if motif not in D:
                D[motif]    = list()
            D[motif].append([float(x) for x in line.split("|")[1].split(",")])
        elif motif and line[0]!="#":
            P[motif].append([float(x) for x in line.strip("\n").split(",")])
    for motif in D:
        D[motif]=D[motif][0]
    return P,D
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

