class uni_frame:
    def __init__(self, 
                    path_to_motif_bed_files="",
                    path_to_chrom_sizes=""):
        self.pmbf   = path_to_motif_bed_files
        self.pcs    = path_to_chrom_sizes
        self.c_sizes= {}
        self.M      = {}
        self.EXPS    = {}
        self.load_pmbf()
    def load_pmbf(self):
        S   = 0
        for line in open("../data_azofeifa_2016/hg19_chromosome_sizes.tsv"):
            chrom,size  = line.strip("\n").split('\t')
            size        = float(size)
            self.c_sizes[chrom]=S
            S+=size
        print "loading motif bed files",

        for f in os.listdir(self.pmbf):
            m   = list()
            for line in open(self.pmbf+f):
                try:
                    chrom,start,stop    = line.split("\t")[:3]
                    x = 0.5*(float(start) + float(stop))
                    m.append(x+self.c_sizes[chrom])
                except:

                    pass
            m.sort()
            self.M[f]   = m
        print "done"
    def load_tfit_file(self, FILE,exp):
        if exp not in self.EXPS:
            self.EXPS[exp]  = list()
        m   = list()
        for i,line in enumerate(open(FILE,"r")):
            '''
                is this a comma separated one?
            '''
            try:
                line_array=line.strip("\n").split(",")
            except:
                '''
                    is this a tab separated one (bed)?
                '''
                try:
                    line_array=line.strip("\n").split("\t")
                except:
                    print "cant understand this file at line:",i
                    return -1
                if len(line_array)!=4:
                    print "cant understand this file at line:",i
                    return -1
            chrom,start,stop, tss   = line_array

            if chrom in self.c_sizes and tss=="0":
                x = 0.5*(float(start) + float(stop))
                m.append(x+self.c_sizes[chrom])
        m.sort()
        self.EXPS[exp].append(m)
    def unqiue(self,A,B,W,overlap=False):
        j,N,U     = 0, len(B),list()
        for a in A:
            while j < N and B[j]+W < a-W:
                j+=1
            if j < N and B[j]-W < a+W and overlap:
                U.append(0.5*(a+B[j]))
            elif j < N and B[j]-W > a+W and not overlap:
                U.append(a)
        return U
    def two_tail(self, x):
        if x < 0.5:
            return max(x*2,pow(10,-40))
        return max((1.0-x)*2,pow(10,-40))
    def enrichment(self,A,B,W=50):
        n,m             = len(self.EXPS[A]),len(self.EXPS[B])
        if n <= m:
            comparisons     =  [(i,j) for i in range(n) for j in range(i,m)  ]
        else:
            comparisons     =  [(i,j) for j in range(m) for i in range(j,n)  ]
        mN      = len(self.M)
        motifs  = self.M.keys()
        FCS,PVS = [np.zeros(mN*2) for k in range(2)]
        for i,j in comparisons:
            '''
                unique to A
            '''
            UA  = self.unqiue(self.EXPS[A][i],self.EXPS[B][j], 1500)
            '''
                unqiue to B
            '''
            UB  = self.unqiue(self.EXPS[B][j],self.EXPS[A][i], 1500)
            '''
                overlap to AB
            '''
            O   = self.unqiue(self.EXPS[B][j],self.EXPS[A][i], 1500,overlap=True)
            '''
                iterate over motif models...this will take some time I imagine
            '''
            unqiue_to_A_N, unqiue_to_B_N    = len(UA),len(UB)
            all_of_A_N,all_of_B_N           = len(self.EXPS[A][i]),len(self.EXPS[B][j])
            N           = unqiue_to_B_N + len(O) 
            for mi,m in enumerate(motifs):
                '''
                    how many hits are there in general?
                '''
                all_A_hits      = self.unqiue(self.EXPS[A][i],self.M[m],W,overlap=True)
                all_B_hits      = self.unqiue(self.EXPS[B][j],self.M[m],W,overlap=True)

                AN      = self.unqiue(UA,self.M[m],W,overlap=True)
                BN      = self.unqiue(UB,self.M[m],W,overlap=True)

                hg1     = ss.hypergeom(all_of_A_N, len(all_A_hits), unqiue_to_A_N )
                hg2     = ss.hypergeom(all_of_B_N, len(all_B_hits), unqiue_to_B_N )

                print m.split("_")[1], all_of_B_N, len(all_B_hits), unqiue_to_B_N , len(BN), hg2.mean()

                pv1     = self.two_tail(hg1.cdf(len(AN)))
                pv2     = self.two_tail(hg2.cdf(len(BN)))
                
                FCS[mi]   +=len(AN)/(hg1.mean()+1)
                FCS[mi+mN]+=len(BN)/(hg2.mean()+1)
                PVS[mi]   +=-2*math.log(pv1)
                PVS[mi+mN]+=-2*math.log(pv2)
        FCS/=len(comparisons)
        direct  = np.array([0 for i in range(mN)]+[1 for i in range(mN)])

        '''
            fishers method for pvalue combination
        '''

        pvals   = 1.0-ss.chi2(len(comparisons)*2).cdf(PVS)
        '''
            adjusting for multiple hypothesis testing
        '''
        reject, pvals_corrected,alphacSidak,alphacBonf  =  multitest.multipletests(pvals)
        '''
            sort by the most significant
        '''
        idx                     = sorted(range(len(pvals)), key=lambda k: FCS[k],reverse=True)
        motifs                  = np.array(motifs+motifs)[idx]
        pvals,pvals_corrected   = pvals[idx],pvals_corrected[idx]
        direct,FCS              = direct[idx],FCS[idx]

        '''
            through this into a pandas dataframe
        '''
        df      = pd.DataFrame(columns=("motif", "treatment_unique","fold_over_expectation", "pvalue", "adjpvalue"  ) )

        for i,m in enumerate(motifs):
            df.loc[i]   = m, direct[i], FCS[i], pvals[i],pvals_corrected[i]
        '''
            write out csv file
        '''
        df.to_csv("motif_enrichment_results.csv", index=False)
