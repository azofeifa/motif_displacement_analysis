from math import exp,log,sqrt,pi
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import differential_MD as dMD
class uniform:
	def __init__(self, a,b, w):
		self.a 	= a
		self.b 	= b
		self.l 	= abs(b-a)
		self.w 	= w
		self.type 	= "uniform"
	def __str__(self):
		return "uniform (" + str(self.a) + "," + str(self.b) + "," + str(self.w)+ ")"

	def pdf(self, x):
		return int(self.a<=x<=self.b)*(self.w/self.l)

class normal:
	def __init__(self, mu,sigma, w):
		self.mu 	= mu
		self.b 	= sigma
		self.w 	= w
		self.type 	= "normal"
	def __str__(self):
		return "norm (" + str(self.mu) + "," + str(self.b) + "," + str(self.w)+ ")"

	def pdf(self, x):
		try:
			return (self.w / sqrt(2*pow(self.b,2)*pi) )*exp(-pow(x-self.mu,2)/(2*pow(self.b,2)))
		except:
			return 0.0


class laplace:
	def __init__(self, mu, b, w ):
		self.mu 	= float(mu)
		self.b 		= float(b)
		self.w 		= w
		self.type 	= "laplace"
	def __str__(self):
		return "laplace (" + str(self.mu) + "," + str(self.b) + "," + str(self.w) + ")"
	def pdf(self, x):
		try:
			return (self.w / (2.0*self.b))*exp(-abs(x-self.mu)/self.b)
		except:
			return 0.0
class EM:
	def __init__(self, k=1, T=300, ct=pow(10,-5)):
		self.k 	= k
		self.T 	= T
		self.ct = ct
		self.rvs=list()
		self.K 	= {}
		self.N 	= 0.0
	def run_all_three(self, X):
		for k in range(3):
			rvs, ll 	= self.fit(X,k)
			self.K[k] 	= (ll, rvs)
	def get_bic(self):
		val 	= min([ (-2*self.K[k][0] + log(self.N)*pow(k+1,1),k) for k in self.K])
		return val[1]
	def fit(self,X,k,norm=True):
		a,b 	= X[0,0], X[-1,0]
		self.N 	= sum(X[:,1])
		if  k == 1:
			if norm:
				CC 	= normal(0,50,0.5)
			else:
				CC 	= laplace(0,150,0.5)
		elif k == 2:
			if norm:
				CC 	= normal(-30,50,0.33),normal(30,50,0.33)
			else:
				CC 	= laplace(-30,150,0.33),laplace(30,150,0.33)


		if k == 0:
			rvs 	= [uniform(a,b,0.5)]
			return rvs, log(1.0/(b-a))*sum(X[:,1])
		elif k == 1:
			rvs 	= [uniform(a,b,0.5) , CC]
		elif k == 2:
			rvs 	= [uniform(a,b,0.33) , CC[0], CC[1]]
		t 			= 0
		self.rvs 	= rvs
		prevll 		= -np.inf
		ll 	= 0.0
		while t < self.T and k > 0:
			EX 	= np.zeros((k+1,))
			mus 	= np.zeros((k,))
			bs 	= np.zeros((k,))
			ll 	= 0.0
			for i in range(X.shape[0]):
				norms 	= np.zeros((k+1,))
				for j,rv in enumerate(rvs):
					norms[j] 	= rv.pdf(X[i,0])
				if sum(norms) > 0:

					ll+=log(sum(norms))*X[i,1]
					norms[:]/=sum(norms)
					for j,rv in enumerate(rvs):
						if rv.type != "uniform":
						 	if k > 1:
								mus[j-1]+=norms[j]*X[i,0]*X[i,1]
							if rv.type=="laplace":
								bs[j-1]+=norms[j]*abs(X[i,0]-rv.mu)*X[i,1]
							else:
								bs[j-1]+=norms[j]*pow(X[i,0]-rv.mu,2)*X[i,1]

					for j in range(k+1):
						EX[j]+=norms[j]*X[i,1]
			for j,rv in enumerate(rvs):
				rv.w 	= (EX[j] ) / (sum(EX) )
				if rv.type!="uniform":# and k ==1:
					rv.b 	= bs[j-1] / EX[j]
					if norm:
						rv.b=sqrt(b)
			if k > 1:
				b 	= ( (bs[0]+1) / (EX[1]+1)  ) + ( (bs[1]+1) / (EX[2]+1))
				b 	/= 2.0
				if norm:
					b 	= sqrt(b)
				v 	= abs(mus[0]/EX[1]) + abs(mus[1]/EX[2])
				v 	/= 2.0
				rvs[1].mu 	=mus[0]/EX[1] #mus[0]/EX[1]
				rvs[2].mu 	=mus[1]/EX[2] #mus[1]/EX[2]
				rvs[1].b, rvs[2].b 	= b,b

				W 	= (rvs[1].w + rvs[2].w )/2.0
				rvs[1].w 	= W
				rvs[2].w 	= W
				S 			= sum([rv.w for rv in rvs])
				for rv in rvs:
					rv.w/=S 	
			if abs(ll - prevll) <self.ct:
				break

			prevll 	= ll
			t+=1
		self.rvs 	= rvs
		return rvs,ll
	def get_best(self,Y,norm=True):
		X 		= np.array([np.linspace(-1500,1500,len(Y) ),Y]).T

		n,BIC 			= np.sum(X[:,1]), np.inf
		arg_rvs,arg_ll = list(), -np.inf
		disp,weight, b = 0,0,0 
		for k in range(3):
			rvs,ll 	= self.fit(X,k,norm=False)
			current 	= -2*ll + log(n+1)*(1+5)*k
			if current < BIC:
				BIC 	= current
				arg_rvs,arg_ll 	= rvs, ll
				if k > 0:
					disp 	= rvs[1].mu
					weight= rvs[1].w*k
					b 	   = rvs[1].b
		return disp,weight,b

	def draw(self,X):
		ax 	= plt.gca()
		ax.hist(X[:,0], weights=X[:,1],normed=1.0,edgecolor="white",bins=X.shape[0],alpha=0.75)
		xs 	= np.linspace(X[0,0],X[-1,0], 200)
		ax.plot(xs, map(self.pdf, xs))
		plt.show()
	def pdf(self, x):
		return sum([rv.pdf(x) for rv in self.rvs])
	def density(self, x,mu,weight, b,l,norm=True):
		if abs(mu) < 1:
			if norm:
				return ((1.0-weight)/l) + (weight / sqrt(2*pow(b,2)*pi) )*exp(-pow(x-mu,2)/(2*pow(b,2)))
			return ((1.0-weight)/l) + ( weight)*exp(-abs(x-mu)/b)/(2.0*b)
		else:
			if norm:
				return ((1.0-weight)/l) + (0.5*weight / sqrt(2*pow(b,2)*pi) )*exp(-pow(x-mu,2)/(2*pow(b,2)))  + (0.5*weight / sqrt(2*pow(b,2)*pi) )*exp(-pow(mu+x,2)/(2*pow(b,2)))  
			return ((1.0-weight)/l) + ( weight*0.5)*exp(-abs(x-mu)/b)/(2.0*b)+ ( weight*0.5)*exp(-abs(x+mu)/b)/(2.0*b)

	def get_ll(self, X):

		return sum([log(self.pdf(X[i,0]))*X[i,1]  for i in range(X.shape[0])])



