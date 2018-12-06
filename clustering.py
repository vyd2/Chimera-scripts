#%matplotlib inline
from __future__ import print_function
import mdtraj as md
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import numpy as np
from scipy.cluster.vq import kmeans
from scipy.spatial.distance import cdist,pdist
from sklearn import datasets
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn import mixture 
import scipy.stats
from itertools import compress
import pandas as pd


sns.set(font_scale=2) 
sns.set_style("ticks")
sns.set_palette("husl")

#############################################

# ALL
import time 

s_time = time.time()

def get_dihedrals(top,traj):
	'''
	Get dihedral angles
	'''

	phi_list=[]
	psi_list=[]

	# load traj
	t = md.load(traj,top)

	# align frames to the first frame
	t.superpose(t, 0)

	exec("ind,phi = md.compute_phi(t%s)" %i)
	exec("ind2,psi = md.compute_psi(t%s)" %i)

	phi_list.append(phi)
	phi_list.append(phi)

	X = [np.concatenate(x) for index, x in enumerate(zip(np.concatenate(phio),np.concatenate(psio))) 
	if index in og_list]

	return phi,psi,X


def PCA_reduction(X):
	'''
	X = merged phi and psi dataset 
	'''

	# PCA
	pca = PCA(0.99,whiten=True)
	X2 = pca.fit(X).transform(X)
	print(dihedrals, X2.shape)

	
	# plot AICs and BIC
	n_components = np.arange(1, 300, 5)
	models = [mixture.GaussianMixture(n, covariance_type='full', random_state=0)
	          for n in n_components]
	aics = [model.fit(X2).aic(X2) for model in models]
	bics = [model.fit(X2).bic(X2) for model in models]
	minbic= dict(zip(n_components,bics))
	minaic= dict(zip(n_components,aics))
	bic_min = min(minbic,key=minbic.get)
	aic_min = min(minaic,key=minaic.get)
	print(aic_min,bic_min)

	# Plot of find ideal cluster size for later GMM
	plt.plot(n_components, aics,'d-',label='AIC')
	plt.plot(n_components, bics,'d-',label='BIC')
	plt.legend()

	# check if converged
	gmm = mixture.GaussianMixture(bic_min, covariance_type='full', random_state=0)
	gmm.fit(X2)
	cIdx = gmm.predict(X2)
	print(gmm.converged_)

	return cIdx,gmm,X2,bic_min

def get_cluster_centers(cIdx,gmm,X2,bic_min):

	ids =[]
	rep_c =[]
	avgdist =[]
	perc_fr =[]
	
	# representative of clusters
	centers = np.empty(shape=(gmm.n_components, X2.shape[1]))

	for i in range(gmm.n_components):
	    density = scipy.stats.multivariate_normal(cov=gmm.covariances_[i], mean=gmm.means_[i
	    centers[i, :] = X2[np.argmax(density)]

	# calc. euc. dist. b/n cent and X2 array, each column represents dist. calc. for a singl
	D_k = [cdist(X2, cent, 'euclidean') for cent in [centers]] 

	# generate list of indices of minimum values of each ROW, for each euc. dist. for K=1..1
	# analogous distances that correspond to cIdx
	dist = [np.min(D,axis=1) for D in D_k]

	kIdx=bic_min # cluster number 
	cluster = []

	for i in range(kIdx):

	    # get list of ones of bool values
	    ind = (cIdx==i)          

	    # to get matching indices into cluster
	    cluster.append(list(compress(range(len(ind)), ind)))  

	     # x=index, and corresponding index, each lst_c is new, with dist from cen
	    lst_c = { x : dist[0][x] for x in cluster[i]}     
	  
	  	# get min. val from dictionary lst_c, and gives the index of the value in the list 
	    rep = min(lst_c, key=lst_c.get)     

	    # get avg. dist. in each cluster
	    Avgdist = np.mean(lst_c.values())   

	    perc = len(cluster[i])/100000
	    
	    rep_c.append(all_perc[rep]+1)
	    avgdist.append(Avgdist)
	    perc_fr.append(perc)

	    print('Cluster %s:' %i,og_list[rep]+1,Avgdist,perc) # print out info

	result = pd.DataFrame({'Rep_fr':rep_c,'Avgdist':avgdist_%s,'percent_fr':perc_fr})
	result.sort_values(by='percent_fr',inplace = True,ascending=False)
	print(result.iloc[0:10])

	return bics, cluster
	
	
	




