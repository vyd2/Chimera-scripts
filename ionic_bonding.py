from __future__ import print_function
#%matplotlib inline
import matplotlib.pyplot as plt
from itertools import compress
import mdtraj as md
import mdtraj.testing
import numpy as np
import pandas as pd


####################

def load_traj(parm,traj,sframe,eframe):
	'''
	Load the traj to get atom indices
	'''

	# load traj
	t = md.load(traj, top)

	# get topology for traj 
	traj_top = t.topology

	# make list pairs for distance eval
	list_basicH  = [i for i in range(0,t.n_atoms) 
                if any(name in traj_top.atom(i).name for name in ('NE2'))
                and ('HIS' in str(traj_top.atom(i).residue.name))]
	list_basicK  = [i for i in range(0,t.n_atoms)
	                if any(name in traj_top.atom(i).name for name in ('HZ', 'NZ')) 
	                and ('LYS' in traj_top.atom(i).residue.name)]
	list_basicR  = [i for i in range(0,t.n_atoms) 
	            if any(name in traj_top.atom(i).name for name in ('NE', 'NH')) 
	            and ('ARG' in traj_top.atom(i).residue.name)]

	# get indexes of acidic atoms involved in salt bridges 
	list_pos = list_basicR + list_basicK + list_basicH

	# get indexes of acidic atoms involved in salt bridges 
	list_acidE  = [i for i in range(0,t.n_atoms) 
	            if ('OE' in traj_top.atom(i).name) and ('GLU' in str(traj_top.atom(i).residue.name))]
	list_acidD  = [i for i in range(0,t.n_atoms) 
	            if ('OD' in traj_top.atom(i).name) and ('ASP' in traj_top.atom(i).residue.name)]
	list_neg = list_acidE + list_acidD 

	pairs = list(itertools.product(list_acid, list_basic))

	return traj_top,pairs 


def count_ibond(traj_top, pairs,list=None):
	'''
	Calculate distances b/c pairs of atoms 
	'''

	if list is None: 
		list = []
	
	conts = md.compute_distances(traj_top, pairs)
	list.append(conts)

	return list


def get_ionbds(traj_top,list,pairs):
	'''
	This is to get a final pandas dataframe of the frequency of pairs of ionic bonding atoms 
	'''

	resn = lambda x: traj_top.atom(x).residue
	big = pd.DataFrame()
	newrang = range(len(list))
	print(ff,len(list))

	#loop through length of each conts_perframe
	for i in newrang:
	    
	    # first: generate bool of T/False 
	    # second: mult 10 to get angstroms, bools based on if distance less than 4
	    # was 5.5 Angstroms cutoff
	    singarr = ((list[i]*10)<6)[0]
	    
	    # there's at least one True in the boolarr
	    if np.sum(singarr) > 1:
	        
	        # pair together pairs of residues and the arrays of singarr 
	        dtest = pd.DataFrame(list(compress(pairs,singarr)),columns=['donor','accept']) 
	        dtestn = dtest.applymap(resn) # apply resn lambda
	        dtestnn = dtestn.drop_duplicates(subset=['donor', 'accept'], keep='first') # remove duplicates
	        big = pd.concat([big,dtestnn])

	final_df = big.groupby(big.columns.tolist(),as_index=False).size()

	return final_df
