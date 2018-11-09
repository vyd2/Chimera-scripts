from __future__ import print_function
##%matplotlib inline
import matplotlib.pyplot as plt
import itertools
import mdtraj as md
import mdtraj.testing
from matplotlib.pyplot import cm 
import numpy as np
import pandas as pd

def get_salt_pairs(traj): 
	'''
	traj = trajectory
	'''

    traj_top = traj.topology

    # get indexes of basic atoms involved in salt bridges 
    # for ARG: changed NH to HH*
    # for LYS: changed NZ to HZ*
    list_basicR  = [i for i in range(0,traj.n_atoms) 
                    if any(name in traj_top.atom(i).name for name in ('NE', 'NH','HH','HE')) 
                    and ('ARG' in traj_top.atom(i).residue.name)]
    list_basicK  = [i for i in range(0,traj.n_atoms)
                    if any(name in traj_top.atom(i).name for name in ('HZ', 'NZ')) 
                    and ('LYS' in traj_top.atom(i).residue.name)]

    # get indexes of acidic atoms involved in salt bridges 
    list_basic = list_basicR + list_basicK 
    
    # get indexes of acidic atoms involved in salt bridges 
    #exec("list_acid  = [i for i in range(0,ff%s1.n_atoms) if ('OP' in traj_top.atom(i).name)]" %ff)
    list_acidE  = [i for i in range(0,traj.n_atoms) 
                if ('OE' in traj_top.atom(i).name) and ('GLU' in traj_top.atom(i).residue.name)]
    list_acidD  = [i for i in range(0,traj.n_atoms) 
                if ('OD' in traj_top.atom(i).name) and ('ASP' in traj_top.atom(i).residue.name)]
    
    list_acid = list_acidE + list_acidD 
    
    
    pairs = list(itertools.product(list_acid, list_basic))
    print(pairs)

    return traj_top,pairs

def pair_distances(sframe,eframe,pairs)
	'''
	sframe = starting frame
	eframe = ending frame
	'''

	conts_perframe=[]

	for frame in range(sframe,eframe):
		conts = md.compute_distances(traj[frame],pairs)
		conts_perframe.append(conts)

	return 