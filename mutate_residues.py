# Purpose of this script is to use Chimera to automate mutating residues 
# Open chimera, File > Open, Select this python script 

import os
from chimera import runCommand as rc
from chimera import replyobj

# change to folder with data files
os.chdir("/Users/vyduong/cloudv/WW/YP-WW")

file = ["YP.pdb"] # name of your pdb file 

#lst = [4,5,7,9,20,23,27,33]
lst = ['4','5'] # note these need to be a string not floats 

for i in lst: 
	for fn in file: 
		replyobj.status("Processing " + fn) # show what file we're working on
		rc("open " + fn) # need these plus signs to insert variable
		rc("swapaa ala :"+i+".a") #swapaa command 'swaps' residue e.g. 4 into an alanine 
		rc("write #0 ~/cloudv/WW/YP-WW/YP_"+i+"A.pdb") # write out pdbfile to location
		rc("close all") #close session

# rc(close all)