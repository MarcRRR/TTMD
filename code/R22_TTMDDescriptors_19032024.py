import os
import sys
import MDAnalysis as mda
import oddt
from oddt import fingerprints
import numpy as np
import sklearn.metrics.pairwise as pair
import subprocess
from GET_RES_from_top import GET_RES_from_top

def desc(ligname, pdb):
	scores = []
	print('....', pdb)
	u = mda.Universe(pdb)
	prot = u.select_atoms('protein')
	protfile = "protein.pdb"
	with mda.Writer(protfile, prot.n_atoms) as W:
		W.write(prot)
	if params["HMR"] == True:
		topfile = "OneLig.NoWat.{}_HMR_vdWRep.top".format(ligname)
	else:
		topfile = "OneLig.NoWat.{}_vdWRep.top".format(ligname)

	num_res_prot,num_lig,name_res_lig,pos_LIG,pos_lig_NoWat,num_res_COFTOT = GET_RES_from_top(topfile)

	lig = u.select_atoms('resname {}'.format(name_res_lig))
	ligfile = "ligand.pdb"
	with mda.Writer(ligfile, lig.n_atoms) as W:
		W.write(lig)

	mols_p = next(oddt.toolkit.readfile('pdb', protfile))
	mols_l = next(oddt.toolkit.readfile('pdb', ligfile))


	ref = fingerprints.InteractionFingerprint(mols_l, mols_p, strict=True)
	l_plif_temp = []
	l_plif_temp.append(ref)
	print('0 ',l_plif_temp)
	sims = []

	u = mda.Universe(topfile, "{}.nc".format(pdb[:-10]))
	
	global lu
	lu = len(u.trajectory)
	print('lenght =',lu)


	for i, ts in enumerate(u.trajectory):
		scores.append([u, i])

	u_adapted = []
	for i in range(0, len(scores), len(scores)//150):
		u_adapted.append(scores[i])
	lu_step = len(u_adapted)
	print('lenght adapted =',lu_step)
	count = 0

	for u_n,i_n in u_adapted:
		
		u_n.trajectory[i_n]		
		prot = u_n.select_atoms('protein')
		protfile = "protein_new_{}.pdb".format(count)
		with mda.Writer(protfile, prot.n_atoms) as W:
			W.write(prot)

		lig = u_n.select_atoms('resname {}'.format(name_res_lig))
		ligfile = "ligand_new_{}.pdb".format(count)
		with mda.Writer(ligfile, lig.n_atoms) as W:
			W.write(lig)

		mols_p = next(oddt.toolkit.readfile('pdb', protfile))
		mols_l = next(oddt.toolkit.readfile('pdb', ligfile))
		fp = fingerprints.InteractionFingerprint(mols_l, mols_p, strict=True)
		l_plif_temp.append(fp)
		print('{} {}'.format(count,l_plif_temp))

		os.system("rm protein.pdb")
		os.system("rm ligand.pdb")

		mat = np.stack(l_plif_temp, axis=0)
		idx = np.argwhere(np.all(mat[..., :] == 0, axis=0))
		mat_dense = np.delete(mat, idx, axis=1)
		x = mat_dense[0].reshape(1,-1)
		y = mat_dense[1].reshape(1,-1)
		if y.size > 0:
			sim = round(float(pair.cosine_similarity(x, y)) * -1, 2)
		else:
			sim = 0
		sims.append(sim)
		print("similarity:", sim)
		l_plif_temp = []
		l_plif_temp.append(ref)
		count += 1

	return sims


ligs = []
DIR      = os.getcwd()
list_DIR = os.listdir(DIR)
for files in list_DIR :
        exist_dir = os.path.isdir(files)
        if exist_dir and files != "graphs" and files != "__pycache__":
                ligs.append(files)
ligs.sort()
print(ligs)

global params
params = eval(open("input_file", "r").read())

for lig in ligs:

	if os.path.exists("out_{}_sim.txt".format(lig)):
		os.system("rm out_{}_sim.txt".format(lig))

	vals = []

	f = open("out_{}_sim.txt".format(lig), "x")

	f.write("##### SIMILARITY COEFFICIENTS #####\n\n")

	os.chdir("{}".format(lig))
	dyn_num = int(subprocess.check_output("ls | grep 'dyn_*' | wc -l", shell=True))

	f.write("frames")
	for dyn in range(1, dyn_num+1):

		for i in range(0, 10000, 10000//150):
			f.write("\t{}".format(i))
		f.write("\n")
		f.write("\t##### DYN {} #####\n\n".format(dyn))
		os.chdir("dyn_{}/TTMD".format(dyn))

		files_to_open = subprocess.check_output("ls | grep _First.pdb", shell=True).decode("utf-8").strip().split()

		for pdb in files_to_open:
			f.write("{}\t".format(pdb[:-10]))
			lis = desc(lig,pdb)
			vals.append(lis)

			for i in range(len(lis)):
				f.write("{}\t".format(lis[i]))
			f.write("\n")
	
			f.write("\n")
			
		os.chdir("../..")
	f.close()
	os.chdir("..")	

