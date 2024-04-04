import os
import sys
import subprocess
from math import sqrt
from GET_RES_from_top import GET_RES_from_top


def check_files(lig, dyn_num):
	if not os.path.exists("TTMD"):
		print("no TTMD directory has been found.")
		exit()
	
	os.chdir("TTMD")
	needed_files = ["in", "rst", "gamd_rst"]
	for suffix in needed_files:
		numfiles = int(subprocess.check_output("ls {}*ttmd.{} | wc -l".format(lig,suffix), shell=True))
		if numfiles == 0:
			print("no {} file has been found".format(suffix))
	numfiles = int(subprocess.check_output("ls {}*vdWRep.top | wc -l".format(lig), shell=True))
	if numfiles == 0:
		print("no {} file has been found".format(suffix))
	os.chdir("..")
	return 

def cpptraj_rmWat(lig, clean=False):

	name_traj = "{}_1_ttmd.nc".format(lig)
	parm_wat = "[Wat]"


	file_cpptraj_inp = "Ainptrj_nowat"
	file_cpptraj_out = "Ainptrj_nowat.out"
	file_cpptraj = open(file_cpptraj_inp, "w")

	file_cpptraj.write('cpptraj > {} << EOF '.format(file_cpptraj_out) + ter)
	file_cpptraj.write(' parm      {}   {}   '.format(top_file,parm_wat)+ ter)
	file_cpptraj.write(' '+ ter)

	file_cpptraj.write(' trajin {}  1  last  1 parm {} '.format(name_traj,parm_wat)+ter)
	file_cpptraj.write(' '+ ter)
	file_cpptraj.write(' center  :1-{}  origin '.format(num_res_PROT)+ ter) # number of resids of protein
	file_cpptraj.write(' image origin '+ ter)
	file_cpptraj.write(' strip :WAT,Na+,Cl-  outprefix  NoWat  '+ ter)
	file_cpptraj.write(' rms first out RMS_first.dat :{}-{}@CA '.format("1", num_res_PROT)+ ter) # ini_res_sup,ifi_res_sup

	outfile = "RMSD_FIRST_NoWat.nc"
	file_cpptraj.write(' trajout {} '.format(outfile)+ ter)
	file_cpptraj.write(' run '+ ter)
	file_cpptraj.write(' '+ ter)
	file_cpptraj.write('EOF'+ ter)
	file_cpptraj.close()

	run_cpptraj = "./" + file_cpptraj_inp
	os.system('chmod u+x {} '.format(run_cpptraj))
	os.system(run_cpptraj)

	if clean :
		os.system  ('rm {} '.format(file_cpptraj_inp) )
		os.system  ('rm {} '.format(file_cpptraj_out) )

	return outfile



def cpptraj_LigByLig(traj_file, clean=False):

	parm_NoWat = "[NoWat]"
	count_ligand = num_res_PROT + 1
	top_OneLigNoWat = "OneLig.NoWat." + top_file
	lig_trajs = []

	if num_lig != 1 :
		for ligand in range ( num_lig ) :

			name_nc_lig = "lig_" + str(count_ligand) + ".nc"

			lig_trajs.append(name_nc_lig)

			file_cpptraj_inp  = "Binptrj_ligands"
			file_cpptraj_out  = "Binptrj_ligands.out"
			file_cpptraj      = open (file_cpptraj_inp,"w")
	
			file_cpptraj.write('cpptraj > {} << EOF'.format(file_cpptraj_out)+ter)
			file_cpptraj.write(' parm NoWat.{} {}'.format(top_file, parm_NoWat)+ter)
			file_cpptraj.write(' '+ter)


			file_cpptraj.write(' trajin {} 1 last 1 parm {}'.format(traj_file, parm_NoWat)+ter) #missing arg
	
			if ligand == 0 :
				file_cpptraj.write(' strip !:1-{},{} outprefix OneLig'.format(num_res_PROT, count_ligand)+ter)

			else:

				file_cpptraj.write(' strip !:1-{},{} outprefix OneLig'.format(num_res_PROT, count_ligand)+ter)

			file_cpptraj.write(' trajout {} {}'.format(name_nc_lig, parm_NoWat)+ter)

			file_cpptraj.write(' run '+ter)
			file_cpptraj.write(' '+ter)
			file_cpptraj.write('EOF'+ter)

			file_cpptraj.close()

			run_cpptraj = "./" + file_cpptraj_inp
			os.system  ('chmod u+x {} '.format(run_cpptraj) )
			os.system (run_cpptraj)

			count_ligand += 1

	if clean :
		os.system  ('rm {} '.format(file_cpptraj_inp) )
		os.system  ('rm {} '.format(file_cpptraj_out) )

	return lig_trajs


def cpptraj_OutPDBFirst(traj_file, clean=False):

	parm_OneLigNoWat = "[OneLigNoWat]"
	count_ligand = num_res_PROT + 1
	top_OneLigNoWat = "OneLig.NoWat." + top_file
	count = 1
	lig_pdbs = []	
	for i in traj_file:

		name_nc_lig = "lig_{}_First.pdb".format(num_res_PROT+count)
		lig_pdbs.append(name_nc_lig)

		file_cpptraj_inp  = "Cinptrj_pdbFirst"
		file_cpptraj_out  = "Cinptrj_pdbFirst.out"
		file_cpptraj      = open (file_cpptraj_inp,"w")
	
		file_cpptraj.write(' cpptraj > {} << EOF'.format(file_cpptraj_out)+ter)
		file_cpptraj.write(' parm OneLig.NoWat.{} {}'.format(top_file, parm_OneLigNoWat)+ter)
		file_cpptraj.write(' '+ter)


		file_cpptraj.write(' trajin {} 1 1 1 parm {}'.format(i, parm_OneLigNoWat)+ter)
	

		file_cpptraj.write(' trajout {} {}'.format(name_nc_lig, parm_OneLigNoWat)+ter)

		file_cpptraj.write(' run '+ter)
		file_cpptraj.write(' '+ter)
		file_cpptraj.write('EOF'+ter)

		file_cpptraj.close()

		run_cpptraj = "./" + file_cpptraj_inp
		os.system  ('chmod u+x {} '.format(run_cpptraj) )
		os.system (run_cpptraj)
		count += 1

	if clean :
		os.system  ('rm {} '.format(file_cpptraj_inp) )
		os.system  ('rm {} '.format(file_cpptraj_out) )

	return lig_pdbs


def read_prot (pdb):
	ABC = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
	abc = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
	protein = ['SER','THR','GLN','ASN','TYR','CYS','CYX','CYM','GLY',  \
                   'ALA','VAL','LEU','ILE','MET','PRO','PHE','TRP',        \
                   'GLU','GLH','ASP','ASH','LYS','ARG','HIE','HID','HIP',  \
                   'PHE','TYR','TRP','ACE','NME','NHE','GLP','NHE' ]
	membrane  = ['PA','PC','OL']
	cofactor  = ['GTP','ATP']
	ionscofac = ['MG','CA']
# print ( '    ... Reading     : {}  '.format(pdb) ) 
	file_pdb  = open (pdb)
	num_lines = 0
	for lines in file_pdb :
		num_lines += 1
# print ( "    ... There are = ",num_lines," Lines " )
	file_pdb  = open (pdb)

	num_atoms     = 0
	num_atoms_lig = 0



	at_CA    = []

	for i in range ( num_lines )  :

		line = file_pdb.readline().strip()
		if 'ATOM' in line  :
			values=line.split()
			res_name_dum   = values [3]
			if res_name_dum in protein :
				num_atoms += 1
				is_CA  = values [2]
				if is_CA == 'CA' :
					at_CA.append(num_atoms)
					at_name.append(at_CA)
				res_name.append (values [3])
				chain_prot = values [4]
				if (chain_prot in ABC) or (chain_prot in abc) :
					res_num.append (int (values [5]))
					x_prot.append (float(values [6]))
					y_prot.append (float(values [7]))
					z_prot.append (float(values [8]))
				else:
					res_num.append (int (values [4]))
					x_prot.append (float(values [5]))
					y_prot.append (float(values [6]))
					z_prot.append (float(values [7]))
			else :
				values=line.split()
				print ( values )
				atom_name_lig.append (values [2])
				resi_name_lig   = values [3]
				chain_lig = values [4]
				if (chain_lig in ABC) or (chain_lig in abc) :
					x_lig.append (float(values [6]))
					y_lig.append (float(values [7]))
					z_lig.append (float(values [8]))
				else :
					x_lig.append (float(values [5]))
					y_lig.append (float(values [6]))
					z_lig.append (float(values [7]))
				num_atoms_lig += 1

	return num_atoms, num_atoms_lig, resi_name_lig, at_CA




def calc_BindingSite (lig_num,num_at_prot,num_at_lig,dis_lig_prot_min,at_CA) :

	lig_BS = {}
	xm = ym = zm = 0
	for at in range(num_at_lig):
		xm += x_lig[at]
		ym += y_lig[at]
		zm += z_lig[at]
	x_gc = xm / num_at_lig
	y_gc = ym / num_at_lig
	z_gc = zm / num_at_lig

	dis_intra_MAX = 0
	for atA in range(num_at_lig-1):
		xA = x_lig[atA]
		yA = y_lig[atA]
		zA = z_lig[atA]
		for atB in range(atA+1,num_at_lig):
			xB = x_lig[atB]
			yB = y_lig[atB]
			zB = z_lig[atB]
			disAB = sqrt( (xA-xB)*(xA-xB) + (yA-yB)*(yA-yB) + (zA-zB)*(zA-zB) )
			if disAB > dis_intra_MAX :
				dis_intra_MAX = disAB
	print ( "    ... Ligand GC   {:7.3f} {:7.3f} {:7.3f} Max Intra-Dist {:7.3f} ".format(x_gc,y_gc,z_gc,dis_intra_MAX))

	res_inBS = []
	dis_LIM = dis_intra_MAX + dis_lig_prot_min
	for atP in range(num_at_prot):
		xP = x_prot [atP]
		yP = y_prot [atP]
		zP = z_prot [atP]
		res_BS = res_num [atP]
		if res_BS in res_inBS :
			continue
		disPL = sqrt( (xP-x_gc)*(xP-x_gc) + (yP-y_gc)*(yP-y_gc) + (zP-z_gc)*(zP-z_gc) )
		if disPL < dis_LIM :
			for atL in range(num_at_lig):
				xL = x_lig[atL]
				yL = y_lig[atL]
				zL = z_lig[atL]
				disPL = sqrt( (xP-xL)*(xP-xL) + (yP-yL)*(yP-yL) + (zP-zL)*(zP-zL) )
				if disPL < dis_lig_prot_min :
					res_inBS.append(res_BS)
					break

	num_resBS = len(res_inBS)

	print ( "    ... There are  {:4} Protein Residues in this Binding Site ".format(num_resBS))

	if num_resBS != 0:

		fileBS.write( str(num_resBS) + '\n' )

	llista = ""
	line_wrt = ''
	for res in range (num_resBS) :
		line_wrt += str(res_inBS[res]) + ' '
		llista  += (str(res_inBS[res])) + ","

	print("lig_" + str(lig_num))
	print(llista[:-1])
	lig_BS["lig_" + str(lig_num)] = llista[:-1].split(",")
	print("ARA")

	if num_resBS != 0:
		fileBS.write( line_wrt + '\n' )

	print ( "    ... RESIDUES : " )
	print ( line_wrt )
	print ( "... ... " )

	for at in range (num_resBS) :
		res_pos = res_inBS [at]
		CA_pos  = at_CA    [res_pos-1]
		x_CA    = x_prot   [CA_pos-1]
		y_CA    = y_prot   [CA_pos-1]
		z_CA    = z_prot   [CA_pos-1]
		print ( "    ... CA({}) :  {:7.3f} {:7.3f} {:7.3f} ".format(res_pos,x_CA,y_CA,z_CA))

	return num_resBS,res_inBS



def binding_sites_ALL(lig, dyn_num):

	if os.path.exists("{}_BS_all.txt".format(lig)):
		os.system("rm {}_BS_all.txt".format(lig))

	os.system("touch {}_BS_all.txt".format(lig))


	with open("{}_BS_all.txt".format(lig), "w") as aBS:

		count = 1
		bss = []
		BS_mat = [[False for i in range(dyn_num)]]

		for dyn in range(1, dyn_num+1):
			os.chdir("{}/dyn_{}/TTMD".format(lig, dyn))
			numlines = int(subprocess.check_output("cat reactive_trajs_First.txt | wc -l", shell=True).decode("utf-8"))

			if numlines != 0:
				file_to_open = subprocess.check_output("ls | grep *_BS.dat", shell=True).decode("utf-8").strip()
				BS_file = open(file_to_open, "r")

				for line in range(0, numlines, 2):
					num_prot = BS_file.readline()			
					resids = BS_file.readline()
					nums = set(resids.strip().split())

					if len(bss) == 0:
						aBS.write('BS{}\n'.format(count)+num_prot+resids)
						bss.append(nums)
						count += 1
						BS_mat[0][dyn-1] = True
						
					elif len(nums) != 0:
						is_in_list = False

						for bs in range(len(bss)):
							inter = len(nums & bss[bs])
							union = len(nums | bss[bs])

							if union == 0 or inter / union > 0.25:
								is_in_list = True
								BS_mat[bs][dyn-1] = True
								break


						if not is_in_list:
							aBS.write('BS{}\n'.format(count)+num_prot+resids)
							bss.append(nums)
							BS_mat.append([False for i in range(dyn_num)])
							BS_mat[-1][dyn-1] = True
							count += 1
				BS_file.close()
			os.chdir("../../..")
		app_bs = {}
		numbs = 1
		fn = "{}_BS_appearance.txt".format(lig)
		if os.path.exists(fn):
			os.system("rm "+fn)


		with open(fn, "x") as f:
			for i in BS_mat:
				f.write("BS{} ".format(numbs) + str(sum(i)) + "\n")
				numbs += 1


	bscount = 1
	print("\tdyn_1\tdyn_2")
	for i in BS_mat:
		print("BS{}\t".format(bscount), end='')
		for j in i:
			print(j, end = "\t")
		bscount += 1
		print()
	return



##### MAIN PROGRAM #####

global ter, top_file, num_res_prot, num_lig, name_res_lig, pos_LIG, pos_lig_NoWat, num_res_COFTOT, num_res_PROT
global  x_prot, y_prot, z_prot, at_name, res_name, res_num
global  x_lig, y_lig, z_lig, atom_name_lig, resi_name_lig
global fileBS, countlig

ter = "\n"



ligs = []
DIR      = os.getcwd()
list_DIR = os.listdir(DIR)
for files in list_DIR :
        exist_dir = os.path.isdir(files)
        if exist_dir :
                ligs.append(files)
ligs.sort()
ligs = ligs[1:]


for lig in ligs:

	os.chdir(lig)
	print('\n :::WORKING ON {}:::\n'.format(lig))
	dyn_num = int(subprocess.check_output("ls | grep 'dyn*' | wc -l", shell=True))

	for i in range(1, dyn_num+1):	

		os.chdir("dyn_{}".format(i))
		print('\t:::WORKING ON DYN_{}:::\n'.format(i))
		check_files(lig, i)
		

		os.chdir("TTMD")
		if i == 1:
			top_file = subprocess.check_output("ls *top | grep -v NoWat",shell=True).decode("utf-8").strip()


		print("\t... Reading top file")
		num_res_prot,num_lig,name_res_lig,pos_LIG,pos_lig_NoWat,num_res_COFTOT = GET_RES_from_top (top_file)
		num_res_PROT = num_res_prot + num_res_COFTOT

		print("\t... cpptrajs")

		print("\t\t... removing waters")
		outfile = cpptraj_rmWat(lig)

		print("\t\t... separating ligands")
		solo_trajs = cpptraj_LigByLig(outfile)

		print("\t\t... extracting PDB from first image")
		outPDBs = cpptraj_OutPDBFirst(solo_trajs)

		print("\t\t... checking reactive trajs")
		dis_prot_lig_min = 9
		countlig = num_res_PROT + 1
		if_reactive = open("reactive_trajs_First.txt", "w")
		fileBS      = open ('lig_' + str(num_lig) + '_BS.dat',"w")
		for pdb in outPDBs:

			x_prot = []
			y_prot = []
			z_prot = []

			x_lig = []
			y_lig = []
			z_lig = []

			at_name = []
			res_name = []

			atom_name_lig = []
			resi_name_lig = "UNK"

			res_num = []

			num_atoms, num_atoms_lig, resi_name_lig, at_CA = read_prot(pdb)
			num_resBS, res_inBS = calc_BindingSite(num_lig, num_atoms, num_atoms_lig, dis_prot_lig_min, at_CA)

			if num_resBS == 0:
				if_reactive.write("0-{} --> Non-Reactive\n".format(countlig))
			else:
				if_reactive.write("1-{} --> Reactive\n".format(countlig))

			countlig += 1
		if_reactive.close()
		fileBS.close()

		os.chdir("../..")
	os.chdir("..")
	binding_sites_ALL(lig, dyn_num)	


	
