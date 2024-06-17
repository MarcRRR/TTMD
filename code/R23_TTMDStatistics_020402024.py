import os
import subprocess
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import sys


def graph_affinity(list_x, list_y, title, label, name, slope_start=0):
	


		#### this function plots average IFPcs vs time

	fig, axs = plt.subplots(nrows=1, ncols=1)

	first_last_t = [list_x[0], list_x[-1]]
	axs.scatter(list_x, list_y)
	first_last_score = [slope_start, list_y[-1]]

	try:
		f = np.poly1d(np.polyfit(list_x, list_y, 1))
		slope = f[0]
		axs.plot(list_x, f(list_x), ls='--', color='red', label="MS = {:.5f}".format(slope))


	except Exception:
		slope = 'Impossible to calculate MS: required at least 2 TTMD steps'

	axs.set_title(title)
	axs.set_xlabel('Temperature (K)')
	axs.set_ylabel(label)

	axs.set_xlim(first_last_t)
	axs.legend()

	fig.savefig(name , dpi=300)

	plt.close()
	return slope



def graph_rmsd(ligand_coords, receptor_coords, temp, title, label, name):
	

	
	#### this function plots average IFPcs vs time

	plt.figure()
	plt.rc("font", size=15)
	plt.plot(temp, lig_coords, label = "ligand")
	plt.plot(temp, rec_coords, label = "backbone")
	
	plt.title(title)
	plt.xlabel('Temperature')
	plt.ylabel(label)
	plt.legend()
	plt.savefig(name , dpi=300)

	plt.close()
	return

ligs = []
DIR      = os.getcwd()
list_DIR = os.listdir(DIR)

for files in list_DIR :
	exist_dir = os.path.isdir(files)
	if exist_dir and files not in ["__pycache__", "graphs"]:
		ligs.append(files)
ligs.sort()

if os.path.exists("graphs"):
	os.system("rm -r graphs")
os.system("mkdir graphs")

os.system("mkdir graphs/affinities")
os.system("mkdir graphs/rmsd")

for lig in ligs:
	
		

	os.system("mkdir graphs/affinities/{}".format(lig))
	os.system("mkdir graphs/rmsd/{}".format(lig))

	with open("{}_BS_all.txt".format(lig), "r") as BSf:
		info = {}
		title = BSf.readline().strip().split("\t")
		title.append("MS")
		for line in BSf:
			line = line.strip().split()
			lignum, dynnum, BS = line[0], int(line[1]), line[2]
			if lignum not in info:
				info[lignum] = {}
			info[lignum][dynnum] = BS
	print(info)
		
	if os.path.exists("{}_final_table.tsv".format(lig)):
		os.system("rm {}_final_table.tsv".format(lig))

	final_BS = open("{}_final_table.tsv".format(lig), "x")

	with open("out_{}_sim.txt".format(lig), "r") as f:
		energies = []
		slopes = []
		for line in f:
			line = line.strip().split()
			if line != [] and line[0][:3] == "lig":
				energies.append(line)

		prev = 0
		dyn = 1

		temp = []
		for i in range(300, 452, 1):
			temp.append(i)

		for line in energies:
			lignum = line[0]
			if int(lignum[-3:]) < prev:
				dyn += 1
			ifps = line[1:]
			new_ifps = []
			for num in ifps:
				new_ifps.append(float(num))

			title = 'Affinity of {} in dyn {}'.format(lignum, dyn)
			label = 'Similarity'
			name = "{}_D{}.png".format(lignum, dyn)
			slope_start = -1

			slope = graph_affinity(temp, new_ifps, title, label, name, slope_start)
			final_BS.write("{}\t{}\t{}\t{}\n".format(lignum, dyn, info[lignum][dyn], slope))
			os.system("mv {} graphs/affinities/{}".format(name, lig))
			
			prev = int(lignum[-3:])



	final_BS.close()

	os.chdir(lig)
	dyn_num = int(subprocess.check_output("ls | grep '^dyn_*' | wc -l", shell=True).decode("utf-8"))
	for dyn in range(1, dyn_num+1):
		os.chdir("dyn_{}/TTMD".format(dyn))
		files = subprocess.check_output('ls | grep "^rmsd_"', shell=True).decode("utf-8").split()
		for rmsd in files:

			lig_rmsd, rec_rmsd = [], []

			f = open(rmsd, "r")
			bb = "bb"+rmsd[4:]
			bb_f = open(bb, "r")
			
			lig_coords = []
			rec_coords = []

			for line in f:
				line = line.strip().split()
				try:
					lig_rmsd.append(float(line[1]))
				except:
					continue

			for line in bb_f:
				line = line.strip().split()
				try:
					rec_rmsd.append(float(line[1]))
				except:
					continue
			f.close()
			bb_f.close()

			for i in range(0, len(lig_rmsd), 66):
				
				lig_coords.append(lig_rmsd[i])

				rec_coords.append(rec_rmsd[i])


			title = "rmsd of {} in TTMD".format(rmsd[5:12])
			label = "RMSD (\305)"
			name = "{}_D{}.png".format(rmsd[5:12], dyn)
			graph_rmsd(lig_coords, rec_coords, temp, title, label, name)
			
			os.system("mv {} ../../../graphs/rmsd/{}".format(name, lig))

		os.chdir("../..")

	os.chdir("..")
			
