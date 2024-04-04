import os
import sys
import subprocess


def create_amber(lig, params):		

# creation of the file that will be run by amber
	
	with open(lig+"_1_ttmd.in", "w") as f:
		
		ntpr = int(params["nstlim"]/1000)
		ntwr = int(params["nstlim"]/2)
		ntwe = ntwr
# first we write the header

		f.write(''' TTMD ''' +'['+str(params["Ti"])+' '+str(params["Tf"])+' '+\
str(params["nsteps"] )+' '+str(params["nstlim"])+']' ''' 
 &cntrl
   imin = 0, irest = 0, ntx = 5, nmropt = 1
   ntc = 2, ntf = 2, dt = '''+str(params["dt"])+''',
   cut = 11.0, fswitch = 8.0, iwrap=1,
   nstlim = '''+str(params["nstlim"])+''',
   ntb=1, ntp = 0,
   ntpr = '''+str(ntpr)+''',
   ntt = 3, gamma_ln = 3.0, ig = -1,
   temp0 = '''+str(params["Ti"])+''',
   ntwv = 0, ntwe = '''+str(ntwe)+''', ntwx = 5000, ntwr = '''+str(ntwr)+''',
   igamd = 3, iE = 1, irest_gamd = 1,
   ntcmdprep = 0, ntcmd = 0,
   ntebprep  = 0, nteb  = 0,
   ntave     =  126312,
   sigma0P =  6.0, sigma0D = 6.0,
 &end\n''')

# then we write the gradual temperature raising. The time variable is used to take the intervals of time for the temperature augment and MD analysis.

		Ti = params["Ti"]
		Tf = params["Tf"]

		step = int((Tf-Ti)/params["nsteps"])
		time = 0
		timestep = int(params["nstlim"]/params["nsteps"])

		for i in range(Ti, Tf, step):
			if time == 0:
				f.write(" &wt type='TEMP0', istep1="+str(time)+",istep2="+str(time+timestep)+",value1="+str(i)+".,value2="+str(i+step)+"., &end\n")
				time += timestep
			else:
				f.write(" &wt type='TEMP0', istep1="+str(time+1)+",istep2="+str(time+timestep)+",value1="+str(i)+".,value2="+str(i+step)+"., &end\n")
				time += timestep

		f.write(" &wt type='END'  &end\n")




def check(lig, params):	# checks all needed files are present and runs create_amber()


# checking if the ligand name exists and how many dynamics we have processed

	if os.path.exists(lig):
		os.chdir(lig)
		print("\nNow working on ligand {}\n".format(lig))
		dyn_num = int(subprocess.check_output("ls | grep 'dyn*' | wc -l", shell=True))
		for i in range(1, dyn_num+1):
			

# we create a TTMD folder for our analyses inside of the dynamic folder. We remove TTMD if it is already existing for a new analysis
			os.chdir("dyn_"+str(i))
			print("\n ... Now working on dyn_{}\n".format(i))
			exist_dir = os.path.isdir("TTMD")
			if exist_dir :
				try:
					print(" ... removing TTMD")
					os.system("rm -f -r TTMD")
					os.mkdir ('TTMD')
				except:
					os.mkdir ('TTMD')

# checking if a dyn folder is present to take needed files. if not, the files are present in the working directory

			if params["HMR"]:
				top_file = lig+"_HMR_vdWRep.top"
			else:
				top_file = lig+"_vdWRep.top"



			if os.path.exists("dyn_"+str(i)+"/dyn"):

				try:
					os.system("cp dyn/"+top_file+" TTMD")
				except:
					print("no topology file found in dyn_"+str(i))
				
					
				rst_files = subprocess.check_output("ls *.rst").decode("utf-8").strip().split()
				maxnum_rst = 0
				maxnum_gam = 0
				for f in rst_files:
					if f[-7:] == "dyn.rst":
						num_file = f.split("_")[1][0]
						if int(num_file) > maxnum_rst:
							maxnum_rst = int(num_file)
					if f[-12:] == "dyn.gamd_rst":
						num_file = f.split("_")[1][0]
						if int(num_file) > maxnum_gam:
							maxnum_gamd = int(num_file)
				if maxnum_rst != maxnum_gamd :
					print(" >>> Error in rst/gamdrst files : EXIT ")
					exit ()

				rst_file = "{}_{}_dyn.rst".format(lig, maxnum_rst)
				gam_file = "{}_{}_dyn.gamd_rst".format(lig, maxnum_gamd)


				try:
					os.system("cp dyn/"+rst_file+" TTMD")
					os.system("cp dyn/"+gam_file+" TTMD")
				except:
					print("no rst file found in dyn_"+str(i))



			else:

				try:
                                        os.system("cp "+top_file+" TTMD")
				except:
                                        print("no topology file found in dyn_"+str(i))
                                

				rst_files = os.listdir("./")
				maxnum_rst = 0
				maxnum_gam = 0
				for f in rst_files:
					if f[-7:] == "dyn.rst":
						num_file = f.split("_")[1][0]
						if int(num_file) > maxnum_rst:
							maxnum_rst = int(num_file)
					if f[-12:] == "dyn.gamd_rst":
						num_file = f.split("_")[1][0]
						if int(num_file) > maxnum_gam:
							maxnum_gamd = int(num_file)
				if maxnum_rst != maxnum_gamd :
					print(" >>> Error in rst/gamdrst files : EXIT ")
					exit ()

				rst_file = "{}_{}_dyn.rst".format(lig, maxnum_rst)
				gam_file = "{}_{}_dyn.gamd_rst".format(lig, maxnum_gamd)

				try:
                                        os.system("cp "+rst_file+" TTMD")
                                        os.system("cp "+gam_file+" TTMD")
				except:
                                        print("no rst file found in dyn_"+str(i))

# we go in TTMD and create the amber file with create_amber(). Then we go out of dyn_n to perform another iteration onto the next dynamic if present
			os.chdir("TTMD")
			create_amber(lig, params)
			os.chdir("../..")

# this line makes us go to the ligands directory so we can use a for loop iterating over different ligands to call this function 
		os.chdir("..")


def GEN_batch_g6gpu03_GaMD0 (name_frag,filename,top_file_name,rst_file,NUM_dyn) :
  send_fout = open('en_{}_DYN_g6gpu03'.format(name_frag),'w')
  send_fout.write('#!/bin/bash                            '+ter)
  send_fout.write('#######################################'+ter)
  send_fout.write("# Opcions i parametres de l'SGE        "+ter)
  send_fout.write('#######################################'+ter)
  send_fout.write('# (1) Nom del treball                  '+ter)
  send_fout.write('#$ -N {}D{} '.format(name_frag,NUM_dyn)+ter)
  send_fout.write('# (2) Recursos sol.licitats            '+ter)
  send_fout.write('##$ -l h_rt                            '+ter)
  send_fout.write('##$ -l mem_free                        '+ter)
  send_fout.write('#$ -pe gpu 1                           '+ter)
  send_fout.write('##$ -l exclusive=true                  '+ter)
  send_fout.write('# (3) Fitxers de sortida               '+ter)
  send_fout.write('#$ -cwd                                '+ter)
  send_fout.write('#$ -o {}.out          '.format(name_frag)+ter)
  send_fout.write('#$ -e {}.err          '.format(name_frag)+ter)
  send_fout.write('#$  -S /bin/bash                       '+ter)
  send_fout.write('# (4) Envia un mail                    '+ter)
  send_fout.write('##$ -m e                               '+ter)
  send_fout.write('##$ -M jaime.rubio@ub.edu              '+ter)
  send_fout.write('#######################################'+ter)
  send_fout.write('# Entorn de usuari                     '+ter)
  send_fout.write('#######################################'+ter)
  send_fout.write('# Es carreguen els moduls              '+ter)
  send_fout.write('source /etc/profile                    '+ter)
  send_fout.write('env                                    '+ter)
  send_fout.write('module load amber/20_cuda9.0_ompi_gcc-5.5.0 '+ter)
  send_fout.write('#######################################'+ter)
  send_fout.write('echo "SGE_O_WORKDIR : $SGE_O_WORKDIR"  '+ter)
  send_fout.write('echo "NSLOTS : $NSLOTS"                '+ter)
  send_fout.write('echo "TMP DIR : $TMPDIR"               '+ter)
  send_fout.write('echo " AMBERHOME   : $AMBERHOME"       '+ter)
  send_fout.write('#######################################'+ter)
  send_fout.write('#######################################'+ter)
  send_fout.write('# Calcul                               '+ter)
  send_fout.write('#######################################'+ter)
  send_fout.write('export RUN=/home/g6jaime/jaime/run     '+ter)
  send_fout.write('echo "RUN DIR : $RUN "                 '+ter)
  send_fout.write('cd $TMPDIR                             '+ter)
  send_fout.write('export CUDA_VISIBLE_DEVICES=`cat $TMPDIR/.gpus` '+ter)
  send_fout.write('#                                      '+ter)
  restarI_file = rst_file + '.gamd_rst'
  send_fout.write(' cp $SGE_O_WORKDIR/{}  $TMPDIR/gamd-restart.dat  '.format(restarI_file)+ter)
  send_fout.write(' $RUN/run_pmemd20_gpu_2021_nc    {}  {}  {}   rst'.format(filename,top_file_name,rst_file))
  send_fout.close()
  os.system('chmod u+x en_{}_DYN_g6gpu03'.format(name_frag))

  return


def GEN_batch_g6gpu09_GaMD0 (name_frag,filename,top_file_name,rst_file,NUM_dyn) :
  send_fout = open('en_{}_DYN_g6gpu09'.format(name_frag),'w')
  send_fout.write('#!/bin/bash                            '+ter)
  send_fout.write('#######################################'+ter)
  send_fout.write("# Opcions i parametres de l'SGE        "+ter)
  send_fout.write('#######################################'+ter)
  send_fout.write('# (1) Nom del treball                  '+ter)
  send_fout.write('#$ -N {}D{} '.format(name_frag,NUM_dyn)+ter)
  send_fout.write('# (2) Recursos sol.licitats            '+ter)
  send_fout.write('##$ -l h_rt                            '+ter)
  send_fout.write('##$ -l mem_free                        '+ter)
  send_fout.write('#$ -pe gpu 1                           '+ter)
  send_fout.write('##$ -l exclusive=true                  '+ter)
  send_fout.write('# (3) Fitxers de sortida               '+ter)
  send_fout.write('#$ -cwd                                '+ter)
  send_fout.write('#$ -o {}.out          '.format(name_frag)+ter)
  send_fout.write('#$ -e {}.err          '.format(name_frag)+ter)
  send_fout.write('#$  -S /bin/bash                       '+ter)
  send_fout.write('# (4) Envia un mail                    '+ter)
  send_fout.write('##$ -m e                               '+ter)
  send_fout.write('##$ -M jaime.rubio@ub.edu              '+ter)
  send_fout.write('#######################################'+ter)
  send_fout.write('# Entorn de usuari                     '+ter)
  send_fout.write('#######################################'+ter)
  send_fout.write('# Es carreguen els moduls              '+ter)
  send_fout.write('source /etc/profile                    '+ter)
  send_fout.write('env                                    '+ter)
  send_fout.write('module load amber/20_cuda9.0_ompi      '+ter)
  send_fout.write('#######################################'+ter)
  send_fout.write('echo "SGE_O_WORKDIR : $SGE_O_WORKDIR"  '+ter)
  send_fout.write('echo "NSLOTS : $NSLOTS"                '+ter)
  send_fout.write('echo "TMP DIR : $TMPDIR"               '+ter)
  send_fout.write('echo " AMBERHOME   : $AMBERHOME"       '+ter)
  send_fout.write('#######################################'+ter)
  send_fout.write('#######################################'+ter)
  send_fout.write('# Calcul                               '+ter)
  send_fout.write('#######################################'+ter)
  send_fout.write('export RUN=/home/g6jaime/jaime/run     '+ter)
  send_fout.write('echo "RUN DIR : $RUN "                 '+ter)
  send_fout.write('cd $TMPDIR                             '+ter)
  send_fout.write('export CUDA_VISIBLE_DEVICES=`cat $TMPDIR/.gpus` '+ter)
  send_fout.write('#                                      '+ter)
  restarI_file = rst_file + '.gamd_rst'
  send_fout.write(' cp $SGE_O_WORKDIR/{}  $TMPDIR/gamd-restart.dat  '.format(restarI_file)+ter)
  send_fout.write(' $RUN/run_pmemd20_gpu_2021_nc    {}  {}  {}   rst'.format(filename,top_file_name,rst_file))

  send_fout.close()
  os.system('chmod u+x en_{}_DYN_g6gpu09'.format(name_frag))

  return

def GEN_batch_iqtc10_GaMD0 (name_frag,filename,top_file_name,rst_file,NUM_dyn) :
  send_fout = open('en_{}_DYN_iqtc10'.format(name_frag),'w')
  send_fout.write('#!/bin/bash                            '+ter)
  send_fout.write('#######################################'+ter)
  send_fout.write("# Opcions i parametres de l'SGE        "+ter)
  send_fout.write('#######################################'+ter)
  send_fout.write('# (1) Nom del treball                  '+ter)
  send_fout.write('#$ -N {}D{} '.format(name_frag,NUM_dyn)+ter)
  send_fout.write('# (2) Recursos sol.licitats            '+ter)
  send_fout.write('#$ -l iqtcgpu=1                        '+ter)
  send_fout.write('# (3) Fitxers de sortida               '+ter)
  send_fout.write('#$ -cwd                                '+ter)
  send_fout.write('#$ -o {}.out          '.format(name_frag)+ter)
  send_fout.write('#$ -e {}.err          '.format(name_frag)+ter)
  send_fout.write('#$  -S /bin/bash                       '+ter)
  send_fout.write('# (4) Envia un mail                    '+ter)
  send_fout.write('##$ -m e                               '+ter)
  send_fout.write('##$ -M jaime.rubio@ub.edu              '+ter)
  send_fout.write('#######################################'+ter)
  send_fout.write('# Entorn de usuari                     '+ter)
  send_fout.write('#######################################'+ter)
  send_fout.write('# Es carreguen els moduls              '+ter)
  send_fout.write('source /etc/profile                    '+ter)
  send_fout.write('env                                    '+ter)
  send_fout.write('module load amber/20_cuda11.0_ompi     '+ter)
  send_fout.write('#######################################'+ter)
  send_fout.write('echo "SGE_O_WORKDIR : $SGE_O_WORKDIR"  '+ter)
  send_fout.write('echo "NSLOTS : $NSLOTS"                '+ter)
  send_fout.write('echo "TMP DIR : $TMPDIR"               '+ter)
  send_fout.write('echo " AMBERHOME   : $AMBERHOME"       '+ter)
  send_fout.write('#######################################'+ter)
  send_fout.write('#######################################'+ter)
  send_fout.write('# Calcul                               '+ter)
  send_fout.write('#######################################'+ter)
  send_fout.write('export RUN=/home/g6jaime/jaime/run     '+ter)
  send_fout.write('echo "RUN DIR : $RUN "                 '+ter)
  send_fout.write('cd $TMPDIR                             '+ter)
  send_fout.write('export CUDA_VISIBLE_DEVICES=`cat $TMPDIR/.gpus` '+ter)
  send_fout.write('#                                      '+ter)
  restarI_file = rst_file + '.gamd_rst'
  send_fout.write(' cp $SGE_O_WORKDIR/{}  $TMPDIR/gamd-restart.dat  '.format(restarI_file)+ter)
  send_fout.write(' $RUN/run_pmemd20_gpu_2021_nc    {}  {}  {}   rst'.format(filename,top_file_name,rst_file))

  send_fout.close()
  os.system('chmod u+x en_{}_DYN_iqtc10'.format(name_frag))

  return


def GEN_batch_mfg_GaMD0 (name_frag,filename,top_file_name,rst_file,NUM_dyn) :
  send_fout = open('en_{}_DYN_mfg'.format(name_frag),'w')
  send_fout.write('#!/bin/bash                             '+ter)
  send_fout.write('#SBATCH --job-name={}D{} '.format(name_frag,NUM_dyn)+ter)
  send_fout.write('#SBATCH --gres=gpu:1                           '+ter)
  send_fout.write('#SBATCH --nodes=1                  '+ter)
  send_fout.write('#SBATCH --ntasks=1               '+ter)
  send_fout.write('#SBATCH -e {}.err          '.format(name_frag)+ter)
  send_fout.write('#SBATCH -o {}.out          '.format(name_frag)+ter)
  send_fout.write('echo $CUDA_VISIBLE_DEVICES > cuda  '+ter)
  send_fout.write('export SCRATCH=scratch$CUDA_VISIBLE_DEVICES           '+ter)
  send_fout.write('export PATH=$PATH:$PWD              '+ter)
  send_fout.write('cd $SLURM_SUBMIT_DIR/              '+ter)
  send_fout.write('source /programas/amber20/amber20/amber.sh            '+ter)
  send_fout.write('export PATH=$PATH:/usr/local/cuda-11.2/bin:/programas/amber20/amber20/bin/         '+ter)
  send_fout.write('export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda-11.2/lib64:/programas/amber20/amber20/lib/      '+ter)
  send_fout.write('export CUDA_HOME=/usr/local/cuda-11.2          '+ter)
  send_fout.write('rm -f -r  /$SCRATCH/$SLURM_JOB_NAME/               '+ter)
  send_fout.write('mkdir     /$SCRATCH/$SLURM_JOB_NAME/               '+ter)
  send_fout.write('rsync -avz   $SLURM_SUBMIT_DIR/* /$SCRATCH/$SLURM_JOB_NAME/    '+ter)
  send_fout.write('export RUN=/home/chema/covid/programs/       '+ter)
  send_fout.write('cd /$SCRATCH/$SLURM_JOB_NAME/                '+ter)
  send_fout.write('#                                            '+ter)
  restarI_file = rst_file + '.gamd_rst'
  send_fout.write(' cp $SLURM_SUBMIT_DIR/{}  ./gamd-restart.dat  '.format(restarI_file)+ter)
  send_fout.write(' $RUN/run_pmemd20    {}  {}  {}   rst'.format(filename,top_file_name,rst_file))

  send_fout.close()
  os.system('chmod u+x en_{}_DYN_mfg'.format(name_frag))

  return


def run(lig): # Generates executable files for calculation sending to different queues

	os.chdir(lig)
	dyn_num = int(subprocess.check_output("ls | grep 'dyn*' | wc -l", shell=True))
	
	for i in range(1, dyn_num+1):

		print("  Working on dyn_{}\n".format(i))
		filename = "{}_1_ttmd".format(lig)
		os.chdir("dyn_{}/TTMD".format(i))
		top_file = subprocess.check_output("ls *.top", shell=True).decode("utf-8")
		rst_file = subprocess.check_output("ls *.rst", shell=True).decode("utf-8")

		top_file = top_file[0:len(top_file)-5]
		rst_file = rst_file[0:len(rst_file)-5]

		print ( '    ... ... ... ... Running Gen Batch g6gpu03 ' )
		GEN_batch_g6gpu03_GaMD0 (lig,filename,top_file,rst_file,i)
		print ( '    ... ... ... ... Running Gen Batch g6gpu09 ' )
		GEN_batch_g6gpu09_GaMD0 (lig,filename,top_file,rst_file,i)
		print ( '    ... ... ... ... Running Gen Batch iqtc10  ' )
		GEN_batch_iqtc10_GaMD0  (lig,filename,top_file,rst_file,i)
		print ( '    ... ... ... ... Running Gen Batch mfg     ' )
		GEN_batch_mfg_GaMD0     (lig,filename,top_file,rst_file,i)


		print(  '    ... ... ... ... Now generating batch files\n')
		


		file_en_batch_all_g6gpu03.write('cd {}/dyn_{}/TTMD'.format(lig,i)+ter)
		file_en_batch_all_g6gpu03.write('qsub -q g6gpu03.q en_{}_DYN_g6gpu03'.format(lig)+ter)
		file_en_batch_all_g6gpu03.write('cd ../../../'+ter)
		file_en_batch_all_g6gpu03.write('sleep 5'+ter)


		file_en_batch_all_g6gpu09.write('cd {}/dyn_{}/TTMD'.format(lig,i)+ter)
		file_en_batch_all_g6gpu09.write('qsub -q g6gpu03.q en_{}_DYN_g6gpu09'.format(lig)+ter)
		file_en_batch_all_g6gpu09.write('cd ../../../'+ter)
		file_en_batch_all_g6gpu09.write('sleep 5'+ter)

		file_en_batch_all_iqtc10.write('cd {}/dyn_{}/TTMD'.format(lig,i)+ter)
		file_en_batch_all_iqtc10.write('qsub -q iqtc10.q  en_{}_DYN_iqtc10'.format(lig)+ter)
		file_en_batch_all_iqtc10.write('cd ../../../'+ter)
		file_en_batch_all_iqtc10.write('sleep 5'+ter)

		file_en_batch_all_mfg.write('cd {}/dyn_{}/TTMD'.format(lig,i)+ter)
		file_en_batch_all_mfg.write('sbatch  en_{}_DYN_mfg'.format(lig)+ter)
		file_en_batch_all_mfg.write('cd ../../../'+ter)
		file_en_batch_all_mfg.write('sleep 5'+ter)


		os.chdir("../..")
	os.chdir("..")
	return



def config(fil):
	
	params = {"nstlim":'', "dt" : '', "Ti" : '',"Tf" : '', "nsteps" : ''}
	ite = 0

	with open(fil, "r") as ini:
		for line in ini:
			line = line.strip().split()
			if len(line) != 0 and line[-1].replace(".", "").isnumeric():
				try:
					num = int(line[-1])
				except:
					num = float(line[-1])
				params[line[0]] = num
				ite += 1

	return params
			





#TESTING


# variables to use

global ter

ter = '\n'

# assuming only the program and the ligands directories are present in the working directory

ligs = []
DIR      = os.getcwd()
list_DIR = os.listdir(DIR)
for files in list_DIR :
	exist_dir = os.path.isdir(files)
	if exist_dir :
		ligs.append(files)
ligs.sort()


#ligs = subprocess.check_output("ls | grep -v R20_TTMD_01022024.py | grep -v input_file", shell=True).strip().split()
params = eval(open('input_file', 'r').read())

# execution for every ligand present in the folder. The ligand is decoded since the ouput from subprocess.check_output is coded

file_en_batch_all_g6gpu03  = open ('en_batch_all_g6gpu03_ttmd', "w")
file_en_batch_all_g6gpu09  = open ('en_batch_all_g6gpu09_ttmd', "w")
file_en_batch_all_iqtc10   = open ('en_batch_all_iqtc10_ttmd' , "w")
file_en_batch_all_mfg      = open ('en_batch_all_mfg_ttmd'    , "w")

for i in ligs:
	check(i, params)
	run(i)

file_en_batch_all_g6gpu03.close()
file_en_batch_all_g6gpu09.close()
file_en_batch_all_iqtc10.close()
file_en_batch_all_mfg.close()

os.system('chmod u+x en_batch_all_g6gpu03_ttmd')
os.system('chmod u+x en_batch_all_g6gpu09_ttmd')
os.system('chmod u+x en_batch_all_iqtc10_ttmd ')
os.system('chmod u+x en_batch_all_mfg_ttmd    ')

