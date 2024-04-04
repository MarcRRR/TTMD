

def GET_RES_from_top (top_name) :

  protein = ['SER','THR','GLN','ASN','TYR','CYS','CYX','CYM','GLY',  \
             'ALA','VAL','LEU','ILE','MET','PRO','PHE','TRP',        \
             'GLU','GLH','ASP','ASH','LYS','ARG','HIE','HID','HIP',  \
             'PHE','TYR','TRP','ACE','NME','NHE','GLP','NHE' ]
  solvent   = ['WAT','OPC' ]
  ions      = ['Cl-','Na+' ]
  membrane  = ['PA','PC','OL']
  cofactor  = ['GTP','ATP']
  ionscofac = ['MG','CA']

  space = ' '
  pos_LIG_line = ''

  sys_top       = open (top_name,"r")

  for line in sys_top :
    if "%FLAG POINTERS" in line :
      line_read     = sys_top.readline ().replace("\n","")
      line_read1    = sys_top.readline ().replace("\n","")
      line_read2    = sys_top.readline ().replace("\n","")
      line_read_sp1 = line_read1.split()
      line_read_sp2 = line_read2.split()
      num_atoms_SYS = int (line_read_sp1 [0])
      num_res_SYS   = int (line_read_sp2 [1])
      num_lines     = num_res_SYS  / 20
      num_lines_int = num_res_SYS // 20
      if (num_lines - num_lines_int) != 0 :
        num_lines_int += 1
      break

  print ( ' ... Number of TOTAL atoms      : {:6d} '.format(num_atoms_SYS ) )
  print ( ' ... Number of TOTAL residues   : {:6d} '.format(num_res_SYS   ) )

  sys_top       = open (top_name,"r")
  residues = []
  pos_LIG  = []
  num_res_prot   = 0
  num_res_LIG    = 0
  num_res_SOL    = 0
  num_res_ION    = 0
  num_res_MEM    = 0
  num_res_COFTOT = 0
  num_res_COFAC  = 0
  num_res_IonsCOFAC = 0

  for line in sys_top :
    if "FLAG RESIDUE_LABEL" in line :
      line_read = sys_top.readline ().replace("\n","")
      for values in range (num_lines_int) :
        line_read    = sys_top.readline ().replace("\n","")
        line_read_sp = line_read.split()
        num_val      = len (line_read_sp)
        for num in range (num_val) :
          residues.append(line_read_sp[num])
      break

  indx = 0
  for res in residues :
    indx += 1
    if res in protein :
      num_res_prot += 1
    elif res in solvent :
      num_res_SOL += 1
    elif res in ions :
      num_res_ION += 1
    elif res in cofactor :
      num_res_COFAC += 1
    elif res in ionscofac :
      num_res_IonsCOFAC += 1
    elif res in membrane :
      num_res_MEM += 1
    else:
      num_res_LIG += 1
      name_res_LIG = res
      pos_LIG.append(indx)
      pos_LIG_char = str(indx)
      pos_LIG_line += pos_LIG_char + space

  if num_res_COFAC != 0 :
    num_res_COFTOT += num_res_COFAC
  if num_res_IonsCOFAC != 0 :
    num_res_COFTOT += num_res_IonsCOFAC

  pos_lig_NoWat  = []
  pos_LIG_LNoWat = ''
  for lig in range (num_res_LIG) :
    indx            = num_res_prot + num_res_COFTOT + lig + 1
    pos_LIG_char    = str(indx)
    pos_LIG_LNoWat += pos_LIG_char + space
    pos_lig_NoWat.append(indx)
  pos_LIG_LNoWat +=  '\n'

  print ( ' ... Number of PROTEIN residues : {:6d} '.format(num_res_prot ) )
  print ( ' ... Number of LIGAND  residues : {:6d} '.format(num_res_LIG  ) )
  print ( ' ... Number of IONS    residues : {:6d} '.format(num_res_ION  ) )
  print ( ' ... Number of SOLVENT residues : {:6d} '.format(num_res_SOL  ) )
  if num_res_MEM != 0 :
    print ( ' ... Number of MEMBRANE residues : {:6d} '.format(num_res_MEM   ) )
  if num_res_COFAC != 0 :  
    print ( ' ... Number of COFACTS residues : {:6d} '.format(num_res_COFAC ) )
  if num_res_IonsCOFAC != 0 :  
    print ( ' ... Number of IonsCOF residues : {:6d} '.format(num_res_IonsCOFAC ) )
  print ( ' ... LIGAND residue name        : {:>6} '.format(name_res_LIG ) )
  print ( ' ... LIGAND Positions           : ' )
  print ( ' ...  {:>6} '.format(pos_LIG_line   ) )
  print ( ' ...  {:>6} '.format(pos_LIG_LNoWat ) )

  return num_res_prot,num_res_LIG,name_res_LIG,pos_LIG,pos_lig_NoWat,num_res_COFTOT 

