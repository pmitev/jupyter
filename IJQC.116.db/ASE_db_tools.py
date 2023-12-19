#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:set foldmarker=#-#,#.# foldmethod=marker:
import numpy as np
from ase.visualize import view
import sys,os

e0=0.00552635;              # Vacuum permittivity constant C/V
ang2bohr= 1./0.52917721;    # Angsrom to bohr
au2D= 1./0.393430307;       # au to Debye
eV_a2au=0.0194469057312645; # force conversion


core_charges= {"H":1., "O":6., "W":-2., "Al":3., "S":6., "N":5.}


#================================================================
def db_read_data(dbname, row_id):
#-#
  import sqlite3
  from ase.io.jsonio import encode, decode
  from ase import Atom, Atoms

  db=  sqlite3.connect(dbname)
  cursor= db.cursor()
  
  # Fetch "data" from the database
  cursor.execute("SELECT data FROM systems WHERE id=?",row_id)
  row= cursor.fetchone()
  data= decode(row[0])
  db.commit()
  db.close()

  return data
#.#

#================================================================
def db_write_data(dbname, row_id, data):
#-#
  import sqlite3
  from ase.io.jsonio import encode, decode
  from ase import Atom, Atoms

  db=  sqlite3.connect(dbname)
  cursor= db.cursor()
  
  # Commit changes to the database
  cursor.execute('UPDATE systems SET data = ? WHERE id = ?',(encode(data),row_id[0]) )
  db.commit()
  db.close()
  
#.#

#================================================================
def add_VASP_structure(dbname, folder):
#-#
  from ase.io import read
  from ase.io.vasp import read_vasp_out
  import ase.db
  import sys, os, re

  start_dir= os.path.realpath(".")
  db=     ase.db.connect(dbname)
  os.chdir(folder)

  structure= read_vasp_out("OUTCAR")

  fin= open("INCAR");   INCAR=   fin.readlines(); fin.close()
  fin= open("KPOINTS"); KPOINTS= fin.readlines(); fin.close()
  fin= open("POTCAR");  POTCARf= fin.read();      fin.close(); POTCAR= re.findall('TITEL  = (.*)',POTCARf)
  fin= open("CONTCAR"); CONTCAR= fin.readlines(); fin.close()

  datain= {"INCAR"         : INCAR,
           "KPOINTS"       : KPOINTS,
           "POTCAR_TITEL"  : POTCAR,
           "CONTCAR"       : CONTCAR
          }
  str_path=os.path.realpath(".")
  os.chdir(start_dir)

  row_id= db.write(structure,data= datain, structure_path= str_path)
  print ("Structure added with id:",row_id)
  return row_id
#.#

#================================================================
def add_VASP_wannier(dbname, row_id, folder):
#-#
  from ase import Atom, Atoms

  start_dir= os.path.realpath(".")         
  row_id= (row_id,)
  os.chdir(folder)

  #structure_w= read("wannier90_centres.xyz")
  #del structure_w[[atom.index for atom in structure_w if atom.symbol!='X']] # Delete all atoms and leave the wannier centers "X"

  # Old wannier format is brocken - can't use the lines above
  structure_w= Atoms()
  fin= open("wannier90_centres.xyz")
  while True:
    line= fin.readline()
    if not line:
      break
    lsplit= line.split()
    if lsplit[0] == "X":
      structure_w.append(Atom(lsplit[0],lsplit[1:]))
  
  with open("wannier90.win") as file:
    data= file.read()
  tmps= data.partition("begin unit_cell_cart")[2].rpartition("end unit_cell_cart")[0].strip().split()
  #tmpf=[float(i) for i in tmps]
  tmpf= list(map(float,tmps))
  cell= np.array(tmpf).reshape((3,3))
  structure_w.set_cell(cell); structure_w.set_pbc(True)
  WANNIER_frac= structure_w.get_scaled_positions()
  
  datain={"WANNIER_frac"  : WANNIER_frac,
          "WANNIER_xyz"   : structure_w.get_positions(),
          "WANNIER_path"  : os.path.realpath(".")
          }
  
  realpath= os.path.realpath(".")
  os.chdir(start_dir)
  
  # Add/update data
  data= db_read_data(dbname, row_id)
  data.update(datain)
  db_write_data(dbname, row_id, data)
  
  import ase.db
  db= ase.db.connect(dbname)
  db.update(row_id[0],WANNIER_centers=len(WANNIER_frac),WANNIER_path= realpath)
  print ("[wannier]: Wannier data added for id:",row_id[0])
#.#

#================================================================
def add_VASP_bader(dbname, row_id, folder):
#-#
  from ase import Atom, Atoms
  
  start_dir= os.path.realpath(".")         
  row_id= (row_id,)
  os.chdir(folder)
  
  try:
    fin= open("INCAR");   INCAR=   fin.readlines(); fin.close()
  except:
    try:
      fin= open("../INCAR");   INCAR=   fin.readlines(); fin.close()
      print ("[bader]: WARNING!!! - INCAR file from parent folder")
    except:  
      print ("[bader]: WARNING!!! - Missing INCAR file")
      INCAR=""
  
  try:
    fin= open("KPOINTS"); KPOINTS= fin.readlines(); fin.close()
  except:
    try:
      fin= open("../KPOINTS"); KPOINTS= fin.readlines(); fin.close()
      print ("[bader]: WARNING!!! - KPOINTS file from parent folder")
    except:
      print ("[bader]: WARNING!!! - Missing KPOINTS file")
      KPOINTS=""
  
  # Read the data from DIP.dat
  bader_coord = np.genfromtxt("DIP.dat",skip_header=5, skip_footer=3, usecols=(1,2,3) ); bader_coord= bader_coord /ang2bohr
  bader_chg   = np.genfromtxt("DIP.dat",skip_header=5, skip_footer=3, usecols=(4)     )
  #bader_dipole= np.genfromtxt("DIP.dat",skip_header=5, skip_footer=1, usecols=(5,6,7) )
  bader_dipole= np.genfromtxt("DIP.dat",skip_header=5, skip_footer=3, usecols=(5,6,7), invalid_raise=False )
  
  datain={"BADER_INCAR"   : INCAR,
          "BADER_KPOINTS" : KPOINTS,
          "BADER_charges" : bader_chg,
          "BADER_xyz"     : bader_coord,
          "BADER_dipoles" : bader_dipole,
          "BADER_path"    : os.path.realpath(".")
          }
  
  realpath= os.path.realpath(".")
  os.chdir(start_dir)

  # Add/update data
  data= db_read_data(dbname, row_id)
  data.update(datain)
  db_write_data(dbname, row_id, data)
  
  import ase.db
  db= ase.db.connect(dbname)
  db.update(row_id[0],BADER_path= realpath)
  
  print ("[bader]: Bader data added for id:",row_id[0])
#.#

#================================================================
def add_VASP_DDEC6(dbname, row_id, folder):
#-#
  
  start_dir= os.path.realpath(".") 
  row_id= (row_id,)
  os.chdir(folder)
  
  try:
    fin= open("DDEC6_even_tempered_net_atomic_charges.xyz");   DDEC6=   fin.readlines(); fin.close()
  except:
    print ("[DDEC6]: missing DDEC6_even_tempered_net_atomic_charges.xyz file")
    sys.exit()
  
  try:
    fin= open("INCAR");   INCAR=   fin.readlines(); fin.close()
  except:  
    print ("[DDEC6]: WARNING!!! - Missing INCAR file")
    INCAR=""
  
  try:
    fin= open("KPOINTS"); KPOINTS= fin.readlines(); fin.close()
  except:
    print ("[DDEC6]: WARNING!!! - Missing KPOINTS file")
    KPOINTS=""
  
  # Read the data from DDEC6_even_tempered_net_atomic_charges.xyza
  for i in range(len(DDEC6)):
    if "dipole_x" in DDEC6[i]:
      ii= 0; charges=[]; dipoles=[]
      while True:
        ii= ii+1
        line= DDEC6[i+ii].split()
        if len(line) != 0 :
          charges.append(float(line[5]))
          dipoles.append(list(map(float,line[6:9])))
        else:
          break
  
  datain={"DDEC6_charges" : charges,
          "DDEC6_dipoles" : dipoles,
          "DDEC6_path"    : os.path.realpath(".")
          }
  
  realpath= os.path.realpath(".")
  os.chdir(start_dir)
 
  # Add/update data
  data= db_read_data(dbname, row_id)
  data.update(datain)
  db_write_data(dbname, row_id, data)
  
  import ase.db
  db= ase.db.connect(dbname)
  db.update(row_id[0],DDEC6_path= realpath)
  
  print ("[DDEC6]: DDEC6 data added for id:",row_id[0])
#.#

#================================================================
def add_OH_PES(dbname, row_id, PES_file, data_key):
#-#
  row_id= (row_id,)
  
  try:
    fin= open(PES_file);   tmp= fin.readlines(); fin.close()
  except:
    print ("ERROR: File",PES_file,"not found")
    sys.exit()
  
  # Read the data 
  PES= np.genfromtxt(PES_file)
  
  datain={ data_key : PES }
  
  # Add/update data
  data= db_read_data(dbname, row_id)
  data.update(datain)
  db_write_data(dbname, row_id, data)
  
  import ase.db
  db= ase.db.connect(dbname)
  db.update(row_id[0])
  
  print ("[PES]: data added for id:",row_id[0])
#.#

#================================================================
def add_DVR_results(dbname, row_id, DVR_res, Hind):
#-#
  from ase.io import read
  row_id= (row_id,)
  
  try:
    fin= open(DVR_res);   lines= fin.readlines(); fin.close()
  except:
    print ("ERROR: File",DVR_res,"not found")
    sys.exit()
 
  harm= list(filter(lambda x: "Harmonic frequency=" in x, lines))[0].split()[2]
  E0=   list(filter(lambda x: "E0=" in x, lines))[0].split()[1]
  E1=   list(filter(lambda x: "E1=" in x, lines))[0].split()[1]
  E2=   list(filter(lambda x: "E2=" in x, lines))[0].split()[1]
  
  # Fetch natoms
  natoms= read(dbname+"@"+str(row_id[0]-1))[0].get_number_of_atoms()
  
  # Fetch "data" from the database
  data= db_read_data(dbname, row_id)

  if "DVR_data" in data.keys():
    DVR_dat= data["DVR_data"]
  else:
    DVR_dat= np.zeros((natoms,4))
  
  DVR_dat[Hind]= [harm,E0,E1,E2]
  datain= {"DVR_data" : DVR_dat }
  
  # Add/update data
  data.update(datain)
  db_write_data(dbname, row_id, data)
  
  import ase.db
  db= ase.db.connect(dbname)
  db.update(row_id[0])
  
  print ("[DVR]: data added for id:",row_id[0])
#.#

#================================================================
#-#

#.#




#================================================================
def db_read_row(dbname,id):
#-# 
  import ase.db
  db= ase.db.connect(dbname)
  # Read all information for one row in the database
  row= db.get(id= id)
  
  # Create Atoms object and glu the data from the row for convenience
  structure= row.toatoms()
  structure.row= row
  return structure
#.#

#================================================================
def get_indexes_Ow(structure):
#-# 
  import numpy as np
  from ase.data import covalent_radii
  from ase.neighborlist import NeighborList

  covalent_radii[11]=1.4 # decrease a bit the Na radious
  
  # Identify molecules based on the covalenr radius * xxxx  =========================================================
  nl= NeighborList(covalent_radii[structure.get_atomic_numbers()]*1.15,skin=0, self_interaction=False,bothways=True);  nl.update(structure)
  
  # Make list for elements of interest
  indexO= [i for i,x in enumerate(structure.get_chemical_symbols()) if x == 'O'] # O indexes
  indexH= [i for i,x in enumerate(structure.get_chemical_symbols()) if x == 'H'] # H indexes
  
  indO_nH= [[i,len(np.intersect1d(nl.get_neighbors(i)[0],indexH))] for i in indexO ] # Array with O and corresponding H neighbors
  indexOw= [indxO for (indxO,nH) in indO_nH if nH == 2] # index of Ow
 
  return indexOw
#.#

#================================================================
def db_find_symmetry(dbname, row_id, symprec= 1.e-5):
#-# 
  structure= db_read_row(dbname, row_id)
  spacegroup, _, atom_labels= find_structure_symmetry(structure, symprec= symprec)

  datain= { "ATOM_labels" : atom_labels }

  # Add/update data
  row_id= (row_id,)
  data= db_read_data(dbname, row_id)
  data.update(datain)
  db_write_data(dbname, row_id, data)

  return spacegroup
#.#

#================================================================
def prop_H2O_wannier_dipoles(structure_in, core_charges,verbose= True):
#-# 
  from copy import deepcopy
  from ase import Atom, Atoms
  from ase.neighborlist import NeighborList

  structure= deepcopy(structure_in)
  index_Ow= get_indexes_Ow(structure)

  if "WANNIER_xyz" in structure.row.data.keys():
    nwannier= len(structure.row.data["WANNIER_xyz"])
    structure.extend(Atoms(symbols="W"*nwannier, positions= structure.row.data["WANNIER_xyz"], cell= structure.get_cell()))
  else:  
    print ("ERROR: Missing WANNIER_xyz")
    sys.exit()

  nl=NeighborList([1.2/2]*len(structure),skin=0, self_interaction=False,bothways=True);   nl.update(structure)

  Dipoles= {}
  for iOw in index_Ow:
    vCENTER=(structure.get_cell()[0]+structure.get_cell()[1]+structure.get_cell()[2])/2. # get translation vector to put the O in the center of the box
    structure.translate(-structure.get_positions()[iOw] + vCENTER)    # center the O in the box
    structure.set_scaled_positions(structure.get_scaled_positions()) # put back everything in the box
    
    D= [0.,0.,0.]
    chg_sum= core_charges[structure.get_chemical_symbols()[iOw]]
    if verbose:
      print ("  Origin at atom idx:", iOw," core_charge:", chg_sum)
      print ("="*44)
    for ip in nl.get_neighbors(iOw)[0]:
      v_OwP= structure.get_positions()[ip] - structure.get_positions()[iOw]
      chg= core_charges[structure.get_chemical_symbols()[ip]]; chg_sum= chg_sum + chg;
      D= D + chg*v_OwP*ang2bohr
      if verbose:
        print ("idx:{0:5}   Chg:{1:5}   dist: {2:6f}".format(ip, chg, np.linalg.norm(v_OwP)) )
    if verbose: 
      print ("="*44)
      print ("  Total Dipole= {0:.4f} [Debye]\n".format(np.linalg.norm(D*au2D)) )

    Dipoles[iOw]= [D*au2D]
    if chg_sum != 0:
      print ("ERROR: Total charge of the H2O molecule:", chg_sum,"\n       Make sure you have provided the correct correct core_charges")

  return Dipoles
#.#

#================================================================
def prop_H2O_bader_dipoles(structure_in, core_charges, verbose=True):
#-# 
  from copy import deepcopy
  from ase import Atom, Atoms
  from ase.neighborlist import NeighborList

  structure= deepcopy(structure_in)
  index_Ow= get_indexes_Ow(structure)
  natoms= len(structure)
  
  if "BADER_dipoles" in structure.row.data.keys():
    BADER_xyz=     structure.row.data["BADER_xyz"]
    BADER_dipoles= structure.row.data["BADER_dipoles"]
    BADER_charges= structure.row.data["BADER_charges"]
    nbdipoles= len(BADER_dipoles)
    structure.extend(Atoms(symbols="X"*nbdipoles, positions= BADER_xyz, cell= structure.get_cell()))
  else:
    print ("ERROR: Missing BADER_dipoles")
    sys.exit()

  nl=NeighborList([1.2/2]*len(structure),skin=0, self_interaction=False,bothways=True);   nl.update(structure)

  Dipoles= {}; Charges={}
  for iOw in index_Ow:
    vCENTER=(structure.get_cell()[0]+structure.get_cell()[1]+structure.get_cell()[2])/2. # get translation vector to put the O in the center of the box
    structure.translate(-structure.get_positions()[iOw] + vCENTER)    # center the O in the box
    structure.set_scaled_positions(structure.get_scaled_positions()) # put back everything in the box
    
    #D= [0.,0.,0.] # THIS WAS WRONG IN THE OLD CODES
    D= BADER_dipoles[iOw]
    chg_sum= core_charges[structure.get_chemical_symbols()[iOw]]
    if verbose:
      print ("  Origin at atom idx:", iOw, "core_charge:",chg_sum)
      print ("="*42)
    for ip in nl.get_neighbors(iOw)[0]:
      v_OwP= structure.get_positions()[ip] - structure.get_positions()[iOw]
      if ip < natoms:
        chg= core_charges[structure.get_chemical_symbols()[ip]]
        D= D + chg*v_OwP*ang2bohr 
      else:
        chg= BADER_charges[ip-natoms]
        D= D + chg*v_OwP*ang2bohr + BADER_dipoles[ip-natoms]
      chg_sum= chg_sum + chg
      if verbose:
        print ("idx:{0:5}   Chg:{1:8.4f}   dist: {2:6f}".format(ip, chg, np.linalg.norm(v_OwP)) )
    if verbose:
      print ("="*42)
      print ("  Total Dipole= {0:7.4f} [Debye]".format(np.linalg.norm(D*au2D)) )
      print ("  Total Charge= {0:7.4f} [e]\n".format(chg_sum) )
    Dipoles[iOw]= D*au2D; Charges[iOw]= chg_sum

  return Dipoles, Charges
#.#

#================================================================
def prop_H2O_DDEC6_dipoles(structure_in, verbose=True):
#-# 
  from copy import deepcopy
  from ase import Atom, Atoms
  from ase.neighborlist import NeighborList

  structure= deepcopy(structure_in)
  index_Ow= get_indexes_Ow(structure)
  
  if "DDEC6_dipoles" in structure.row.data.keys():
    DDEC6_charges= structure.row.data["DDEC6_charges"]
    DDEC6_dipoles= structure.row.data["DDEC6_dipoles"]
  else:
    print ("Missing DDEC6_dipoles")
    sys.exit()

  nl=NeighborList([1.2/2]*len(structure),skin=0, self_interaction=False,bothways=True);   nl.update(structure)

  Dipoles= {}
  for iOw in index_Ow:
    vCENTER=(structure.get_cell()[0]+structure.get_cell()[1]+structure.get_cell()[2])/2. # get translation vector to put the O in the center of the box
    structure.translate(-structure.get_positions()[iOw] + vCENTER)    # center the O in the box
    structure.set_scaled_positions(structure.get_scaled_positions()) # put back everything in the box
    
    D= DDEC6_dipoles[iOw]
    if verbose:
      print ("  Origin at atom idx:", iOw,"  Chg:", DDEC6_charges[iOw],"       |D|:",np.linalg.norm(DDEC6_dipoles[iOw]),"[Debye]")
      print ("="*68)
    for ip in nl.get_neighbors(iOw)[0]:
      v_OwP= structure.get_positions()[ip] - structure.get_positions()[iOw]
      chg= DDEC6_charges[ip]
      D= D + chg*v_OwP*ang2bohr + DDEC6_dipoles[ip]
      if verbose:
        print ("idx:{0:5}   Chg:{1:8.4f}   dist: {2:6f} [Å]  |D|: {3:.6f} [Debye]".format(ip, chg, np.linalg.norm(v_OwP),np.linalg.norm(DDEC6_dipoles[ip]) ) )
    if verbose:    
      print ("="*68)
      print ("  Total Dipole= {0:.4f} [Debye]\n".format(np.linalg.norm(D*au2D)) )
    Dipoles[iOw]= D*au2D

  return Dipoles
#.# 

#================================================================
def EF_H2O_bader(structure_in, core_charges, verbose=True):
#-# 
  from copy import deepcopy
  from ase.neighborlist import NeighborList
  from pymatgen.core import Lattice, Structure
  from pymatgen.analysis import ewald

  structure= deepcopy(structure_in)
  index_Ow= get_indexes_Ow(structure)

  if "BADER_charges" in structure.row.data.keys():
    BADER_charges= structure.row.data["BADER_charges"]
  else:
    print ("Missing BADER_charges")
    sys.exit()

  nl=NeighborList([1.2/2]*len(structure),skin=0, self_interaction=False,bothways=True);  nl.update(structure)
  indexH= [i for i,x in enumerate(structure.get_chemical_symbols()) if x == 'H'] # H indexes

  # Run Ewald
  print ("Running Ewald ...")
  lattice= Lattice(structure.get_cell())
  coords= structure.get_scaled_positions()
  structure.charges= [core_charges[x] for x in structure.get_chemical_symbols()] + BADER_charges
  struct= Structure(lattice,structure.get_chemical_symbols(), coords, site_properties={"charge" : structure.charges})
  Ew= ewald.EwaldSummation(struct, compute_forces=True)
  print ("Done")
  
  for iOw in index_Ow:
    print (" Ow idx:", iOw)
    print ("="*50)
  
    vCENTER=(structure.get_cell()[0]+structure.get_cell()[1]+structure.get_cell()[2])/2. # get translation vector to put the O in the center of the box
    structure.translate(-structure.get_positions()[iOw] + vCENTER)    # center the O in the box
    structure.set_scaled_positions(structure.get_scaled_positions()) # put back everything in the box
  
    neighbors= nl.get_neighbors(iOw)[0]
    indexH_loc= np.intersect1d(neighbors,indexH)
    vOH1= structure.get_positions()[indexH_loc[0]] - structure.get_positions()[iOw]
    vOH2= structure.get_positions()[indexH_loc[1]] - structure.get_positions()[iOw]
    vMid= (vOH1/np.linalg.norm(vOH1) + vOH2/np.linalg.norm(vOH2)); vMid= vMid/np.linalg.norm(vMid)
  
    for j in indexH_loc:
      vOH= structure.get_positions()[j] - structure.get_positions()[iOw]
      uOH= vOH/np.linalg.norm(vOH)
      OH= structure.get_distance(iOw,j,mic=True)
  
      # Total EF at H position 
      EF_H= Ew.forces[j]/structure.charges[j]*eV_a2au;
      if BADER_charges[j] > -0.01:
        print (" WARNING: H atom with idx:",j,"has no associated bader charge")
      print (" {:>2}{:<4}{:>2}{:<4} dist: {:.6f}   chg{:<4} {:7.4f}".format(structure.get_chemical_symbols()[j],j, structure.get_chemical_symbols()[j], j, 0, j, structure.charges[j]) )
  
      # O contribution at the H
      EF_OatH= structure.charges[iOw]/(4.*np.pi*e0*(OH)**2)*(vOH/OH)*eV_a2au;
      print (" {:>2}{:<4}{:>2}{:<4} dist: {:.6f}   chg{:<4} {:7.4f}".format(structure.get_chemical_symbols()[j],j, structure.get_chemical_symbols()[iOw], iOw, OH, iOw, structure.charges[iOw]) )
  
      # Contribution from the other H
      indexHother= np.setdiff1d(indexH_loc,j) # the index of the other H
      vHH= structure.get_positions()[j] - structure.get_positions()[indexHother[0]]
      HH=  structure.get_distance(indexHother[0],j,mic=True)
    
      EF_HatH= structure.charges[indexHother[0]]/(4.*np.pi*e0*(HH)**2)*(vHH/HH)*eV_a2au;
      if BADER_charges[indexHother[0]] > -0.01:
        print (" WARNING: H atom with idx:",indexHother[0],"has no associated bader charge")
      print (" {:>2}{:<4}{:>2}{:<4} dist: {:.6f}   chg{:<4} {:7.4f}".format(structure.get_chemical_symbols()[j],j, structure.get_chemical_symbols()[indexHother[0]], indexHother[0], HH, indexHother[0], structure.charges[indexHother[0]]) )
  
      EF_H= EF_H - EF_OatH - EF_HatH;
  
      print ("-"*50)
      print (" H{0:<4}  rOH= {1:9.6f}  EF@H_along_bisector= {2:3f} a.u.".format(j,OH,np.dot(EF_H,vMid)) ) 
      print ("                        EF@H_along_OH_bond=  {0:3f} a.u.\n".format(np.dot(EF_H,uOH)) ) 

#.#

#================================================================
def EF_H2O_DDEC6(structure_in, verbose=True):
#-# 
  from copy import deepcopy
  from ase.neighborlist import NeighborList
  from pymatgen.core import Lattice, Structure
  from pymatgen.analysis import ewald

  structure= deepcopy(structure_in)
  index_Ow= get_indexes_Ow(structure)

  if "DDEC6_charges" in structure.row.data.keys():
    DDEC6_charges= structure.row.data["DDEC6_charges"]
  else:
    print ("Missing DDEC6_charges")
    sys.exit()

  nl=NeighborList([1.2/2]*len(structure),skin=0, self_interaction=False,bothways=True); nl.update(structure)
  indexH= [i for i,x in enumerate(structure.get_chemical_symbols()) if x == 'H'] # H indexes

  # Run Ewald
  print ("Running Ewald ...")
  lattice= Lattice(structure.get_cell())
  coords= structure.get_scaled_positions()
  structure.charges= DDEC6_charges
  struct= Structure(lattice,structure.get_chemical_symbols(), coords, site_properties={"charge" : structure.charges})
  Ew= ewald.EwaldSummation(struct, compute_forces=True)
  print ("Done")
  
  for iOw in index_Ow:
    print (" Ow idx:", iOw)
    print ("="*50)
  
    vCENTER=(structure.get_cell()[0]+structure.get_cell()[1]+structure.get_cell()[2])/2. # get translation vector to put the O in the center of the box
    structure.translate(-structure.get_positions()[iOw] + vCENTER)    # center the O in the box
    structure.set_scaled_positions(structure.get_scaled_positions()) # put back everything in the box
  
    neighbors= nl.get_neighbors(iOw)[0]
    indexH_loc= np.intersect1d(neighbors,indexH)
    vOH1= structure.get_positions()[indexH_loc[0]] - structure.get_positions()[iOw]
    vOH2= structure.get_positions()[indexH_loc[1]] - structure.get_positions()[iOw]
    vMid= (vOH1/np.linalg.norm(vOH1) + vOH2/np.linalg.norm(vOH2)); vMid= vMid/np.linalg.norm(vMid)
  
    for j in indexH_loc:
      vOH= structure.get_positions()[j] - structure.get_positions()[iOw]
      uOH= vOH/np.linalg.norm(vOH)
      OH= structure.get_distance(iOw,j,mic=True)
  
      # Total EF at H position 
      EF_H= Ew.forces[j]/structure.charges[j]*eV_a2au;
      print (" {:>2}{:<4}{:>2}{:<4} dist: {:.6f}   chg{:<4} {:7.4f}".format(structure.get_chemical_symbols()[j],j, structure.get_chemical_symbols()[j], j, 0, j, structure.charges[j]) )
  
      # O contribution at the H
      EF_OatH= structure.charges[iOw]/(4.*np.pi*e0*(OH)**2)*(vOH/OH)*eV_a2au;
      print (" {:>2}{:<4}{:>2}{:<4} dist: {:.6f}   chg{:<4} {:7.4f}".format(structure.get_chemical_symbols()[j],j, structure.get_chemical_symbols()[iOw], iOw, OH, iOw, structure.charges[iOw]) )
  
      # Contribution from the other H
      indexHother= np.setdiff1d(indexH_loc,j) # the index of the other H
      vHH= structure.get_positions()[j] - structure.get_positions()[indexHother[0]]
      HH=  structure.get_distance(indexHother[0],j,mic=True)
    
      EF_HatH= structure.charges[indexHother[0]]/(4.*np.pi*e0*(HH)**2)*(vHH/HH)*eV_a2au;
      print (" {:>2}{:<4}{:>2}{:<4} dist: {:.6f}   chg{:<4} {:7.4f}".format(structure.get_chemical_symbols()[j],j, structure.get_chemical_symbols()[indexHother[0]], indexHother[0], OH, indexHother[0], structure.charges[indexHother[0]]) )
  
      EF_H= EF_H - EF_OatH - EF_HatH;
  
      print ("-"*50)
      print (" H{0:<4} rOH= {1:9.6f}  EF@H_along_bisector= {2:3f} a.u.".format(j,OH,np.dot(EF_H,vMid)) )
      print ("                       EF@H_along_OH_bond=  {0:3f} a.u.\n".format(np.dot(EF_H,uOH)) )
#.#

#================================================================
def EF_H2O_wannier(structure_in, core_charges, verbose=True):
#-# 
  from copy import deepcopy   
  from ase.neighborlist import NeighborList
  from ase import Atom, Atoms
  from pymatgen.core import Lattice, Structure
  from pymatgen.analysis import ewald

  structure= deepcopy(structure_in)
  index_Ow= get_indexes_Ow(structure)

  if "WANNIER_xyz" in structure.row.data.keys():
    nwannier= len(structure.row.data["WANNIER_xyz"])
    structure.extend(Atoms(symbols="W"*nwannier, positions= structure.row.data["WANNIER_xyz"], cell= structure.get_cell()))
  else:
    print ("ERROR: Missing WANNIER_xyz")
    sys.exit()
  
  nl=NeighborList([1.2/2]*len(structure),skin=0, self_interaction=False,bothways=True);  nl.update(structure)
  
  indexH= [i for i,x in enumerate(structure.get_chemical_symbols()) if x == 'H'] # H indexes
  indexW= [i for i,x in enumerate(structure.get_chemical_symbols()) if x == 'W'] # W indexes

  # Run Ewald
  print ("Running Ewald ...")
  lattice= Lattice(structure.get_cell())
  coords= structure.get_scaled_positions()
  structure.charges= [core_charges[x] for x in structure.get_chemical_symbols()]
  struct= Structure(lattice,structure.get_chemical_symbols(), coords, site_properties={"charge" : structure.charges})
  Ew= ewald.EwaldSummation(struct, compute_forces=True)
  print ("Done")
  
  
  for iOw in index_Ow:
    print (" Ow idx:", iOw)
    print ("="*50)
  
    vCENTER=(structure.get_cell()[0]+structure.get_cell()[1]+structure.get_cell()[2])/2. # get translation vector to put the O in the center of the box
    structure.translate(-structure.get_positions()[iOw] + vCENTER)    # center the O in the box
    structure.set_scaled_positions(structure.get_scaled_positions()) # put back everything in the box
  
    neighbors= nl.get_neighbors(iOw)[0]
    indexH_loc= np.intersect1d(neighbors,indexH)
    vOH1= structure.get_positions()[indexH_loc[0]] - structure.get_positions()[iOw]
    vOH2= structure.get_positions()[indexH_loc[1]] - structure.get_positions()[iOw]
    vMid= (vOH1/np.linalg.norm(vOH1) + vOH2/np.linalg.norm(vOH2)); vMid= vMid/np.linalg.norm(vMid)
  
    for j in indexH_loc:
      vOH= structure.get_positions()[j] - structure.get_positions()[iOw]
      uOH= vOH/np.linalg.norm(vOH)
      OH= structure.get_distance(iOw,j,mic=True)
  
      # Total EF at H position 
      EF_H= Ew.forces[j]/structure.charges[j]*eV_a2au;
      print (" {:>2}{:<4}{:>2}{:<4} dist: {:.6f}   chg{:<4} {:7.4f}".format(structure.get_chemical_symbols()[j],j, structure.get_chemical_symbols()[j], j, 0, j, structure.charges[j]) )
  
      # O contribution at the H
      EF_OatH= structure.charges[iOw]/(4.*np.pi*e0*(OH)**2)*(vOH/OH)*eV_a2au;
  
      # Contribution from the other H
      indexHother= np.setdiff1d(indexH_loc,j) # the index of the other H
      vHH= structure.get_positions()[j] - structure.get_positions()[indexHother[0]]
      HH=  structure.get_distance(indexHother[0],j,mic=True)
    
      EF_HatH= structure.charges[indexHother[0]]/(4.*np.pi*e0*(HH)**2)*(vHH/HH)*eV_a2au;
      print (" {:>2}{:<4}{:>2}{:<4} dist: {:.6f}   chg{:<4} {:7.4f}".format(structure.get_chemical_symbols()[j],j, structure.get_chemical_symbols()[indexHother[0]], indexHother[0], HH, indexHother[0], structure.charges[indexHother[0]]) )
  
      EF_H= EF_H - EF_OatH - EF_HatH;
  
      # Wannier centers
      indexW_loc= np.intersect1d(neighbors,indexW)
      for k in indexW_loc:
        vXH= structure.get_positions()[j] - structure.get_positions()[k]
        XH= structure.get_distance(j,k,mic=True)
        print (" {:>2}{:<4}{:>2}{:<4} dist: {:.6f}   chg{:<4} {:7.4f}".format(structure.get_chemical_symbols()[j],j, structure.get_chemical_symbols()[k], k, XH, k, structure.charges[k]) )
  
        EF_XatH= structure.charges[k]/(4.*np.pi*e0*(XH)**2)*(vXH/XH)*eV_a2au;
        EF_H= EF_H - EF_XatH;
  
      print ("-"*50)
      print (" H{0:<4}  rOH= {1:9.6f}  EF@H_along_bisector= {2:3f} a.u.".format(j,OH,np.dot(EF_H,vMid)) )
      print ("                        EF@H_along_OH_bond=  {0:3f} a.u.\n".format(np.dot(EF_H,uOH)) )
#.# 

#================================================================
def find_structure_symmetry(structure, symprec= 1.e-5, verbose= True):
#-#
  import spglib
  from operator import concat
  natoms= len(structure)

  dataset= spglib.get_symmetry_dataset(structure, symprec= symprec)
  equivalent_atoms= dataset['equivalent_atoms'].copy()
  spacegroup= dataset["international"]

  # "re-align" the equivalent labels
  atom_types= set(structure.get_atomic_numbers())
  for iat in atom_types:
    atom_list= [i for i,x in enumerate(structure.get_atomic_numbers()) if x == iat]
    first_index= equivalent_atoms[atom_list].min() -1
    equivalent_atoms[atom_list]= equivalent_atoms[atom_list] - first_index
  
  atom_labels= map(concat, structure.get_chemical_symbols(), map(str,equivalent_atoms.tolist()))
  
  if verbose:
    print ("Spacegroup:", spacegroup)
  return spacegroup, equivalent_atoms, list(atom_labels)
#.#

#================================================================
def Ewald_energy(structure_in, core_charges, verbose=True):
#-#
  from copy import deepcopy
  from pymatgen.core import Lattice, Structure
  from pymatgen.analysis import ewald
  
  structure= deepcopy(structure_in)
  lattice= Lattice(structure.get_cell())
  coords= structure.get_scaled_positions()
  structure.charges= [core_charges[x] for x in structure.get_chemical_symbols()]
  struct= Structure(lattice,structure.get_chemical_symbols(), coords, site_properties={"charge" : structure.charges})

  Ew= ewald.EwaldSummation(struct, compute_forces=True)
  #Ew= ewald.EwaldSummation(struct, eta=0.057887, acc_factor= 12,compute_forces=True)

  if verbose:
    print ("="*50)
    print ("Real space          energy:",Ew.real_space_energy,"[eV]")
    print ("Reciprocal space    energy:",Ew.reciprocal_space_energy + Ew.point_energy,"[eV]")
    print ("="*50)
    print ("Total electrostatic energy:",Ew.total_energy,"[eV]")
  
  return Ew
#.#

#================================================================
def write2res(structurei, fileout):
#-#
  cell= structure.get_cell()
  frac= structure.get_scaled_positions()
  xyz=  structure.get_positions()
  
  if "ATOM_labels" in structure.row.data.keys():
    labels= structure.row.data["ATOM_labels"]
  else:
    labels= structure.get_chemical_simbols()
  
  if "WANNIER_xyz" in structure.row.data.keys():
    WANNIER_xyz= structure.row.data["WANNIER_xyz"]
   
  f=open(fileout,"w") 
  f.write("single\n\nvectors\n")
  for i in range(3):
      f.write("{:16} {:16} {:16}\n".format(cell[i,0], cell[i,1], cell[i,2] ) )
  
  #print "fractional"
  #for i in range(len(structure)):
  #  print ("{:<7} {:16.10f} {:16.10f} {:16.10f}".format(labels[i], frac[i,0], frac[i,1], frac[i,2] ))
  
  f.write("\ncartesian\n")
  for i in range(len(structure)):
    f.write("{:<7} {:16.10f} {:16.10f} {:16.10f}\n".format(labels[i], xyz[i,0], xyz[i,1], xyz[i,2] ) )
  if "WANNIER_xyz" in structure.row.data.keys():
    for i in range(len(WANNIER_xyz)):
      f.write("{:<7} {:16.10f} {:16.10f} {:16.10f}\n".format("W", WANNIER_xyz[i,0], WANNIER_xyz[i,1], WANNIER_xyz[i,2] ) )
  f.close()
#.#



#================================================================
def Main():
#-#
  import argparse
  import json

  parser = argparse.ArgumentParser(description='ASE database tools')
  parser.add_argument('dbname', help='database name')
  parser.add_argument('-id', type= int, dest='id', metavar="", help='number of the row in the database')

  group= parser.add_mutually_exclusive_group(required=True)
  group.add_argument('--add-VASP-tructure',         action= "store_true", help= 'Add VASP structre to the database')
  group.add_argument('--add-VASP-bader'   ,         action= "store_true", help= 'Add VASP Bader calculation results to the database')
  group.add_argument('--add-VASP-wannier' ,         action= "store_true", help= 'Add VASP Wannier calculation results to the database')
  group.add_argument('--add-VASP-DDEC6'   ,         action= "store_true", help= 'Add VASP DDEC6 calculation results to the database')
  group.add_argument('--add-symmetry'     ,         action= "store_true", help= 'Find and add symmetry for the structure in the database')
  group.add_argument('--add-OH-PES'       ,         action= "store_true", help= 'Add OH stretch PES in the database')
  group.add_argument('--add-DVR-results'  ,         action= "store_true", help= 'Add DVR results for the sytructure in the database')

  group.add_argument('--prop-H2O-bader-dipoles'   ,     action= "store_true", help= 'Calculate H2O Bader dipoles for the structure in the database')
  group.add_argument('--prop-H2O-wannier-dipoles' ,     action= "store_true", help= 'Calculate H2O Wannier dipoles for the structure in the database')
  group.add_argument('--prop-H2O-DDEC6-dipoles'   ,     action= "store_true", help= 'Calculate H2O DDEC6 dipoles for the structure in the database')

  group.add_argument('--EF-H2O-bader'      ,     action= "store_true", help= 'Calculate Bader EF at H positions for H2O molecules for the structure in the database')
  group.add_argument('--EF-H2O-wannier'    ,     action= "store_true", help= 'Calculate Wannier EF at H positions for H2O molecules for the structure in the database')
  group.add_argument('--EF-H2O-DDEC6'      ,     action= "store_true", help= 'Calculate DDEC6 EF at H positions for H2O molecules for the structure in the database')

  parser.add_argument("-d","--directory"  ,         type=str, metavar= "", help= "Specify folder input")
  parser.add_argument("-f","--file"       ,         type=str, metavar= "", help= "Specify file   input")
  parser.add_argument("-l","--label"      ,         type=str, metavar= "", help= "Specify label  input")

  parser.add_argument("-i","--index"      ,         type=int, metavar= "", help= "Specify index  input")
  parser.add_argument("-s","--symprec"    ,         type=float, metavar= "", default=1.e-5,  help= "Specify symmetry precision input")
  
  parser.add_argument("-chg","--core_charges"    ,         type=json.loads , metavar= "", help= "Specify core_charges in json format.\n \'{ \"H\" : 1.0 , \"O\" : 6.0 }\'  ")
  
  parser.add_argument("-v","--verbose"    ,         action= "store_true", help= "Verbose output")

  args= parser.parse_args()
  if args.verbose:
    print(args)

  # Check for missing input ====
  if not args.add_VASP_tructure and args.id is None:
    parser.error("Argument requires complimentary -id input")
  
  if (args.add_VASP_tructure \
      or args.add_VASP_bader \
      or args.add_VASP_DDEC6 \
      or args.add_VASP_wannier) and args.directory is None:
    parser.error("Argument requires complimentary --folder input")
  
  if (args.add_OH_PES or args.add_DVR_results) and args.file is None:
    parser.error("Argument requires complimentary --file input")
  
  if (args.add_OH_PES) and args.label is None:
    parser.error("Argument requires complimentary --label input")
  
  if (args.add_DVR_results) and args.index is None:
    parser.error("Argument requires complimentary --index input")
  
  if (args.prop_H2O_wannier_dipoles \
      or args.prop_H2O_bader_dipoles \
      or args.EF-H2O-bader \
      or args.EF-H2O-wannier) and args.core_charges is None:
    parser.error("Argument requires complimentary --core_charges input")
   
  
  # Call requested routines =====
  if args.add_VASP_tructure:
    add_VASP_structure(args.dbname, args.folder)


  if args.add_VASP_bader:
    add_VASP_bader(args.dbname, args.id, args.folder)
  
  if args.add_VASP_wannier:
    add_VASP_wannier(args.dbname, args.id, args.folder)
  
  if args.add_VASP_DDEC6:
    add_VASP_DDEC6(args.dbname, args.id, args.folder)

  if args.add_OH_PES:
    add_OH_PES(args.dbname, args.id, args.file, args.label)

  if args.add_DVR_results:
    add_DVR_results(args.dbname, args.id, args.file, args.index)

  if args.add_symmetry:
    db_find_symmetry(args.dbname, args.id, args.symprec)


  if args.prop_H2O_wannier_dipoles:
    structure= db_read_row(args.dbname, args.id)
    prop_H2O_wannier_dipoles(structure, args.core_charges)
  
  if args.prop_H2O_bader_dipoles:
    structure= db_read_row(args.dbname, args.id)
    prop_H2O_bader_dipoles(structure, args.core_charges)

  if args.prop_H2O_DDEC6_dipoles:
    structure= db_read_row(args.dbname, args.id)
    prop_H2O_DDEC6_dipoles(structure)


  return args
#.#

if __name__ == "__main__":
  args= Main()  




