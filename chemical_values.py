import sys
import itertools
import numpy as np
from string import digits
import weakref
import shutil
import os

def get_magmom(M,redox):
    if M=='Fe':
        if redox == True:
            magmom = 5 # Fe3+
        else:
            magmom = 4 #Fe2+
    elif M=='Mn':
        if redox == True:
            magmom = 4 #Mn3+
        else:
            magmom = 5 #Mn2+
    elif M=='Ni':
        if redox == True:
            magmom = 3 #Ni3+ # was 1 but it was wrong due to the anions spin
        else:
            magmom = 2 #Ni2+
    elif M=='Co':
        if redox == True:
            magmom = 4 #Co3+
        else:
            magmom = 3 #Co2+
    elif M=='Ga':
        magmom= 0
    else:
        raise ValueError(f'Magmom is not known for {M}')
    return magmom

def get_U_value(M):
    if M=='Mn':
        U_val =3.9
    elif M=='Fe':
        U_val=5.3
    elif M=='Ni':
        U_val=6.2
    elif M=='Co':
        U_val=3.32
    else:
        sys.exist(f'U value is not known for {M}')
    return U_val

def get_tot_magmom(M_ion,n_M,n_Na):
    if M_ion=='Mn':
        mag_tot = 4*n_M +n_Na
    elif M_ion=='Fe':
        mag_tot = 5*n_M -n_Na
    elif M_ion=='Ni':
        mag_tot = 3*n_M -n_Na
    elif M_ion=='Co':
        mag_tot = 4*n_M -n_Na
    else:
        raise ValueError(f'Magmom is not known for {M_ion}')
    return mag_tot
def get_tot_magmom_silicate(M_ion,n_M,n_Na):
        if M_ion=='Mn':
            mag_tot = 4*n_M +n_Na-4
        elif M_ion=='Fe':
            mag_tot = 5*n_M -n_Na+4
        elif M_ion=='Ni':
            mag_tot = 3*n_M -n_Na+4
        elif M_ion=='Co':
            mag_tot = 4*n_M -n_Na+4
        else:
            sys.exit(f'Magmom is not known for {M_ion}')
        return mag_tot
def get_tot_magmom_alluadite(M_ion,n_M,n_Na):
        if M_ion=='Mn':
            mag_tot = 4*n_M +n_Na-3
        elif M_ion=='Fe':
            mag_tot = 5*n_M -n_Na+3
        elif M_ion=='Ni':
            mag_tot = 3*n_M -n_Na+3
        elif M_ion=='Co':
            mag_tot = 4*n_M -n_Na+3
        else:
            sys.exit(f'Magmom is not known for {M_ion}')
        return mag_tot

def redox_sort(metal_ion):
    # sort function based on where the redox will happen first
    # This sort knowlegde is based on their voltage profile
    if metal_ion =='Fe':
        return 0
    elif metal_ion == 'Mn':
        return 1
    elif metal_ion =='Co':
        return 2
    elif metal_ion =='Ni':
        return 3
    else:
        raise ValueError(f'Redox_sort is not known for {metal_ion}')

def redox_combination(M_ion_list_str,number_redox):
    # Sort the M_ion after where the redox will happen
    sorted_list = M_ion_list_str.copy()
    sorted_list.sort(key=redox_sort)
    combi_list = []
    for i in range(number_redox):
        M_ion_redox = sorted_list[i]
        combi_list.append(list(np.argwhere(np.array(M_ion_list_str)==M_ion_redox).T[0]))
    return list(set([tuple(sorted(k)) for k in itertools.product(*combi_list ) if len(tuple(set(k))) == len(k)]))

def redox_cheater(metal_ion):
    # Give each of the metal ion their own redox atom
    if metal_ion =='Fe':
        return 'Ga'
    elif metal_ion == 'Mn':
        return 'In'
    elif metal_ion =='Ni':
        return 'Ti'
    elif metal_ion =='Co':
        return 'Al'
    else:
        raise ValueError(f'Redox_sort is not known for {metal_ion}')

def remove_redox_cheater(metal_ion):
    # Give each of the metal ion their own redox atom
    if metal_ion =='Ga':
        return 'Fe'
    elif metal_ion == 'In':
        return 'Mn'
    elif metal_ion =='Ti':
        return 'Ni'
    elif metal_ion =='Al':
        return 'Co'
    else:
        raise ValueError(f'Redox_sort is not known for {metal_ion}')

def redox_cheater_replacement(metal_ion):
    # Replace the redox atom such that they are all Ga
    if metal_ion =='In':
        return 'Ga'
    elif metal_ion == 'Ti':
        return 'Ga'
    elif metal_ion =='Al':
        return 'Ga'
    else:
        return metal_ion

def create_name_atom(atom_formula):
    # create name for the atom based on ase db formula
    # Remove digits from name
    #remove_digits = str.maketrans('', '', digits)
    #a_name = atom_formula.translate(remove_digits)
    #Keep the digits
    a_name = atom_formula
    # Remove cheaters
    a_name = a_name.replace('Ga','Fe')
    a_name = a_name.replace('In','Mn')
    a_name = a_name.replace('Ti','Ni')
    a_name = a_name.replace('Al','Co')
    # Remove vacancy and add Na
    a_name = a_name.replace('X','')
    a_name = a_name.replace('Na','')
    a_name = 'Na'+a_name
    return a_name
def create_name_db(database_row,m_metal_list):
    # create name for the atom based on ase db formula
    # Remove digits from name
    atom_formula = database_row.formula
    remove_digits = str.maketrans('', '', digits)
    a_name = atom_formula.translate(remove_digits)
    # Remove cheaters
    a_name = a_name.replace('Ga','')
    a_name = a_name.replace('In','')
    a_name = a_name.replace('Ti','')
    a_name = a_name.replace('Al','')
    # Remove vacancy and add Na
    a_name = a_name.replace('X','')
    a_name = a_name.replace('Na','')
    
    for m_str in m_metal_list:
        a_name = m_str+a_name
    a_name = 'Na'+a_name
    return a_name

def get_Na_index(atoms):
    Na_index = [a.index for a in atoms if a.symbol == 'Na']
    return Na_index
def get_M_index(atoms):
    M_index = [a.index for a in atoms if a.symbol != 'O' and a.symbol != 'Na' and a.symbol != 'P' and a.symbol != 'S' and a.symbol != 'Si']
    return M_index
def redox_sort_func(x):
    if isinstance(x,tuple):
        x = x[0]
    if x == 'Fe':
        return 1
    elif x == 'Mn':
        return 2
    elif x == 'Co':
        return 3
    elif x == 'Ni':
        return 4
    elif x == 'Na':
        return 5
    else:
        return 10
def formula_sort_func(x):
    if isinstance(x,tuple):
        x = x[0]
    if x == 'Na':
        return 0
    if x == 'Fe':
        return 1
    if x == 'Mn':
        return 2
    if x == 'Co':
        return 3
    if x == 'Ni':
        return 4
    if x == 'O':
        return 10
    if x == 'P':
        return 5
    if x == 'S':
        return 7
    if x == 'Si':
        return 8
def sort(atoms, tags=None,key=None):
    if tags is None:
        tags = atoms.get_chemical_symbols()
    else:
        tags = list(tags)
    deco = sorted([(tag, i) for i, tag in enumerate(tags)],key=key)
    indices = [i for tag, i in deco]
    return atoms[indices], indices

class MD_Saver():
    def __init__(self, dyn, calc,nelm,max_unconverged=100, oszicar=True,procar=False,chgcar=False,outcar=False, root_dir="chgcars",scratch_dir="chgcars",backup_name="mart"):
        self.dyn = weakref.proxy(dyn)
        self.calc = weakref.proxy(calc)
        self.root_dir = root_dir
        self.scratch_dir = scratch_dir
        self.backup_name = backup_name
        self.oszicar = oszicar
        self.procar = procar
        self.chgcar = chgcar
        self.outcar = outcar
        self.counter = 0
        self.max_unconverged = max_unconverged
        self.nelm = nelm

    def __call__(self):

        step_num = self.dyn.nsteps
        
        # MD paths
        oszicar_path = self.calc.directory + "/OSZICAR"
        chgcar_path = self.calc.directory + "/CHGCAR"
        outcar_path = self.calc.directory + "/OUTCAR"
        procar_path = self.calc.directory + "/PROCAR"

        oszicar_dir = self.root_dir + "/OSZICAR_FILES" 
        chgcar_dir = self.scratch_dir + "/CHGCAR_FILES" 
        outcar_dir= self.scratch_dir + "/OUTCAR_FILES" 
        procar_dir = self.scratch_dir + "/PROCAR_FILES" 

        oszicar_dest = oszicar_dir +"/%04d.OSZICAR" % step_num
        chgcar_dest = chgcar_dir + "/%04d.CHGCAR" % step_num
        outcar_dest= outcar_dir + "/%04d.OUTCAR" % step_num
        procar_dest = procar_dir +"/%04d.PROCAR" % step_num
        

        # MD limiter
        # open Ozciar file and check if the SCF is finished
        with open(oszicar_path) as file:
            line = [[l] for l in file]
            nelm_scf = len(line)-2 # number of electronic steps, minus two to take account for header and footer
        
        if self.nelm <= nelm_scf:
            print('MD is not converging')
            self.counter +=1
        
        if self.counter == self.max_unconverged:
            raise ValueError("Too many unconverged MD steps!")

        # MD saver
        if self.oszicar:
            if not os.path.isdir(oszicar_dir):
                os.makedirs(oszicar_dir)
            shutil.copy(oszicar_path, oszicar_dest)
        if self.outcar:
            if not os.path.isdir(outcar_dir):
                os.makedirs(outcar_dir)
            shutil.copy(outcar_path, outcar_dest)
        if self.procar:
            if not os.path.isdir(procar_dir):
                os.makedirs(procar_dir)
            shutil.copy(procar_path, procar_dest)
        if self.chgcar:
            if not os.path.isdir(chgcar_dir):
                os.makedirs(chgcar_dir)
            shutil.copy(chgcar_path, chgcar_dest)

import re, os
from ase.io import read

def get_converged_structures(root_dir: str) -> list:
    """
    Get the converged structures from OUTCAR file
    Parameters
    ----------
    root_dir : str
        The root directory of the calculation

    Returns
    -------
    list of ase.Atoms objects 
    """
    # Get nelm from INCAR
    with open(os.path.join(root_dir,'INCAR'), 'r') as f:
        lines = f.readlines()
        for line in lines:
            if 'NELM = ' in line:
                nelm =int(line.split('=')[-1].strip())
                break
    # Get index of converged calculations
    print('NELM:',nelm)
    outcar_index = []
    count = 0
    with open(os.path.join(root_dir,'OSZICAR'), 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if 'F=' in line:
                header = lines[i-1].split(':')[0]
                scf = int(re.search(header+r':\s+(\d+)', lines[i-1]).group(1))
                if scf != nelm:
                    outcar_index.append(count)
                    count += 1
                else:
                    count += 1
    # Get converged structures from OUTCAR
    outcar = read(os.path.join(root_dir, 'OUTCAR'),format='vasp-out', index=':')
    structures = [a for i, a in enumerate(outcar) if i in outcar_index] 
    #print(len(structures), len(outcar_index), len(outcar))
    if len(structures) == 0:
        print('No converged structures found in {}'.format(root_dir))
    if len(outcar) == len(outcar_index):
        print('All structures are converged in {}'.format(root_dir))
    if len(outcar) != len(outcar_index):
        print(' {} structures are not converged in {}'.format(len(outcar)-len(outcar_index),root_dir))  
    return structures