import sys
import itertools
import numpy as np
from string import digits
import weakref
import shutil
import os
from ase.db import connect
def get_magmom(M:str,redox:bool) -> int:
    """
    Get the magnetic moment of the metal ion. Following mapping is used:
        Fe -> 4 (Fe2+)
        Mn -> 5 (Mn2+)
        Ni -> 2 (Ni2+)
        Co -> 3 (Co2+)
        Ga -> 0
        Fe -> 5 (Fe3+)
        Mn -> 4 (Mn3+)
        Ni -> 3 (Ni3+)
        Co -> 4 (Co3+)

    Args:
        M (str): The metal ion
        redox (bool): If the metal ion is in the redox state
    
    Returns:
        int: The magnetic moment of the metal ion
    """
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

def get_U_value(M:str)-> float:
    """ 
    Get the U value for the metal ion
    Ref: https://docs.materialsproject.org/methodology/materials-methodology/calculation-details/gga+u-calculations/hubbard-u-values
    Following mapping is used:
        Fe -> 5.3
        Mn -> 3.9
        Ni -> 6.2
        Co -> 3.32
    
    Args:
        M (str): The metal ion

    Returns:
        float: The U value of the metal ion
    """
    if M=='Mn':
        U_val =3.9
    elif M=='Fe':
        U_val=5.3
    elif M=='Ni':
        U_val=6.2
    elif M=='Co':
        U_val=3.32
    else:
        raise ValueError(f'U value is not known for {M}')
    return U_val

def redox_sort(metal_ion:str)->int:
    """
    Sort the metal ion based on where the redox will happen first.
    Following redox order is used:
        Fe -> 0
        Mn -> 1
        Co -> 2
        Ni -> 3

    Args:
        metal_ion (str): The metal ion
    
    Returns:
        int: The index of the metal ion

    """
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

def redox_combination(M_ion_list_str:list,number_redox:int)-> list:
    """
    Sort the metal ion based on where the redox will happen first

    Args:
        M_ion_list_str (list): The list of metal ions
        number_redox (int): The number of redox reactions
    returns:
        list: The combination of the redox reactions
    """
    # Sort the M_ion after where the redox will happen
    sorted_list = M_ion_list_str.copy()
    sorted_list.sort(key=redox_sort)
    combi_list = []
    for i in range(number_redox):
        M_ion_redox = sorted_list[i]
        combi_list.append(list(np.argwhere(np.array(M_ion_list_str)==M_ion_redox).T[0]))
    return list(set([tuple(sorted(k)) for k in itertools.product(*combi_list ) if len(tuple(set(k))) == len(k)]))

def redox_psudo(metal_ion:str)-> str:
    """
    Give each of the metal ion their own redox atom
    Following mapping is used:
        Fe -> Ga
        Mn -> In
        Ni -> Ti
        Co -> Al

    Args:
        metal_ion (str): The metal ion

    Returns:
        str: The redox atom
    """
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

def remove_redox_psudo(metal_ion:str) -> str:
    """
    Remove the redox atom.
    Following mapping is used:
        Ga -> Fe
        In -> Mn
        Ti -> Ni
        Al -> Co

    Args:
        metal_ion (str): The metal ion
    
    Returns:
        str: The metal ion
    """
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

def redox_psudo_replacement(metal_ion:str)-> str:
    """
    Replace the redox atom with Ga for all metal ions

    Args:
        metal_ion (str): The metal ion
    
    Returns:
        str: The metal ion
    """

    # Replace the redox atom such that they are all Ga
    if metal_ion =='In':
        return 'Ga'
    elif metal_ion == 'Ti':
        return 'Ga'
    elif metal_ion =='Al':
        return 'Ga'
    else:
        return metal_ion

def create_name_atom(atom_formula:str)-> str:
    """
    Create a name for the atom based on the ase database formula

    Args:
        atom_formula (str): The formula of the atom
    
    Returns:
        str: The name of the atom
    """
    # create name for the atom based on ase db formula
    # Remove digits from name
    #remove_digits = str.maketrans('', '', digits)
    #a_name = atom_formula.translate(remove_digits)
    #Keep the digits
    a_name = atom_formula
    # Remove psudos
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
    """
    Create a name for the atom based on the ase database formula

    Args:
        database_row (ase.db.row): The row of the database  
        m_metal_list (list): The list of metal ions 

    Returns:
        str: The name of the atom

    """

    # create name for the atom based on ase db formula
    # Remove digits from name
    atom_formula = database_row.formula
    remove_digits = str.maketrans('', '', digits)
    a_name = atom_formula.translate(remove_digits)
    # Remove psudos
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