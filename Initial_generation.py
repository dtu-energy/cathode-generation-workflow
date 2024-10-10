import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
#from tqdm import tqdm
import ase
from ase.io import read
from ase.io.trajectory import Trajectory
from ase.visualize import view
from ase.db import connect
import itertools
from ase.utils.structure_comparator import SymmetryEquivalenceCheck
from string import digits
import logging
import toml
from tqdm import tqdm
import argparse
import sys
import shutil
from distutils.dir_util import copy_tree
import time
import sys
from chemical_values import*
from typing import Tuple, Dict, Any
from perqueue.constants import DYNAMICWIDTHGROUP_KEY

def main(cif_file:str, M_ion_cif:str,ion_cif:str, symmertry_tol:int, full:bool, mix_without_full:bool,
          min_conc:bool, four_M_ion_mix:bool, M_ion_list:list, deformation_dist:float,
          Add_MD:bool,cfg_name:str = 'config.toml',root_dir:str = '.') -> Tuple[bool, Dict[str, Any]]:
    
    # Read the unit cell for the fully sodiated system
    unit_cell = read(cif_file)

    # Generate all possible structures based on the sodiation level
    # Set up nessery values
    ion_indice = [a.index for a in unit_cell if a.symbol == ion_cif]
    ion_tot = len(ion_indice)

    # Defining the main folder from the zif file and metal ions
    Folder_name = cif_file.split('/')[-1][:-4]
    if not full:
        Folder_name = Folder_name+'_mix_'
    else:
        Folder_name = Folder_name+'_'

    for m in M_ion_list:
        Folder_name = Folder_name+m
    Folder_dir = os.path.join(root_dir,Folder_name)
    try:
        os.makedirs(Folder_dir)
    except:
        raise ValueError(f'Folder {Folder_dir} already exists')
    print(Folder_dir)
    # Set up the database
    if os.path.exists(f'{Folder_dir}/{Folder_name}_initial.db'):
        os.remove(f'{Folder_dir}/{Folder_name}_initial.db')
    db = connect(f'{Folder_dir}/{Folder_name}_initial.db')
   
    # Create the mixture of ions
    # The idea is to first generate a fully M_ion structure and then change M_ion combinatoric with the other mixtures
    if not full:
        M_ion_mix_list = M_ion_list
    else:
        M_ion_mix_list = [M_ion_cif]
    
    # Find the M_ion indices and make the system full of M_ion
    M_indice = [a.index for a in unit_cell if a.symbol == M_ion_cif]
    
    # makes the full system for M_ion
    if M_ion_cif != M_ion_mix_list[0]:
        for ind in M_indice:
            unit_cell.symbols[ind] = M_ion_mix_list[0]
        M_ion = M_ion_mix_list[0]
    else:
        M_ion = M_ion_cif
    
    atom_mixture_list = [unit_cell.copy()] # add the fully M_ion structure
    
    # Create the mixture state 
    for M_ion_mix in M_ion_mix_list: # loop over all the mixture ions
        if M_ion_mix == M_ion: # We do not want to add the fully M_ion structure again
                continue
        # Copy the old atom list and loop over the old generated structures
        atom_mixture_list_old = atom_mixture_list.copy() 
        for atom_old in atom_mixture_list_old:
            M_new_indice = [a.index for a in atom_old if a.symbol == M_ion]

            for M_mix in range(1,len(M_new_indice)+1): # loop over the number of M_ion to change        
                Mix_add_list = list(itertools.combinations(M_new_indice, r=M_mix))
                for Mix_add in Mix_add_list: # loop over the combination 
                    atom_new = atom_old.copy()
                    for Mix_ind in Mix_add: # loop over the indices
                        atom_new.symbols[Mix_ind] = M_ion_mix        
                    atom_mixture_list.append(atom_new)
    print('Number of structures generated:',len(atom_mixture_list))
    print(atom_mixture_list)

    # Remove the fully structures if we only consider mixtures
    if not full:
        remove_ind = []
        count = 0    
        for i, atom_mix in enumerate(atom_mixture_list):
            # Find all M_ion in the system based on the old M_ion index
            M_ion_symbols = [atom_mix[index].symbol for index in M_indice]
            unique_symbols = np.unique(M_ion_symbols)

            # Remove fully M_ion structures with only one M_ion
            if mix_without_full and len(unique_symbols) == 1:
                remove_ind.append(i)
            
            # Remove structures with less than 4 M_ion
            if four_M_ion_mix and len(unique_symbols) == 2:
                remove_ind.append(i)
        # Remove the fully M_ion occupied structures
        remove_ind.reverse()
        for rm in remove_ind:
            atom_mixture_list.pop(rm)
        print('Number of structures generated without fully:',len(atom_mixture_list))
    
    # Loop over all sodiation levels for each of the mixture ion structures:
    for atom_mix in atom_mixture_list:
        print(atom_mix.symbols)
        for ion_rm in range(ion_tot+1):
            count = 0
            if ion_rm == 0: # Fully Sodiated
                db.write(atom_mix,ion_rm=ion_rm,data = {'Ga_convergence': {}})
                continue
            if ion_rm/ion_tot > min_conc:
                continue
            # Combination list of different index to be removed
            ion_rm_list = list(itertools.combinations(ion_indice, r=ion_rm))
            print(f'Number of Na removed: {ion_rm}')
            #print(f'Na index to be removed: {ion_rm_list}')
            # Loop over the different combination 
            for del_list in ion_rm_list:
                # Delete the Na atoms and set vaccancy
                atom_0 = atom_mix.copy()
                del_list = list(del_list) # need to sort the index before del atoms, to aviod index problems deleting M instead of Na
                del_list.sort(reverse = True)
                for del_idx in del_list:
                    atom_0.symbols[del_idx] = 'X'
        
                #  Find M indices and replace it by Ga
                Ga_add_list = redox_combination(list(atom_0.symbols[M_indice]),ion_rm)
                #print(Ga_add_list)
                for Ga_add in Ga_add_list:
                    atom_1 = atom_0.copy()
                    M_ion_to_Ga = {}
                    for i, Ga_ind in enumerate(Ga_add):
                        M_ion_to_Ga[i] = atom_1.symbols[M_indice[Ga_ind]]## IMPORTANT this convergence only works if the redox indicies are sorted (low->high) and when you rm a atom with ASE the atom order do not change 
                        atom_1.symbols[M_indice[Ga_ind]] = redox_psudo(M_ion_to_Ga[i])
                    # Write the total structure with removed Na and added redox (Ga)
                    db.write(atom_1, ion_rm=ion_rm,data = {'Ga_convergence': M_ion_to_Ga} )
                    count+=1
            print(f'Number of structures created: {count}')
    tot_len_db = len(db)
    print(f'Number of structures created in total: {tot_len_db}')

    # Symmertry comparision
    stol = symmertry_tol
    comp = SymmetryEquivalenceCheck(stol=stol)
    sym_id = []
    not_sym_id_tot = []
    not_sym_atoms_tot = []
    print('Symmertry check')
    # Loop over all sodium levels
    for ion_rm in tqdm(range(int(ion_tot*min_conc)+1)): 
        # Loop over all atoms for this sodium level:
        atoms = [row.toatoms() for row in db.select([('ion_rm','=',str(ion_rm))]) ]
        ids =[row.id for row in db.select([('ion_rm','=',str(ion_rm))]) ]
        not_sym_atoms = []
        not_sym_id = []
        for i in range(len(atoms)):
            # Set reference atoms and id
            id_0 = ids[i]
            atoms_0 = atoms[i]
            # The first atom is always not symmetrical
            if len(not_sym_atoms) == 0:
                not_sym_atoms.append(atoms_0)
                not_sym_id.append(id_0)
                continue
            
            # Compare id0 with all the other symmetrical structures
            c1 = comp.compare(atoms_0,not_sym_atoms)
            # if id0 is equal to some of the other non_symmertical structures remove it
            if c1:
                sym_id.append(id_0)
            else:
                # It can be that compare(a,b) != compare(b,a) so we need to check both ways
                c2 = [comp.compare(a, atoms_0) for a in not_sym_atoms]
                if any(c2):
                    sym_id.append(id_0)
                else:
                    not_sym_atoms.append(atoms_0)
                    not_sym_id.append(id_0)
        not_sym_id_tot += not_sym_id
        not_sym_atoms_tot += not_sym_atoms
        print(f'Number of symmetrical inequivarient structures: {len(not_sym_id_tot)}' )
    print(f'Total number of symmetrical inequivarient structures: {len(not_sym_id_tot)}' )
    # Remove the symmertical equivarient structures from database
    sym_del = [i for i in sym_id if i not in not_sym_id_tot]
    db.delete(sym_del)

    
    # Add the metal ions and create databases
    # If we only consider fully structures 
    if full:
        M_ion_list = M_ion_list
    else:
        M_ion_list = [1]
    db_list = []
    # Loop over all metal ions
    for M_ion_i in M_ion_list:
        # name the database for the structures
        if full:
            remove_digits = str.maketrans('', '', digits)
            atom_name = str(db.get(id=1).toatoms().symbols).translate(remove_digits)
            atom_name = atom_name.replace(M_ion,M_ion_i)
            if os.path.exists(f'{Folder_dir}/{atom_name}_psudo_structures.db'):
                os.remove(f'{Folder_dir}/{atom_name}_psudo_structures.db')
            db_new = connect(f'{Folder_dir}/{atom_name}_psudo_structures.db')
            db_list.append(f'{Folder_dir}/{atom_name}_psudo_structures.db')
            print(f'Create structures for: {atom_name}')
        else:
            atom_name = create_name_db(db.get(id=1),M_ion_mix_list)
            if os.path.exists(f'{Folder_dir}/{atom_name}_mix_psudo_structures.db'):
                os.remove(f'{Folder_dir}/{atom_name}_mix_psudo_structures.db')
            db_new = connect(f'{Folder_dir}/{atom_name}_mix_psudo_structures.db')
            db_list.append(f'{Folder_dir}/{atom_name}_mix_psudo_structures.db')

        # Loop over sym. equivarient structures
        index = 0
        for row in db.select():
            # name the database for the structures
            # Clean up after the symmerty comparison
            atom = row.toatoms().copy()
            for ind in M_indice: # loop over all metal ion and give them magmom without redox
                # change the atoms to the new metal ion
                if full and atom.symbols[ind] == M_ion: 
                    atom.symbols[ind] = M_ion_i
                M_ion_old = atom.symbols[ind]
                atom.symbols[ind] = redox_psudo_replacement(M_ion_old) # remove the "fake" psudos if any
                M_ion_new = atom.symbols[ind]
                atom[ind].magmom = get_magmom(M_ion_new,redox=False)
            # Assume that there are 6 oxygen close to the Metal ion
            Ga_indice = [a.index for a in atom if a.symbol == 'Ga']
            O_indice = [a.index for a in atom if a.symbol == 'O']
            for ind in Ga_indice: # Loop over the redox procces
                atom[ind].magmom = get_magmom('Ga',redox=True)  # this can be deleted since Ga always has magmom=0
        
                # Deformation of oxygen atoms
                pos =np.sqrt(np.sum((atom.positions[O_indice] - atom.positions[ind])**2,axis=1))
                min_ind =np.argpartition(pos, 6)[:6] # The closest 6 O atoms
                O_min_ind = np.array(O_indice)[min_ind]
                # Move the closest oxygen atoms closer to the metal ion
                for pos_ind in O_min_ind:
                    # only works with no negative values(there is one but it is approx 0)
                    argu1 = atom.positions[pos_ind]-atom.positions[ind]>0 #decrease distance
                    argu2 = atom.positions[pos_ind]-atom.positions[ind]<0 #incrase distance
                    atom.positions[pos_ind] = atom.positions[pos_ind] - deformation_dist*argu1 +deformation_dist*argu2

            if full:
                new_Ga_convergence = {}
                for key in row.data.Ga_convergence.keys():
                    new_Ga_convergence[key] = M_ion_i
                
                db_new.write(atom,ids = str(index)+'_'+str(atom.symbols),ids_old=str(row.id)+'_'+str(row.toatoms().symbols),name=create_name_atom(row.formula),
                         conc=(1-row.ion_rm)/ion_tot, data = {'Ga_convergence': new_Ga_convergence})
            else:
                db_new.write(atom,ids = str(index)+'_'+str(atom.symbols),ids_old=str(row.id)+'_'+str(row.toatoms().symbols),name=create_name_atom(row.formula),
                         conc=(1-row.ion_rm)/ion_tot, data = {'Ga_convergence': row.data.Ga_convergence})
            index += 1
    
    # Generate the workflow used for relaxation and perform MD
    # write toml file
    data = {}
    if Add_MD == True:
        list_scripts = ['Psudo_relax','Relax','MD']
    else:
        list_scripts = ['Psudo_relax','Relax']
    # set the toml for each of the script we want to run
    for script in list_scripts:
        data[script] = {}
        for i, db_name in enumerate(db_list):
            db_new = connect(db_name)
            for row in db_new.select():
                name = row.ids
                atom = row.toatoms

                if full:
                    name_split = name.split('_')
                    name_split.insert(1,M_ion_list[i])
                    name = '_'.join(name_split) 
                    if script == 'Psudo_relax':
                        data[script][name] = {'db_dir': db_name, 'db_ids':row.ids,'M_ion':M_ion_list[i]}
                    else:
                        data[script][name] = {'db_ids':row.ids,'M_ion':M_ion_list[i]}
                elif not full:
                    if script == 'Psudo_relax':
                        data[script][name] = {'db_dir': db_name, 'db_ids':row.ids,'db_ids_old':row.ids_old,'M_ion':M_ion_mix_list}
                    else:
                        data[script][name] = {'db_ids':row.ids,'M_ion':M_ion_mix_list}
    
    # Copy the config file and add the new data
    with open(cfg_name, 'r') as f:
        params = toml.load(f)
    output_file_name = f"{Folder_dir}/config.toml"

    shutil.copy('config.toml',output_file_name)    
    with open(output_file_name, "w") as toml_file:
        toml.dump({'root_dir':Folder_dir},toml_file)
        toml.dump(params, toml_file) # add the old params
        toml.dump(data, toml_file) # add the new params

    
    # Return the dynamic width group key and the run path
    dmkey = len(data['Psudo_relax'])-1 # should not include the resource key
    return_parameters = {DYNAMICWIDTHGROUP_KEY: dmkey, 'run_path': Folder_dir}
    return True , return_parameters

if __name__ == "__main__":
    main()
