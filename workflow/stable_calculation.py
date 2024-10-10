from ase.db import connect
import numpy as np
import toml
import argparse
import sys
sys.path.append('/home/energy/mahpe/Structure_generation')
from chemical_values import*
from perqueue.constants import DYNAMICWIDTHGROUP_KEY


def main(run_path:str ='.',db_dir:str ='.',MD_run:bool = False,**kwargs):

    # Load the config file
    cfg = run_path +'/config.toml'
    with open(cfg, 'r') as f:
        params = toml.load(f)

    # If MD_run is true we will run MD sim on all all relaxed structures 
    # Else we will only run MD sim for the most stable configuration for each concentration
    if MD_run == True:     
        dmkey = len(params['MD'])-1 # should not include the resource key
        return_parameters = {DYNAMICWIDTHGROUP_KEY: dmkey, 'run_path': run_path}
        return True, return_parameters
    
    # Load the database
    root_dir = params['root_dir']
    ion = params['initial_generation']['ion_cif']
    db = connect(db_dir)

    # Create a dictionary with the lowest energy structure for each composition
    structure_database = {}

    # Loop over all structures in the database
    for row in db.select():

        # Get the atoms object
        atoms = row.toatoms()

        # Initialize the dictionary for this concentration if it does not exist
        if sum(atoms.symbols==ion) not in structure_database.keys():
            structure_database[sum(atoms.symbols==ion)] = {}

        # Get the chemical formula without anion
        formula_config = ''.join((atoms.get_chemical_symbols())).replace(ion,'')
        formula_config = formula_config.replace('O','')
        formula_config = formula_config.replace('P','')
        formula_config = formula_config.replace('Si','')
        formula_config = formula_config.replace('S','')
    
        # Add the structure to the database if it is not already there
        if formula_config not in structure_database[sum(atoms.symbols==ion)].keys():
            structure_database[sum(atoms.symbols==ion)][formula_config] = {'ID': row.id, 'atoms': atoms,'symbols':atoms.symbols,'energy':atoms.get_potential_energy()}
        else: # If it is there compare the energies and keep the lowest one
            if atoms.get_potential_energy() < structure_database[sum(atoms.symbols==ion)][formula_config]['energy']:
                structure_database[sum(atoms.symbols==ion)][formula_config] = {'ID': row.id, 'atoms': atoms,'symbols':atoms.symbols,'energy':atoms.get_potential_energy()}

    Na_tot = np.max([k for k in structure_database.keys()])

    # Print stable structures and test if they are the same as in the previous concentration
    structures_prev = list(structure_database[list(structure_database.keys())[-1]])
    structures_prev.sort()
    params['MD'] = {}
    for key in structure_database.keys():
        structures = list(structure_database[key])
        structures.sort()
        print('Ion concentration: ',key/Na_tot, 'Number of structures: ',len(structures),'   ',structures )
        assert structures == structures_prev
        for config_key in structure_database[key].keys():
            struc = structure_database[key][config_key]
            row_i = db.get(id=struc['ID'])
            params['MD'][row_i.name] = {'db_ids':row_i.ids,'M_ion': params['Relax'][row_i.name]['M_ion'], 'db_dir': db_dir }

    with open(cfg, "w") as toml_file:
        toml.dump(params, toml_file)
    
    # Return the dynamic width group key and the run path
    dmkey = len(params['MD'])-1 # should not include the resource key
    return_parameters = {DYNAMICWIDTHGROUP_KEY: dmkey, 'run_path': run_path}
    return True, return_parameters

if __name__ == "__main__":
    main()