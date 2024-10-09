from ase.calculators.vasp import Vasp
from ase.io import read, write, Trajectory
from ase.db import connect
from shutil import copy
import os, subprocess
import numpy as np
import argparse
import json
import toml
from pathlib import Path
from ase.calculators.calculator import CalculationFailed
import sys
sys.path.append('/home/energy/mahpe/Structure_generation')
from chemical_values import*
from perqueue.constants import INDEX_KW

def get_arguments(arg_list=None):
    parser = argparse.ArgumentParser(
        description="General Active Learning", fromfile_prefix_chars="+"
    )
    parser.add_argument(
        "--db_dir",
        type=str,
        help="Path to the data base",
    )
    parser.add_argument(
        "--db_ids", 
        type=int, 
        help="Id with chemical formula in the database for this particular structure",
    )
    parser.add_argument(
        "--M_ion", 
        type=str, 
        help="Metalic ion or list for this particular structure",
    )
    parser.add_argument(
        "--Job_name",
        type=str,
        help="Name of the releaxation to run",
    )
    parser.add_argument(
        "--cfg",
        type=str,
        help="Directory of the config file",
    )
    return parser.parse_args(arg_list)

def update_namespace(ns, d):
    for k, v in d.items():
        if not isinstance(v, dict):
            ns.__dict__[k] = v

def main(run_path:str='.',**kwargs):
    # Load perqueue index
    idx, *_ =kwargs[INDEX_KW]

    # Load the config file
    cfg = run_path +'/config.toml'
    with open(cfg, 'r') as f:
        params = toml.load(f)
    
    # Find the specific parameters for this relaxation using perqueue index
    Job_name = list(params['Cheater_relax'].keys())[idx]
    args = get_arguments()
    update_namespace(args, params['Cheater_relax'][Job_name])
    print(Job_name,idx)
    # Import arguments
    M_ion = args.M_ion
    db_dir = args.db_dir
    name = Job_name
    # Setting the particlur structure for this relaxation
    db = connect(db_dir)
    row = db.get(ids = args.db_ids)
    atom = row.toatoms()
    
    # Remove the vaccancies
    X_indice = [a.index for a in atom if a.symbol == 'X']
    del atom[X_indice]

    # setting and creating the directory for the saved files
    root_dir = params['root_dir']
    relax_directory = f'{root_dir}/cheater_relax'
    relaxsim_directory =f'{relax_directory}/{name}'
    try:
        os.makedirs(relaxsim_directory)
    except:
        pass

    # Vasp calculator
    vasp_params = params['VASP']
    tot_magmom = 0
    for a in atom:
        tot_magmom += a.magmom # find the total magmom     
    
    vasp_params['nupdown'] = tot_magmom
    print(f'{name} has nupdown {tot_magmom}')
    calc = Vasp(directory=relaxsim_directory,**vasp_params)
    ldau_luj = {'ldau_luj':{}}
    if type(M_ion)==str:
        ldau_luj['ldau_luj'][M_ion] = {'L': 2, 'U': get_U_value(M_ion),'J':0}
        print(f'{name} has L, U, J values: (2, {get_U_value(M_ion)}, 0)')
    else:
        for m in M_ion:
             ldau_luj['ldau_luj'][m] = {'L': 2, 'U': get_U_value(m),'J':0}
             print(f'{name} has L, U, J values: (2, {get_U_value(m)}, 0)')
    calc.set(**ldau_luj)

    # Set th VASP calcualtor
    atom.set_calculator(calc)

    # Start the calculation for structure optimization.
    atom.get_potential_energy()
    
    # Check if the relaxation have reaxhed required accuracy
    with open(relaxsim_directory+'/OUTCAR') as file:
        # read all lines using readline()
        lines = file.readlines()
        try:
            lines.index(' reached required accuracy - stopping structural energy minimisation\n')
            var= True
        except:
            var=False
    #if not var:
    #    sys.exit(f'Relaxation did not converge. Fmax: {max(np.sqrt(np.sum(np.square(atom.get_forces()),axis=1)))}' )
    
    # Update database 
    db_new = connect(relax_directory+'/structures_relax.db')
    db_new.write(atom,name=name, ids=args.db_ids,ids_plot= str(args.db_ids).split('_')[0]+'_'+str(atom.symbols), relaxed=var,conc=row.conc,
                 data = {'Ga_convergence': row.data.Ga_convergence})
    os.remove(relaxsim_directory+'/WAVECAR')
    # Update config.toml
    with open(cfg, 'r') as f:
        params = toml.load(f)
    
    params['Relax'][name]['db_dir'] = relax_directory+'/structures_relax.db'
    
    with open(cfg, "w") as toml_file:
        toml.dump(params, toml_file)

    return_parameters = {'run_path':run_path}
    return True, return_parameters

if __name__ == "__main__":
    main()