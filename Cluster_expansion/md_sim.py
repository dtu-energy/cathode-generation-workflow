import numpy as np
import ase
from ase.visualize import view
from ase.md import Langevin
from ase.io.trajectory import Trajectory
from ase.md import MDLogger
from ase.io import read
from ase import units 
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.calculators.vasp import Vasp
import numpy as np 
import weakref
import shutil
import os
import toml
from ase.db import connect
import argparse
from ase.calculators.calculator import CalculationFailed
from chemical_values import*

def get_arguments(arg_list=None):
    parser = argparse.ArgumentParser(
        description="Ab-initio molecular dynamic", fromfile_prefix_chars="+"
    )
    parser.add_argument(
        "--id",
        type=int,
        help="Id of the structure in the database",
    )
    parser.add_argument(
        "--db_path",
        type=int,
        help="Path to the data base",
    )
    parser.add_argument(
        "--run_path",
        type=int,
        help="Path you want to save the results",
    )
    parser.add_argument(
        "--M_ion", 
        type=str, 
        help="Metalic ion for this particular structure",
    )
    parser.add_argument(
        "--oszicar",
        type=bool,
        default=True,
        help="Save all the oszicar file",
    )
    parser.add_argument(
        "--procar",
        type=bool,
        default=False,
        help="Save all the procar file",
    )
    parser.add_argument(
        "--chgcar",
        type=bool,
        default=False,
        help="Save all the chgcar file",
    )
    parser.add_argument(
        "--outcar",
        type=bool,
        default=False,
        help="Save all the outcar file",
    )
    parser.add_argument(
        "--max_unconverged",
        type=int,
        default=50,
        help="Maximuk number of unconverged steps in MD before it stops",
    )
    return parser.parse_args(arg_list)

def update_namespace(ns, d):
    for k, v in d.items():
        if not isinstance(v, dict):
            ns.__dict__[k] = v

def main():
    # Load parameters
    args = get_arguments()
    with open('config.toml', 'r') as f:
        params = toml.load(f)

    params_md = params['MD_params']

    # Import arguments
    db_id = int(args.id)
    M_ion = str(args.M_ion)

    # Load the database
    db_dir = args.db_name
    db = connect(db_dir)
    row = db.get(id = db_id)
    atom = row.toatoms()
    name = row.name

    # Remove the vaccancies
    X_indice = [a.index for a in atom if a.symbol == 'X']
    if len(X_indice) != 0:
        del atom[X_indice]

    # setting and creating the directory for the saved files
    relax_directory = os.path.join(args.run_path,f'md_{M_ion}')
    relaxsim_directory =os.path.join(relax_directory,str(name))
    if not os.path.exists(relaxsim_directory):
        os.makedirs(relaxsim_directory)
    

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

    n_Mion = np.sum(atom.symbols == M_ion)
    n_Na = np.sum(atom.symbols == 'Na')
    tot_magmom = get_tot_magmom_alluadite(M_ion,n_Mion,n_Na) 
    

    # Vasp calculator
    vasp_params = params['VASP']
    vasp_params['nupdown'] = tot_magmom
    vasp_params['nsw'] = 0
    vasp_params['istart'] = 1 
    
    print(f'{name} has nupdown {tot_magmom}')

    calc = Vasp(**vasp_params)
    ldau_luj = {'ldau_luj':{}}
    if type(M_ion)==str:
        ldau_luj['ldau_luj'][M_ion] = {'L': 2, 'U': get_U_value(M_ion),'J':0}
        print(f'{name} has L, U, J values: (2, {get_U_value(M_ion)}, 0)')
    else:
        for m in M_ion:
             ldau_luj['ldau_luj'][m] = {'L': 2, 'U': get_U_value(m),'J':0}
             print(f'{name} has L, U, J values: (2, {get_U_value(m)}, 0)')
    calc.set(**ldau_luj)
    calc.set(directory=relaxsim_directory)
    # Set th VASP calcualtor
    atom.set_calculator(calc)

    # Set the momenta corresponding to T=1000K
    T = params_md['temperature'] # The temperature of the MD simulation
    T0 = str(T)
    f = params_md['friction_term'] # Frictional term in the Langevin equation

    MaxwellBoltzmannDistribution(atom, temperature_K=T)

    mdsim_name_log= 'md_'+T0+'K_'+name+'.log'
    mdsim_name_traj= 'MD.traj'#'md_'+T0+'K_'+name+'.traj'
    
    md = Langevin(atom, params_md['time_step'] * units.fs,
              temperature_K=T,
              friction=f,
              logfile=os.path.join(relaxsim_directory, mdsim_name_log) )

    traj = Trajectory(os.path.join(relaxsim_directory , mdsim_name_traj),
	    "w",atom)#properties=["magmoms","energy"])
   
   # Set and attach logger to save MD log file
    md.attach(traj.write, interval=params_md['dump_step']) 
    

    # Set and attach MD_saver to save vasp output files each step in MD. It is also used to limit the MD simulation
    if vasp_params['nelm']:
        nelm = vasp_params['nelm']
    else:
        nelm = 60 # vasp default value
    saver = MD_Saver(md, calc,nelm=nelm, max_unconverged=args.max_unconverged,root_dir=relaxsim_directory,backup_name="", 
            oszicar=args.oszicar,procar=args.procar,chgcar=args.chgcar,outcar=args.outcar)
    md.attach(saver, interval=params_md['dump_step'])

    # Start MD simulation
    md.run(params_md['max_step']) # Number of steps we want the simulation to run for 

if __name__ == "__main__":
    main()