from ase.calculators.vasp import Vasp
from ase.calculators.calculator import kptdensity2monkhorstpack
from ase.db import connect
from clease.tools import update_db
from ase.io import read
import os,sys
import shutil
import numpy as np
import argparse
import toml
from chemical_values import*

# Setting the input 
parser = argparse.ArgumentParser()
parser.add_argument('--id', default=None, required=False, type=int, help="ID of the structure in the database")
parser.add_argument('--M_ion', default=None, required=False, type=str, help="Define the transition metal ion")
parser.add_argument('--db_name', default=None, required=False, type=str, help="Path to the database")
parser.add_argument('--run_path', default=None, required=False, type=str, help="Path to the run directory")
cli_args = parser.parse_args()

def get_input(var_name):
    if auto_input := getattr(cli_args, var_name, None):
        print("Auto input:", auto_input)
        return auto_input
    else:
        return input("Manual input: ")

id = int(get_input("id"))

# Load the atom structure
db_name = str(get_input("db_name"))
db = connect(db_name)
row = db.get(id=id)
atom =row.toatoms()

# randomly displacement of the atoms, to break atoms 
atom.rattle(0.01)

# Changing the tranisiton metal ion from the cif file to the new ion
M_ion_old = 'Fe'

M_ion = str(get_input("M_ion"))

M_ion_indice = [a.index for a in atom if a.symbol == M_ion_old]
for i in M_ion_indice:
    atom[i].symbol = M_ion

name = row.name.replace(M_ion_old,M_ion)

# Determine the total magnetic moment
n_Mion = np.sum(atom.symbols == M_ion)
n_Na = np.sum(atom.symbols == 'Na')

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

mag_tot = get_tot_magmom_alluadite(M_ion,n_Mion,n_Na) 

# Remove the vaccancies
X_indice = [a.index for a in atom if a.symbol == 'X']
del atom[X_indice]

# setting and creating the directory for the saved files
run_path = str(get_input("run_path"))
relax_directory = os.path.join(run_path,f'Na2Fe2SO4/relax_{M_ion}')
relaxsim_directory =os.path.join(relax_directory,name)

try:
    os.makedirs(relaxsim_directory)
except:
    pass

# Load the database
db_dir = os.path.join(run_path,f'Na2Fe2SO4/all_{M_ion}_DB')
if not os.path.exists(db_dir):
    os.makedirs(db_dir)

db_new = connect(db_dir)


#ionic relaxation
with open('cathode-generation-workflow/config.toml', 'r') as f:
    params = toml.load(f)
vasp_params = params['VASP']
vasp_params['nupdown'] = mag_tot
print(f'{name} has nupdown {mag_tot}')

print(vasp_params)

calc = Vasp(**vasp_params,
                ldau_luj = {M_ion: {'L': 2, 'U': get_U_value(M_ion),'J':0}},
                directory=relaxsim_directory)

for a in atom:
    if a.symbol == M_ion:
        a.magmom = get_magmom(a.symbol, redox=False)

atom.set_calculator(calc)
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


db_new.write(atom,name=name, converged=var)
#update_db(uid_initial=id, final_struct=atom, db_name=db_name)
   
