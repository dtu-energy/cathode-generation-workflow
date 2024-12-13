# cathode-generation-workflow

This workflow generates combinatorial all possible atomic configurations a specific battery cathode material and perform Denisty Functional Theory (DFT) calculation on them as well as perform Molecular dynamic (MD) simulations on them.

The workflow diagram is vizualized below, using Na ion cathode materials as an example:
![plot](workflow.jpg)

This workflow is utilized with PerQueue and one need to install PerQueue, MyQueue, ASE, toml, itertools and have a python enviorment with a python version of 3.11 or newer. This is also presented in the requirment.txt
To intall required packages:
```bash
pip install -r requirement.txt
```

Afterwards PerQueue needs to be initialized before running the workflow:
```bash
pq init
```

To run the workflow you need to:
- Change the input paramters in config.toml. All settings for the workflow is set in config.toml, including the cathode material cif file to consider
- Change the paths and resources in pq_submit.py. To run the workflow you need to set the paths correct. You also need to change the resources used to run the different tasks. To use PerQueue see (https://gitlab.com/asm-dtu/perqueue). Perqueue uses MyQueue (https://myqueue.readthedocs.io/) to submit jobs to either your local computer or a cluster. MyQueue used the following notation;  n_cores:partition_name:submission_time , to submit jobs. The code is set at the moment to use local computer resources but different HPC cluster settings can be added. A set of HPC configuration files are presented in (https://github.com/dtu-energy/Myqueue-for-HPC). 

When these changes are done you utilize the workflow using the command:
```
python pq_submit.py
```
As an example, use one of the four cif files provided. The NaMPO4_olivine/NaMPO4_maricite phasespace are the smallest and it is recommended to use either one of them.

Note: 
- Only VASP is used as the DFT calculator, but others can be added upon request
- Only four cation (Fe, Mn, Co and Ni) are considered as possible dopants. Others can be added upon request
- Only PO, SiO and SO are considered as anions in stable_calculation.py. Others can be added upon request

Aditionally, in 'Random_sampling/Random.ipynb' we present an example script to random sample cathode materials with substitutions disorder, where Na and metal ions can swap positions. For randdom sampling the 'Random_sampling/relax_ion.py' script is used to relax the structures and 'Random_sampling/md_sim.py' is used to perform Molecular dynamic simulations.

Note:
- This random sample example only works for the alluadite structure. To get the specific cif file one need to download it from ICSD. The ID is 243842. When downloaded rename the cif to 'Na2Fe2SO4.cif' and add it to the '/cif_file' folder.
