# cathode-generation-workflow

This workflow generates combinatorial all possible atomic configurations a specific battery cathode material and perform Denisty Functional Theory (DFT) calculation on them as well as perform Molecular dynamic (MD) simulations on them.

The workflow diagram is vizualized below, using Na ion cathode materials as an example:


This workflow is utilized with PerQueue and one need to install PerQueue, MyQueue, ASE, toml, itertools and have a python enviorment with a python version of 3.11 or newer.

To run the workflow you need to:
- Change the input paramters in config.toml. All settings for the workflow is set in config.toml, including the cathode material cif file to consider
- Change the paths and resources in pq_submit.py. To run the workflow you need to set the paths correct. You also need to change the resources used to run the different tasks. To use PerQueue see (https://gitlab.com/asm-dtu/perqueue). Perqueue uses MyQueue (https://myqueue.readthedocs.io/) to submit jobs to either your local computer or a cluster. MyQueue used the following notation;  n_cores:partition_name:submission_time , to submit jobs. The code is set at the moment to use local computer resources

When these changes are done you utilize the workflow using the command:
```
python pq_submit.py
```
As an example, use one of the four cif files provided. The NaMPO4_olivine/NaMPO4_maricite phasespace are the smallest and it is recommended to use either one of them.

Note: 
- Only VASP is used as the DFT calculator, but others can be added upon request
- Only PO, SiO and SO is considered as anions in stable_calculation.py. Others can be added upon request
