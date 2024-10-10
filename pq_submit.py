from pathlib import Path
import toml
from perqueue import PersistentQueue
from perqueue.task_classes.task import Task
from perqueue.task_classes.task_groups import CyclicalGroup, Workflow, DynamicWidthGroup, StaticWidthGroup
from perqueue.constants import DYNAMICWIDTHGROUP_KEY

### Define relavant Paths ###
root_path = Path('.') # Path to the root directory, where the git repository is cloned must be changed to the correct path
path_workflow = root_path/'workflow/'

### Define resources ###
init_gen_res = '1:local:10m'
relax_res = '1:local:10m'
stable_res = '1:local:10m'
md_sim_res ='1:local:10m'


### Load the parameters ###
with open(root_path/'config.toml', 'r') as f:
    params = toml.load(f)

# Update the parameters
params['initial_generation']['cfg_name'] = str(root_path/'config.toml')
params['initial_generation']['root_dir'] = str(root_path)

### Define the workflow Tasks ###
initial_gen = Task(root_path/'Initial_generation.py',  params['initial_generation'], init_gen_res)
psudo_relax = Task(path_workflow/'psudo_relaxation.py', {},relax_res )
relax = Task(path_workflow/'relaxation.py', {}, relax_res)
stable = Task(path_workflow/'stable_calculation.py', {}, stable_res)
md_sim = Task(path_workflow/'md_simulation.py', {}, md_sim_res)

### Define the workflow groups ###
relaxation_group = DynamicWidthGroup([psudo_relax,relax])
md_simulation = DynamicWidthGroup([md_sim])

### Define the workflow ###
wf = Workflow({initial_gen:[],relaxation_group:[initial_gen],stable:[relaxation_group],md_simulation:[stable]})

### Submit the workflow ###
with PersistentQueue() as pq:
    pq.submit(wf)
