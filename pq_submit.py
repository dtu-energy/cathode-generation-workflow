from pathlib import Path
import toml
from perqueue import PersistentQueue
from perqueue.task_classes.task import Task
from perqueue.task_classes.task_groups import CyclicalGroup, Workflow, DynamicWidthGroup, StaticWidthGroup
from perqueue.constants import DYNAMICWIDTHGROUP_KEY
root_path = Path('/home/energy/mahpe/Structure_generation/')
path = Path('/home/energy/mahpe/Structure_generation/Perqueue_structure_generation/')
path_workflow = Path('/home/energy/mahpe/Structure_generation/Perqueue_structure_generation/workflow/')

with open(root_path/'config.toml', 'r') as f:
    params = toml.load(f)
print(params)

params['initial_generation']['cfg_name'] = str(root_path/'config.toml')
params['initial_generation']['root_dir'] = '/home/energy/mahpe/Playground/solid_state_NEB'#'/home/scratch3/mahpe'
print( params['initial_generation'])
initial_gen = Task(path/'Initial_generation.py',  params['initial_generation'], '24:xeon24el8_test:30m')#'8:sm3090_devel:3h'
cheter_relax = Task(path_workflow/'cheater_relaxation.py', {}, '24:xeon24el8:2d')#'8:sm3090_devel:3h'
relax = Task(path_workflow/'relaxation.py', {}, '24:xeon24el8:2d')#'8:sm3090_devel:3h')
stable = Task(path_workflow/'stable_calculation.py', {}, '24:xeon24el8:2d')#'8:sm3090_devel:3h')
md_sim = Task(path_workflow/'md_simulation.py', {}, '48:xeon24el8:2d')#'8:sm3090_devel:3h')

relaxation_group = DynamicWidthGroup([cheter_relax,relax])
md_simulation = DynamicWidthGroup([md_sim])

wf = Workflow({initial_gen:[],relaxation_group:[initial_gen],stable:[relaxation_group],md_simulation:[stable]})

with PersistentQueue() as pq:
    pq.submit(wf)
