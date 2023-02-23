import yt
import json
import os
import numpy as np 
import functions as func
import load_ins as load

yt.enable_plugins()
yt.enable_parallelism()

### After simulation completes, what outputs did the chosen star change types? 
### This code will also make images and movies.

run_dir = load.run_dir

param_info = load.param_info

fns = func.data_list()

NSM_dict = {}


for sto, ds in fns.piter(storage = NSM_dict, dynamic=True):
    run_name = str(ds.directory).split('/')[-2]
    o = str(str(ds).split('output_')[1])

    ad = ds.all_data()

    temp_dir = {}

    particle_id = param_info[run_name]['id']
    
    if len(np.argwhere(ad['all', 'particle_index'] == particle_id)) != 0:
        chosen_index = np.argwhere(ad['all', 'particle_index'] == particle_id)[0][0]
        chosen_type = ad['all', 'particle_type'][chosen_index]
    
        chosen_pos = ad['all', 'particle_position'][chosen_index]
        temp_dir[o] = [chosen_type.v.tolist(), ds.current_time.to('Myr').v.tolist()] 
        sto.result = temp_dir
        sto.result_id = o

    else:

        temp_dir[o] = [0, 0]
        sto.result = temp_dir
        sto.result_id = o
        print(f'{particle_id} not present in DD{str(o)}! \n')
        continue

    
keys = NSM_dict.keys()
new_NSM_dict = {}

new_keys = []
for i in keys:
    if i != None:
        new_keys.append(i)

new_keys = np.sort(new_keys)
    
for i, j in enumerate(new_keys):
    new_NSM_dict[str(j)] = NSM_dict[j]

NSM_dict = new_NSM_dict

if yt.is_root():
    print(NSM_dict)
    with open('NSM_dict.json', 'w') as file:
        json.dump(NSM_dict, file)
