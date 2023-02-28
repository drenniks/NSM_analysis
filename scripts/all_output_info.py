import yt
import numpy as np
import json
import os
import load_ins as load
import functions as func

yt.enable_plugins()
yt.enable_parallelism()

### Find outputs, time and z -- creates output.txt

run_dir = load.run_dir

runs = load.runs

DD_data = {}

fns = func.data_list()

for sto, ds in fns.piter(storage=DD_data, dynamic=True):
    run_name = str(ds.directory).split('/')[-2]
    output = str(ds).split('output_')[-1]
    
    time = ds.current_time.to('Myr').v.tolist()
    redshift = ds.current_redshift
    
    temp_dir = {}
    temp_dir[output] = {}
    temp_dir[output]['time'] = time
    temp_dir[output]['redshift'] = redshift 

    sto.result_id = run_name + '_' + output
    sto.result = temp_dir  
    #DD_data[run_name][output] = {}
    #DD_data[run_name][output]['time'] = time
    #DD_data[run_name][output]['redshift'] = redshift

if yt.is_root():
    ## Need to reorder the data
    new_DD_data = {}
    keys = DD_data.keys()
    for run in runs:
        new_DD_data[run] = {}

        new_keys = []
        for key in keys:
            if key.startswith(run):
                new_keys.append(key.split('_')[-1])
        
        new_keys = np.sort(new_keys)
        
        for i, j in enumerate(new_keys):
            new_DD_data[run][j] = {}
            new_DD_data[run][j]['time'] = DD_data[run + '_' + j][j]['time']
            new_DD_data[run][j]['redshift'] = DD_data[run+ '_' + j][j]['redshift']

    with open('run_DD_data.json', 'w') as f:
        json.dump(new_DD_data, f)
    
    param_info = {}
    for i, run_name in enumerate(runs):

        if run_name != 'run_original':
            #run_name = j.split('/')[-1]
            param_info[run_name] = {}
            param_file = run_dir + run_name + "/PopIII_params.py" 
            f = open(param_file, "r")
            for line in f.readlines():
                if line.startswith('PopIII_NSMParticleID'):
                    id = line.split('= ')[-1].split('\n')[0]
                    param_info[run_name]['id'] = id
                if line.startswith('PopIII_NSMExplosionEnergy'):
                    param_info[run_name]['energy'] = line.split('= ')[-1].split('\n')[0]
                if line.startswith('PopIII_NSMDelayTime'):
                    param_info[run_name]['delay_time'] = line.split('= ')[-1].split('\n')[0]
                if line.startswith('PopIII_NSMMetalMass'):
                    param_info[run_name]['metal_mass'] = line.split('= ')[-1].split('\n')[0] 
        else:
            param_info[run_name] = {}
            param_info[run_name]['id'] = id
            param_info[run_name]['energy'] = 0
            param_info[run_name]['delay_time'] = 0
            param_info[run_name]['metal_mass'] = 0

    with open('param_info.json', 'w') as outfile:
        json.dump(param_info, outfile)
