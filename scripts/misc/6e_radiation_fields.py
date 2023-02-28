import yt
import numpy as np
import json
from os.path import exists
import functions as func
import load_ins as load

yt.enable_plugins()
yt.enable_parallelism()

run_dir = load.run_dir
data_dir = load.data_dir
image_outputs = load.image_outputs

kph_fields = [('enzo', 'HI_kph'), ('enzo', 'HeI_kph'), ('enzo', 'HeII_kph')]
weights = ['density', 'ones', 'cell_mass']

res = load.res
width = load.width

limits = False
if limits == True:
    out_dir = "limits"
else:
    out_dir = "no_limits"

massive_prog = load.massive_prog

run_halos = load.run_halos

runs = load.runs

fns = func.data_list()

halo_rad = {}
for sto, ds in fns.piter(storage = halo_rad, dynamic=True):
    run_name = str(ds.directory).split('/')[-2]
    o = str(ds).split('output_')[-1]
    
    temp_dict = {}
    temp_dict[run_name] = {}
    temp_dict[run_name][o] = {}
    
    halo_pos, halo_rvir, halo_mass = func.get_halo(ds, run_name, o)

    sp = ds.sphere(halo_pos, halo_rvir)
    
    for field in kph_fields:
        temp_dict[run_name][o][field] = {}
        tot_val = sp.quantities.total_quantity(field).to('1/s').v.tolist()
        temp_dict[run_name][o][field]['total'] = tot_val
        
        img_name = f'{data_dir}projection_{field[1]}_{run_name}_{o}.npy'
        
        if exists(img_name) == False:
            prj = yt.ProjectionPlot(ds, 'x', field, center = halo_pos, weight_field = 'ones')
            prj.set_unit(field, '1/s')
            frb = prj.data_source.to_frb(width, res)
            with open(img_name, 'wb') as f:
                np.save(f, np.array(frb[field]))
            del prj
            del frb
        else:
            print(f'{img_name} created.')
        
        for weight in weights:
            value = sp.quantities.weighted_average_quantity(field, weight).to('1/s').v.tolist()
            temp_dict[run_name][o][field][weight] = value
            
            
    sto.result_id = f'{run_name}_{o}'
    sto.result = temp_dict

if yt.is_root():
    with open('halos_rad.json', 'w') as outfile:
        json.dump(halo_rad, outfile)