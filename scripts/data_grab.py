import yt
import ytree
import numpy as np
import matplotlib.pyplot as plt
import json
import os
import os.path
from os.path import exists
import functions as func
import load_ins as load

yt.enable_plugins()
yt.enable_parallelism()

run_dir = load.run_dir

stars = load.stars

runs = load.runs

## Halo quantities to calculate and their weights, pre- and post-NSM + some radiation fields
halo_fields = load.halo_fields
halo_weights = load.halo_weights
total_fields = load.total_fields

## Profile plot fields
profile_fields = load.profile_fields

## Phase plot fields
phase_fields = load.phase_fields

## Projection plot fields
proj_fields = load.proj_fields

## Star types and fields to get information for
star_types = load.star_types
star_fields = load.star_fields

## Units of all fields
units = load.units

## Resolution of the fixed resoution buffers and width of projection
res = load.res
width = load.width

DD_data = load.DD_data
param_info = load.param_info
massive_prog = load.massive_prog
run_halos = load.run_halos

fns = func.data_list()
quantities = {}

for sto, ds in fns.piter(storage=quantities, dynamic=True):
    run_name = str(ds.directory).split('/')[-2]
    o = str(ds).split('output_')[-1]
    print(run_name)
    print(o)
    
    temp_dict = {}
    temp_dict[o] = {}
    
    hds = yt.load(f'{run_dir}{run_name}/rockstar_halos/halos_DD{o}.0.bin')
    
    for s in stars:
        ds.add_particle_filter(s)
    ds.add_particle_filter('p3')
    
    ad = ds.all_data()
    halos = hds.all_data()
    
    halo_pos, halo_rvir, halo_mass = func.get_halo(ds, run_name, o)

    sp = ds.sphere(halo_pos, halo_rvir)
    reg = ds.sphere(halo_pos, (20, 'kpc'))
    
    ## Adding star info to dictionary
    particular_stars = [param_info[run_name]['id']] 

    temp_dict[o]['star_info'] = func.star_data(ds, ad, star_fields, star_types, particular_stars = particular_stars)
    
    ## Adding halo quantities info to dictionary

    temp_dict[o]['halo_quantities'] = func.get_halo_info(ds, sp, halo_fields, halo_weights, units)
    temp_dict[o]['halo_quantities']['total'] = func.get_total_info(ds, sp, total_fields, units)
    
    ## Profile plots
    print('Making profile plots for ' + run_name + ' output_' + o)
    func.profile_plot_data(ds, sp, total_fields, run_name, o, units)

    ## Phase plots
    print('Making phase plots for ' + run_name + ' output_' + o)
    func.phase_plot_data(ds, sp, phase_fields, units, run_name, o)

    ## Projections
    print('Making projection plots for ' + run_name + ' output_' + o)
    func.projection_plot_data(ds, reg, halo_pos, proj_fields, units, run_name, o)

    sto.result_id = run_name + '_' + o
    sto.result = temp_dict

    ds.index.clear_all_data()


if yt.is_root():

    ## Save raw dictionary
    with open('raw_data_grab.json', 'w') as outfile:
        json.dump(quantities, outfile)

    ## Let's reorganize the dictionaries slightly.
    halo_quantities_unorder = {}
    star_info_unorder = {}

    for run_name in runs:
        halo_quantities_unorder[run_name] = {}
        star_info_unorder[run_name] = {}
        
    for i, key in enumerate(quantities):
        key_split = key.split('_')
        run_name = f'{key_split[0]}_{key_split[1]}'
        o = key_split[-1]
        
        halos = quantities[key][o]['halo_quantities']
        stars = quantities[key][o]['star_info']
        
        halo_quantities_unorder[run_name][o] = halos
        star_info_unorder[run_name][o] = stars

    halo_quantities = {}
    star_info = {}
        
    for i, key in enumerate(halo_quantities_unorder):
        outputs = np.array(list(halo_quantities_unorder[key].keys()))
        ordered = np.sort(outputs)
        
        halo_quantities[key] = {}
        star_info[key] = {}
        for o in ordered: 
            halo_quantities[key][o] = halo_quantities_unorder[key][o]
                
            star_info[key][o] = star_info_unorder[key][o]

        with open('halo_quantities.json', 'w') as outfile:
            json.dump(halo_quantities, outfile)

        with open('star_info.json', 'w') as outfile:
            json.dump(star_info, outfile)
