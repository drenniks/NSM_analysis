import yt
import numpy as np
import matplotlib.pyplot as plt
import json
import os
import os.path
from os.path import exists
from unyt import mass_hydrogen_cgs, eV, erg, Msun, s, Myr
from bisect import bisect
import load_ins as load
import functions as func
yt.enable_plugins()
yt.enable_parallelism()

limits = True
if limits == True:
    out_dir = "limits"
else:
    out_dir = "no_limits"

run_dir = load.run_dir
data_dir = load.data_dir
image_outputs = load.image_outputs

field = ('gas', 'gas_fraction')
weight_field = 'density'
zlimits = [2e-1, 9e-1]

dims = 'xyz'

vel_field = [("gas", "velocity_x"), ("gas", "velocity_y"), ("gas", "velocity_z")]

width = load.width

DD_data = load.DD_data
param_info = load.param_info
massive_prog = load.massive_prog
run_halos = load.run_halos

runs = load.runs

dataset_list = []
for i, run_name in enumerate(runs):
    if run_name == 'run_original':
        for j, o in enumerate(massive_prog):
            path = f'{image_outputs}{run_name}/{out_dir}/projection_gas_fraction_x_{o}.png'
            if os.path.isfile(path) != True:
                dataset_list.append(f'{run_dir}{run_name}/DD{o}/output_{o}')
                
    else:
        for j, o in enumerate(run_halos[run_name]):
            path = f'{image_outputs}{run_name}/{out_dir}/projection_gas_fraction_x_{o}.png'
            if os.path.isfile(path) != True:
                dataset_list.append(f'../{run_name}/DD{o}/output_{o}')

fns = yt.DatasetSeries(dataset_list)

for ds in fns.piter(dynamic=True):
    print(run_name)
    print(o)
    run_name = str(ds.directory).split('/')[-2]
    o = str(ds).split('output_')[-1]

    halo_pos, halo_rvir, halo_mass = func.get_halo(ds, run_name, o)

    sp = ds.sphere(halo_pos, halo_rvir)

    for dim in dims:
        if dim == 'x':
            q_field_x = ("gas", "velocity_y")
            q_field_y = ("gas", "velocity_z")
        elif dim == 'y':
            q_field_x = ("gas", "velocity_z")
            q_field_y = ("gas", "velocity_x")
        elif dim == 'z':
            q_field_x = ("gas", "velocity_x")
            q_field_y = ("gas", "velocity_y")
        try:
            prj = yt.ProjectionPlot(ds, dim, field, center = 'c', width = (10, 'kpc'), weight_field = weight_field)
            prj.annotate_quiver(q_field_x, q_field_y)
            prj.annotate_sphere(halo_pos, halo_rvir)
            prj.annotate_timestamp(corner="upper_left", redshift=True, draw_inset_box=True)
            prj.set_zlim(field, zlimits[0], zlimits[1])
            prj.save(f'{image_outputs}{run_name}/{out_dir}/projection_gas_fraction_{dim}_{o}.png')
        except Exception as e:
            print(f'{run_name} {o} {dim} except: {e}')
