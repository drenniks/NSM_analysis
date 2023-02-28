#import yt
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
import json
import numpy as np
import numpy.ma as ma
import os
import os.path
from os.path import exists
import functions as func
import load_ins as load

#yt.enable_plugins()
#yt.enable_parallelism()

## Set these in load_ins.py
set_limits = load.set_limits
plot_stars = load.plot_stars

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 18}

matplotlib.rc('font', **font)

run_dir = load.run_dir
data_dir = load.data_dir
image_outputs = load.image_outputs

stars = load.stars
types = load.types

DD_data = func.DD_data
massive_prog = func.massive_prog
run_halos = func.run_halos

halo_quantities = func.halo_quantities
star_info = func.star_info

proj_fields = load.proj_fields

phase_fields = load.phase_fields

profile_fields = load.profile_fields

units = load.units

field_labels = load.field_labels

colorbars = load.colorbars

starcolors = load.starcolors

## Projection limits
zlim = load.zlim

## Phase plot limits
zlim_phase = load.zlim_phase
xlim_phase = load.xlim_phase
ylim_phase = load.ylim_phase

## Profile plot limits
xlim_profile = load.xlim_profile
ylim_profile = load.ylim_profile

fns = func.data_list_noyt()

for ds in fns:
    run_name = ds.split('../')[1].split('/')[0]
    o = ds.split('_')[-1]

    time = DD_data[run_name][o]['time']
    redshift = DD_data[run_name][o]['redshift']

    chosen_star = func.param_info[run_name]['id']
    chosen_star_living = func.star_info[run_name][o]['chosen_stars'][chosen_star]['particle_position']
    if chosen_star_living == -1:
        chosen_star = -1

    print(run_name)
    print(o)
    
    ## Halo info
    halo_pos, halo_rvir, halo_mass = func.get_halo_nods(run_name, o)

    print('Making projections')
    func.make_projections(run_name, o, halo_rvir, halo_pos, time, redshift, proj_fields, field_labels, zlim, colorbars, starcolors, set_limits = set_limits)
    
    ## Phase Plots
    print('Making phase plots')
    func.make_phase(run_name, o, time, redshift, phase_fields, field_labels, xlim_phase, ylim_phase, zlim_phase, set_limits = set_limits, chosen_star = chosen_star)
    
    ## Profile Plots
    print('Making profile plots')
    func.make_profile(run_name, o, time, redshift, profile_fields, field_labels, xlim=xlim_profile, ylim = ylim_profile, set_limits = set_limits, chosen_star = chosen_star, halo_pos = halo_pos)
