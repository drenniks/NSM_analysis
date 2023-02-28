import time
startTime = time.time()

import yt
from yt.units import Msun
from ytree_significant import *
import ytree 
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import json
import os
from mpl_toolkits.axes_grid1 import AxesGrid 

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 22}

matplotlib.rc('font', **font)

yt.enable_plugins()
yt.enable_parallelism()

data_dir = "../run_original/"

with open('DD_data_OG.json') as f:
    DD_data = json.load(f)

#min_mass = 20*Msun
#max_mass = 30*Msun

outputs = np.array(list(DD_data.keys()))
final_output = outputs[-1]

final_string = 'DD' + final_output + '/output_' + final_output

redshifts = []
for i, j in enumerate(DD_data):
    redshifts.append(DD_data[j][1])
redshifts = np.array(redshifts)

## Load in most recent dataset to find biggest halo

ds = yt.load(data_dir + final_string)
hds = yt.load(data_dir + 'rockstar_halos/halos_DD' + final_output + '.0.bin')

stars = ['p3', 'p2', 'p3_living']
for s in stars:
    ds.add_particle_filter(s)

ad = ds.all_data()
halos = hds.all_data()

## Biggest halo
print('Finding the biggest halo')
big_index = np.argwhere(halos['particle_mass'] == np.max(halos['particle_mass']))[0][0]
big_mass = halos['particle_mass'][big_index].to('Msun')
big_pos = halos['particle_position'][big_index].to('kpc')
big_rvir = halos['virial_radius'][big_index].to('kpc')
big_ident = halos['particle_identifier'][big_index]

tree = ytree.load(data_dir + 'rockstar_halos/trees/tree_0_0_0.dat')
#fn = tree.save_arbor()
#tree = ytree.load(fn)

## Halo in tree:
tree_index = np.argwhere(tree['halo_id'] == big_ident)[0][0]
my_tree = tree[tree_index]
tree_nodes = list(my_tree['tree'])

## Initializing a dictionary for the most massive progenitor
z = []
for i in tree_nodes:
    redshift = i['redshift']

    if redshift not in z:
        z.append(redshift)

massive_prog = {}
for i in z:
    diff = redshifts - float(i)
    minimum = np.min(np.abs(diff))
    diff_index = np.argwhere(np.abs(diff) == np.min(minimum))[0][0]
    
    massive_prog[str(outputs[diff_index])] = {}
    massive_prog[str(outputs[diff_index])]['halo_id'] = -1
    massive_prog[str(outputs[diff_index])]['mass'] = -1
    massive_prog[str(outputs[diff_index])]['position'] = -1
    massive_prog[str(outputs[diff_index])]['rvir'] = -1
    massive_prog[str(outputs[diff_index])]['T/|U|'] = -1

for i in tree_nodes:
    redshift = i['redshift']
    
    diff = redshifts - redshift
    minimum = np.min(np.abs(diff))
    diff_index = np.argwhere(np.abs(diff) == np.min(minimum))[0][0]
    key = outputs[diff_index]
    
    m = i['mass']
    identifier = i['halo_id']
    pos = i['position'].to('unitary')
    rvir = i['virial_radius'].to('unitary')
    vir_rat = i['T/|U|']
    
    if m > massive_prog[key]['mass']:
        massive_prog[key]['mass'] = m.v.tolist()
        massive_prog[key]['halo_id'] = identifier.tolist()
        massive_prog[key]['position'] = pos.v.tolist()
        massive_prog[key]['rvir'] = rvir.v.tolist()
        massive_prog[key]['T/|U|'] = vir_rat.tolist()

if yt.is_root():
    with open('massive_prog_info_incl_vir.json', 'w') as outfile:
        json.dump(massive_prog, outfile)

if yt.is_root():
    import time
    executionTime = (time.time() - startTime)
    print('Execution time in seconds: ' + str(executionTime))
