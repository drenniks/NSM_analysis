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
from mpl_toolkits.axes_grid1 import AxesGrid 

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 22}

matplotlib.rc('font', **font)

yt.enable_plugins()
yt.enable_parallelism()

data_dir = "../"

with open('DD_data.json') as f:
    DD_data = json.load(f)

min_mass = 20*Msun
max_mass = 30*Msun

outputs = np.array(list(DD_data.keys()))
final_output = outputs[-1]
final_string = 'DD' + final_output + '/output_' + final_output

redshifts = []
for i, j in enumerate(DD_data):
    redshifts.append(DD_data[j][1])
redshifts = np.array(redshifts)

## Load in most recent dataset to find biggest halo

ds = yt.load(data_dir + final_string)
hds = yt.load(data_dir + 'old_rockstar_halos/halos_DD' + final_output + '.0.bin')

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

tree = ytree.load(data_dir + 'old_rockstar_halos/trees/tree_0_0_0.dat')
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
    
    if m > massive_prog[key]['mass']:
        massive_prog[key]['mass'] = m.v.tolist()
        massive_prog[key]['halo_id'] = identifier.tolist()
        massive_prog[key]['position'] = pos.v.tolist()
        massive_prog[key]['rvir'] = rvir.v.tolist()

if yt.is_root():
    with open('massive_prog_info.json', 'w') as outfile:
        json.dump(massive_prog, outfile)
'''
### Loop over the branch to find potential Pop III choices 
### Looping over the most massive progenitor
dataset_series = []

for i, o in enumerate(massive_prog):
    dataset_series.append(data_dir + 'DD' + o + '/output_' + o)

fns = yt.DatasetSeries(dataset_series)

potential_p3 = {}

print('Finding potential P3 stars in branch.')
for sto, ds in fns.piter(storage=potential_p3, dynamic=True):
    o = str(ds).split('output_')[1]
        
    temp_dict = {}

    temp_dict['redshift'] = ds.current_redshift
    temp_dict['time'] = ds.current_time.to('Myr').v.tolist()
    
    stars = ['p3', 'p2', 'p3_living']
    for s in stars:
        ds.add_particle_filter(s)
            
    ad = ds.all_data()
    
    halo_id = massive_prog[o]['halo_id']
    halo_pos = ds.arr(massive_prog[o]['position'], 'unitary').to('kpc')
    halo_rvir = ds.quan(massive_prog[o]['rvir'], 'unitary').to('kpc')
    halo_mass = ds.quan(massive_prog[o]['mass'], 'Msun')
    
    sp = ds.sphere(halo_pos, halo_rvir)
        
    p3_indices = sp['p3', 'particle_index']
    p3_test = sp['p3', 'particle_mass'].to('Msun')
    p3_pos = sp['p3', 'particle_position'].to('kpc')
    p2_pos = sp['p2', 'particle_position'].to('kpc')
    
    dr = np.sqrt(((p3_pos - halo_pos)**2).sum(1))
    huge_number = 1e20
        
    ## Replace with the central halo 
    
    prj = yt.ProjectionPlot(ds, 'x', 'density', center = halo_pos, weight_field='density', width=(10, 'kpc'))
    prj.annotate_sphere([0,0], halo_rvir, coord_system='plot')
    for p3 in p3_pos:
        prj.annotate_marker(p3, plot_args={'color':'white'})
    for p2 in p3_pos:
        prj.annotate_marker(p2, plot_args={'color':'black'})
    prj.annotate_timestamp(corner='upper_left', redshift=True, draw_inset_box=True)
    prj.save(data_dir + 'NSM_analysis/images/' + o + '_massive_density.png')
    
    p3_fixed = []
    for m in p3_test:
        if m < 1e-10:
            p3_fixed.append(m*huge_number)
        else:
            p3_fixed.append(m)
    p3_fixed = np.array(p3_fixed)
        
    any_stars = np.argwhere((min_mass < p3_fixed) & (p3_fixed < max_mass))
    
    if len(any_stars) == 0:
        #print('No P3 stars in desired mass range in this halo.')
        sto.result_id = o
        sto.result = []
        continue
    else:
        p_mass = []
        p_index = []
        p_dist = []
        for n in any_stars:
            p_mass.append(p3_test[n][0].v.tolist())
            p_index.append(p3_indices[n][0].v.tolist())
            p_dist.append(dr[n][0].v.tolist())
        temp_dict['mass'] = p_mass
        temp_dict['index'] = p_index
        temp_dict['dist'] = p_dist
        
        sto.result_id = o
        sto.result = temp_dict

## Reording data
keys = potential_p3.keys()

new_potential_p3 = {}

new_keys = []
for i in keys:
    new_keys.append(i)
new_keys = np.sort(new_keys)
for i, j in enumerate(new_keys):
    new_potential_p3[j] = potential_p3[j]

if yt.is_root():    
    with open('potential_p3_info.json', 'w') as outfile:
        json.dump(new_potential_p3, outfile)

potential_p3 = new_potential_p3

## Find the unique star ids 

identity = [] 
earliest_output = 1000
to_be_deleted = []

for i, j in enumerate(potential_p3):
    try:
        star_identity = potential_p3[j]['index']
        if star_identity[0] not in identity:
            identity.append(star_identity[0])
        if int(j) < earliest_output:
            earliest_output = int(j)
    
    except:
        #print('No stars in dataset ' + j)
        to_be_deleted.append(j)

for i in to_be_deleted:
    del potential_p3[i]

fig = plt.figure(figsize=(10, 8))
for star in identity:
    dist = []
    z = []
    for i, j in enumerate(potential_p3):
        index = np.argwhere(np.array(potential_p3[j]['index']) == star)
        
        if len(index) == 0:
            continue
        else:
            dist.append(potential_p3[j]['dist'][index[0][0]])
            z.append(potential_p3[j]['redshift'])
    plt.plot(dist, z, label = str(star))
plt.legend()
plt.xlabel('Distance to center of halo [kpc]')
plt.ylabel('Redshift')
plt.tight_layout()

if yt.is_root():
    plt.savefig('images/potential_distances.png')

# Grab min and max values of the below fields
print('Grabbing min and max values.')
ds = yt.load(data_dir + final_string)
hds = yt.load(data_dir + 'rockstar_halos/halos_DD' + final_output + '.0.bin')

stars = ['p3', 'p2', 'p3_living']
for s in stars:
    ds.add_particle_filter(s)

ad = ds.all_data()
halos = hds.all_data()

stars = ['p3', 'p2', 'p3_living']
for s in stars:
    ds.add_particle_filter(s)

ad = ds.all_data()

p3_pos = ad['p3', 'particle_position'].to('kpc')
p3_id = ad['p3', 'particle_index']

width = (10, 'kpc')
direction = 'x'
fields = ['density', 'temperature', 'El_fraction', 'metallicity']

min_max = {}

for star in range(len(identity)):
    p3_place = np.argwhere(p3_id == identity[star])[0][0]
    
    min_max[str(identity[star])] = {}

    p3_position = p3_pos[p3_place].to('kpc')

    for field in fields:
        min_max[str(identity[star])][field] = {}

        prj = yt.ProjectionPlot(ds, direction, field, center = p3_position.to('kpc'), width = width, weight_field = 'density', origin='native')

        p = prj.plots[field]
        min_max[str(identity[star])][field]['min'] = p.cb.vmin
        min_max[str(identity[star])][field]['max'] = p.cb.vmax

# Saving output dict
if yt.is_root():
    with open('min_max_' + str(width[0]) + '.json', 'w') as outfile:
        json.dump(min_max, outfile)

# Making images for each potential choice
data_dumps = []
for i, j in enumerate(potential_p3):
    data_dumps.append(data_dir + 'DD' + j + '/output_' + j)

ts = yt.DatasetSeries(data_dumps)

for ds in ts.piter():
    print(ds)
    key = str(ds).split('_')[1]
    
    stars = ['p3', 'p2', 'p3_living']
    for s in stars:
        ds.add_particle_filter(s)

    ad = ds.all_data()
    
    prog_pos = ds.arr(massive_prog[key]['position'],'unitary').to('kpc')
    prog_rvir = ds.quan(massive_prog[key]['rvir'], 'unitary').to('kpc')
    
    sp = ds.sphere(prog_pos, prog_rvir)
    
    all_id = ad['all', 'particle_index']
    all_pos = ad['all', 'particle_position']
    all_type = ad['all', 'particle_type']
    p2_pos = ad['p2', 'particle_position']

    potential_ids = potential_p3[key]['index']
    for star_focus in potential_ids:
        p3_place = np.argwhere(all_id == star_focus)[0][0]
        p3_position = all_pos[p3_place]
        p3_type = all_type[p3_place]
        
        for field in fields:
            min_val = min_max[str(star_focus)][field]['min']
            max_val = min_max[str(star_focus)][field]['max']

            prj = yt.ProjectionPlot(ds, direction, field, center = prog_pos, width = width, weight_field = 'density')
            prj.annotate_marker(p3_position, plot_args={'color':'white'})
            for p2 in p2_pos:
                prj.annotate_marker(p2, plot_args={'color':'black'})
            prj.annotate_sphere([0,0], prog_rvir, coord_system='plot')
            prj.annotate_timestamp(corner='upper_left', redshift=True, draw_inset_box=True)
            prj.annotate_scale(corner='upper_right')
            prj.set_zlim(field, zmin=min_val, zmax = max_val)
            prj.annotate_title(str(p3_type))
            prj.save('images/potential_choices/' + str(ds) + '_' + field + '_' + str(width[0]) + '_' + str(star_focus) + '.png')
'''
if yt.is_root():
    import time
    executionTime = (time.time() - startTime)
    print('Execution time in seconds: ' + str(executionTime))
