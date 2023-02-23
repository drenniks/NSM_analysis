import time
start = time.time()

import yt
import ytree
import numpy as np
import matplotlib.pyplot as plt
import json


yt.enable_plugins()
yt.enable_parallelism()

run_dir = '../'

stars = ['p3_living', 'p3_binary', 'ns_binary', 'bh', 'p2', 'dm']
types = [5, 12, 13, 6, 7, 1]

with open('param_info.json') as f:
    param_info = json.load(f)
    
with open('run_DD_data.json') as f:
    DD_data = json.load(f)

with open('massive_prog_info.json') as f:
    massive_prog = json.load(f)

run_names = list(DD_data.keys())

restart_ds = '0037'

## See if we can match the run_original merger tree until the restart redshift. 

merger_tree_match = {}

for run_name in run_names:
    merger_tree_match[run_name] = {}

potential_ds = {}
matching_ds = {}
for i in run_names:
    potential_ds[i] = []
    matching_ds[i] = []

## Finds the closest dataset in each run for each original dataset

for i, ds in enumerate(massive_prog):
    match_redshift = DD_data['run_original'][ds]['redshift']
    for run_name in run_names:
        if run_name != 'run_original':
            
            match_position = massive_prog[ds]['position']
            match_mass = massive_prog[ds]['mass']

            target_redshifts = np.array([DD_data[run_name][key]['redshift'] for key in DD_data[run_name].keys()])
            target_outputs = np.array([key for key in DD_data[run_name].keys()])
            
            diff = target_redshifts - match_redshift
            minimum = np.min(np.abs(diff))
            diff_index = np.argwhere(np.abs(diff) == np.min(minimum))[0][0]
            
            target_output = target_outputs[diff_index]
            potential_ds[run_name].append([ds, target_output])

## It's possible the timing will be off. I'm going to find the datasets that closely match the original dataset

for i, run_name in enumerate(run_names):
    for j, ds in enumerate(DD_data[run_name]):
        ds_list = []
        for matches in potential_ds[run_name]:
            if ds == matches[1]:
                ds_list.append(matches)
                
        if len(ds_list) == 1:
            matching_ds[run_name].append(ds_list[0])
        
        elif len(ds_list) > 1:
            original_ds = [o[0] for o in ds_list]
            original_redshifts = [DD_data['run_original'][o[0]]['redshift'] for o in ds_list]
            matching_redshift = DD_data[run_name][ds_list[0][1]]['redshift']
            
            diff = np.array(original_redshifts) - matching_redshift
            minimum = np.min(np.abs(diff))
            diff_index = np.argwhere(np.abs(diff) == np.min(minimum))[0][0]
            
            target_output = original_ds[diff_index]
            
            matching_ds[run_name].append([target_output, ds_list[0][1]])
        
        else:
            continue

matching_halo_info = {}
for key in matching_ds.keys():
    matching_halo_info[key] = {}

## Match the merger tree at the beginning up until the restart data set as close as possible 

## Path to original run
b = ytree.load('../run_original/rockstar_halos/trees/tree_0_0_0.dat')

for i, run_name in enumerate(matching_ds):
    a = ytree.load(run_dir + run_name + '/rockstar_halos/trees/tree_0_0_0.dat')
    for outputs in matching_ds[run_name]:
        original_output = outputs[0]
        matching_output = outputs[1]

        output_images = run_dir + run_name + '/analysis/halo/'
        
        if int(matching_output) <= int(restart_ds):
            
            matching_redshift = DD_data[run_name][matching_output]['redshift']

            matching_halo_info[run_name][matching_output] = {}

            target_halo = massive_prog[original_output]['halo_id']
            target_mass = massive_prog[original_output]['mass']
            target_position = massive_prog[original_output]['position']
            target_rvir = massive_prog[original_output]['rvir']

            delta = 0.01
            criteria = "(tree['tree', 'redshift'] > " + str(matching_redshift - delta) + ") & (tree['tree', 'redshift'] < " + str(matching_redshift + delta) + ')'

            halos_original = list(b.select_halos(criteria))

            if len(halos_original) == 0:
                delta = 0.05
                criteria = "(tree['tree', 'redshift'] > " + str(matching_redshift - delta) + ") & (tree['tree', 'redshift'] < " + str(matching_redshift + delta) + ')'

                halos_original = list(b.select_halos(criteria))
                if len(halos_original) != 0:
                    print("You've got halos. Let's find the matching halo.")
                else:
                    print('Refine again!')

            else:
                print('Checking to find the matching halo.')

            ## Check to see if we found the matching halo in the original tree:
            original_positions = []
            original_masses = []
            original_ids = []
            for halo in halos_original:
                original_positions.append(halo['position'].to('unitary'))
                original_masses.append(halo['mass'])
                original_ids.append(halo['halo_id'])

            dr_position = np.sqrt(((np.array(original_positions) - target_position)**2).sum(1))
            dr_mass = np.abs(np.array(original_masses) - target_mass)

            min_pos_index = np.argwhere(dr_position == np.min(dr_position))[0][0]
            min_mass_index = np.argwhere(dr_mass == np.min(dr_mass))[0][0]

            if min_pos_index == min_mass_index:
                print('Found the matching halo! Halo ID #' + str(original_ids[min_pos_index]))
                original_redshift = halos_original[min_pos_index]['redshift']
                print(original_redshift)
            else:
                print(original_output)
                print('Cant find matching halo in original tree...')
                print('Must refine criteria.')


            halos = list(a.select_halos("(tree['tree', 'redshift'] == " + str(original_redshift) + ')'))

            if len(halos) == 0:
                halos = list(a.select_halos(criteria))

            positions = []
            masses = []
            ids = []
            rvirs = []
            vir_ratios = []
            for halo in halos:
                positions.append(halo['position'].to('unitary'))
                masses.append(halo['mass'])
                ids.append(halo['halo_id'])
                rvirs.append(halo['virial_radius'].to('unitary'))
                vir_ratios.append(halo['T/|U|'])

            dr_position = np.sqrt(((np.array(positions) - target_position)**2).sum(1))
            dr_mass = np.abs(np.array(masses) - target_mass)

            min_pos_index = np.argwhere(dr_position == np.min(dr_position))[0][0]
            min_mass_index = [i[0] for i in np.argwhere(dr_mass == np.min(dr_mass))]

            if min_pos_index in min_mass_index:
                print('Found the matching halo! Halo ID #' + str(ids[min_pos_index]))
                redshift = halos[min_pos_index]['redshift']

                #matching_halo_info[run_name][matching_output]['halo_id'] = ids[min_pos_index].v.tolist()
                matching_halo_info[run_name][matching_output]['position'] = positions[min_pos_index].v.tolist()
                matching_halo_info[run_name][matching_output]['mass'] = masses[min_pos_index].v.tolist()
                matching_halo_info[run_name][matching_output]['rvir'] = rvirs[min_pos_index].v.tolist()
                matching_halo_info[run_name][matching_output]['T/|U|'] = vir_ratios[min_pos_index]

            else:
                print('Cant find matching halo in original tree...')
                print('Choosing closest halo that is of the same order mass')
                order_mass = int("{:e}".format(target_mass).split('e+')[1])
                mass_criteria  = (np.array(masses) > 10**(order_mass)) & (np.array(masses) < 9*10**(order_mass))

                index = np.argwhere(dr_position == np.min(dr_position))[0][0]

                #matching_halo_info[run_name][matching_output]['halo_id'] = ids[index].v.tolist()
                matching_halo_info[run_name][matching_output]['position'] = positions[index].v.tolist()
                matching_halo_info[run_name][matching_output]['mass'] = masses[index].v.tolist()
                matching_halo_info[run_name][matching_output]['rvir'] = rvirs[index].v.tolist()
                matching_halo_info[run_name][matching_output]['T/|U|'] = vir_ratios[index]

            print('\n')
        else:
            continue

## After the restart output, we find the largest halo containing our star particle

particle_id = int(param_info['run_A']['id'])

## Identify position of the particle --> see what halo it's in 

print('Finding the halo that the particle is in.')
for i, run_name in enumerate(matching_halo_info):
    a = ytree.load('../' + run_name + '/rockstar_halos/trees/tree_0_0_0.dat')
    for j, o in enumerate(DD_data[run_name]):
        if int(o) > int(restart_ds):
            print(ds)
            ds = yt.load('../' + run_name + '/DD' + o + '/output_' + o)
            redshift = ds.current_redshift
            
            stars = ['p3', 'p2', 'p3_living']
            for s in stars:
                ds.add_particle_filter(s)
            
            ad = ds.all_data()
            
            particle_ids = ad['all', 'particle_index']
            particle_index = np.argwhere(particle_ids == particle_id)[0][0]
            particle_position = ad['all', 'particle_position'][particle_index].to('kpc')
            
            ## This works for now. Wonder about in the future... May have to find the closest redshift.
            delta = 0.01
            criteria = "(tree['tree', 'redshift'] > " + str(redshift - delta) + ") & (tree['tree', 'redshift'] < " + str(redshift + delta) + ')'

            halos_original = list(a.select_halos(criteria))
                
            halos_positions = []
            halos_rvirs = []
            halos_mass = []
            halos_vir_ratios = []
            
            for halo in halos_original:
                halos_positions.append(ds.arr(halo['position'].to('unitary'), 'unitary').to('kpc'))
                halos_rvirs.append(ds.quan(halo['virial_radius'].to('unitary'), 'unitary').to('kpc'))
                halos_mass.append(ds.quan(halo['mass'], 'Msun'))
                halos_vir_ratios.append(halo['T/|U|'])
            
            dr_position = np.sqrt(((halos_positions - particle_position)**2).sum(1))
            inside = dr_position < halos_rvirs
            
            halos_mass = np.array(halos_mass)[inside]
            halos_positions = np.array(halos_positions)[inside]
            halos_rvirs = np.array(halos_rvirs)[inside]
            halos_vir_ratios = np.array(halos_vir_ratios)[inside]
            
            massive_index = np.argwhere(halos_mass == np.max(halos_mass))[0][0]
            
            matching_halo_info[run_name][o] = {}
            matching_halo_info[run_name][o]['position'] = ds.arr(halos_positions[massive_index], 'kpc').to('unitary').v.tolist()
            matching_halo_info[run_name][o]['mass'] = ds.quan(halos_mass[massive_index], 'Msun').v.tolist()
            matching_halo_info[run_name][o]['rvir'] = ds.quan(halos_rvirs[massive_index], 'kpc').to('unitary').v.tolist()
            matching_halo_info[run_name][o]['T/|U|'] = halos_vir_ratios[massive_index]

class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)


file = open('run_halos_incl_vir.json', 'w')
file.write(json.dumps(matching_halo_info, cls=NpEncoder))

#with open('run_halos_incl_vir.json', 'w') as outfile:
#    json.dump(matching_halo_info, outfile)

finish = time.time()
total_time = finish - start
print('Execution time [s]: ' + str(total_time))