import yt
import numpy as np
import json

yt.enable_plugins()
yt.enable_parallelism()

run_dir = '../'
data_dir = 'run_data/'
image_outputs = 'images/'

tiny_number = 1e-30

stars = ['p3_age', 'p3_binary', 'ns_binary', 'bh', 'p2', 'dm']
types = [5, 12, 13, 6, 7, 1]

with open('run_DD_data.json') as f:
    DD_data = json.load(f)
    
with open('star_info.json') as f:
    star_info = json.load(f)
    
with open('run_halos_incl_vir.json') as f:
    run_halos = json.load(f)

dataset_list = []

for i, run_name in enumerate(star_info):
    for j, o in enumerate(star_info[run_name]):
        if star_info[run_name][o]['star_info']['chosen_star']['id'] != -1:

            dataset_list.append(f'{run_dir}{run_name}/DD{o}/output_{o}')

        else:

            star_info[run_name][o]['star_info']['chosen_star']['velocity'] = -1

with open('star_info_with_vel.json', 'w') as f:
    json.dump(star_info, f)

fns = yt.DatasetSeries(dataset_list)

star_info_temp = {}

for sto, ds in fns.piter(storage = star_info_temp, dynamic = True):
    run_name = str(ds.directory).split('/')[-2]
    
    o = str(ds).split('output_')[-1]
    
    temp_dict = {}
    temp_dict[run_name] = {}
    temp_dict[run_name][o] = {}
    
    ident = star_info[run_name][o]['star_info']['chosen_star']['id']
    #ds = yt.load(f'{run_dir}{run_name}/DD{o}/output_{o}')
    ds.add_particle_filter('p3')
    ad = ds.all_data()
    
    p3_id = np.array(ad['all', 'particle_index'])
    where_p3 = np.argwhere(p3_id == ident)[0][0]

    p3_vel = ad['all', 'particle_velocity'][where_p3].v.tolist()
    
    temp_dict[run_name][o]['chosen_vel'] = p3_vel

    sto.result_id = run_name + '_' + o
    sto.result = temp_dict

if yt.is_root():
    with open('star_info_temp.json', 'w') as f:
        json.dump(star_info_temp, f)

    with open('star_info_with_vel.json', 'r') as f:
        star_info_vel = json.load(f)
    
    
    for i, run_name in enumerate(star_info_vel):
        for j, o in enumerate(star_info_vel[run_name]):
            vel_present = 'velocity' in star_info_vel[run_name][o]['star_info']['chosen_star'].keys()
            if vel_present == False:
                key_name = run_name + '_' + o
                star_info_vel[run_name][o]['star_info']['chosen_star']['velocity'] = star_info_temp[key_name][run_name][o]['chosen_vel']

    with open('star_info_vel.json', 'w') as f:
        json.dump(star_info_vel, f)