import yt
import json
import numpy as np 

### After simulation completes, what outputs did the chosen star change types? 

current_type = 5

NSM_dict = {}
NSM_dict['p3_binary'] = {}
NSM_dict['p3_binary']['output'] = []
NSM_dict['p3_binary']['time'] = []

NSM_dict['ns_binary'] = {}
NSM_dict['ns_binary']['output'] = []
NSM_dict['ns_binary']['time'] = []

NSM_dict['bh'] = {}
NSM_dict['bh']['output'] = []
NSM_dict['bh']['time'] = []


stars = ['p3_living', 'p3_binary', 'ns_binary', 'bh', 'p2']
types = [5, 12, 13, 6, 7]

with open('DD_data.json') as f:
    DD_data = json.load(f)

outputs = np.array(list(DD_data.keys())) ### Need to append this to start when the restart happened
data_dumps = []
for i in outputs:
    data_dumps.append(data_dir + 'DD' + i + '/output_' + i

fns = yt.load(data_dumps)

for ds in fns.piter():
     
    for star in stars:
        ds.add_particle_filter(star)

    ad = ds.all_data()
    chosen_index = np.argwhere(ad['all', 'particle_index'] == PopIII_params.PopIII_NSMParticleID)[0][0]
    chosen_type = ad['all', 'particle_type'][chosen_index]
    
    if chosen_type == 6:
        if chosen_type != current_type:
            index = np.argwhere(types == chosen_type)[0][0]
            star_type = stars[index]

            NSM_dict[star_type]['output'].append(str(ds))
            NSM_dict[star_type]['time'].append(ds.current_time.to('Myr').v.tolist())
            current_type = chosen_type
            
            break
    else:
        if chosen_type != current_type:
            index = np.argwhere(types == chosen_type)[0][0]
            star_type = stars[index]

            NSM_dict[star_type]['output'].append(str(ds))
            NSM_dict[star_type]['time'].append(ds.current_time.to('Myr').v.tolist())
            current_type = chosen_type

with open('NSM_dict.json', 'w') as file:
    json.dump(NSM_dict, file)