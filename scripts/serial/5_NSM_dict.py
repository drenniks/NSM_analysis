import yt
import json
import numpy as np 
import PopIII_params

yt.enable_plugins()

### After simulation completes, what outputs did the chosen star change types? 
### This code will also make images and movies.

data_dir = "../"
restart_dir = "../NSM_restart/"

current_type = 5
final_type = 1

NSM_dict = {}
NSM_dict['p3_binary'] = {}
NSM_dict['p3_binary']['output'] = []
NSM_dict['p3_binary']['time'] = []

NSM_dict['ns_binary'] = {}
NSM_dict['ns_binary']['output'] = []
NSM_dict['ns_binary']['time'] = []

NSM_dict['dm'] = {}
NSM_dict['dm']['output'] = []
NSM_dict['dm']['time'] = []


stars = ['p3_living', 'p3_binary', 'ns_binary', 'bh', 'p2', 'dm']
types = [5, 12, 13, 6, 7, 1]

width = (10, 'kpc')
direction = 'x'
fields = ['density', 'temperature', 'El_fraction', 'metallicity', 'NSMRProcess', 'SN_Colour']

movie = './movie.sh'
image_output = 'images/NSM_restart/'
movie_output = image_output + 'movies/'

with open('massive_prog_info.json') as f:
    massive_prog = json.load(f)

data_dumps = []
for i in os.listdir(restart_dir):
    if i.startswith('DD'):
        o = i.split('DD')[1]

        data_dumps.append(restart_dir + 'DD' + o + '/output_' + o)
        
data_dumps = np.sort(data_dumps).tolist()

fns = yt.DatasetSeries(data_dumps)

for ds in fns.piter():
     
    for star in stars:
        ds.add_particle_filter(star)

    ad = ds.all_data()

    o = str(ds).split('output_')[1]
    
    try:
        chosen_index = np.argwhere(ad['all', 'particle_index'] == PopIII_params.PopIII_NSMParticleID)[0][0]
        chosen_type = ad['all', 'particle_type'][chosen_index]
        chosen_pos = ad['all', 'particle_position'][chosen_index]
        p2_pos = ad['p2', 'particle_position'].to('kpc')

        type_name = stars[np.argwhere(types == chosen_type)[0][0]]
        
        halo_id = massive_prog[o]['halo_id']
        halo_pos = ds.arr(massive_prog[o]['position'], 'unitary').to('kpc')
        halo_rvir = ds.quan(massive_prog[o]['rvir'], 'unitary').to('kpc')
        halo_mass = ds.quan(massive_prog[o]['mass'], 'Msun')

        sp = ds.sphere(halo_pos, halo_rvir)

        for field in fields:
            prj = yt.ProjectionPlot(ds, 'x', field, center = halo_pos, weight_field='density', width=width)
            prj.annotate_sphere([0,0], halo_rvir, coord_system='plot')
            for p2 in p2_pos:
                prj.annotate_marker(p2, plot_args={'color':'black'})
            prj.annotate_marker(chosen_pos, plot_args={'color':'white'})
            prj.annotate_timestamp(corner='upper_left', redshift=True, draw_inset_box=True)
            prj.annotate_title(type_name)
            prj.save(data_dir + 'NSM_analysis/images/NSM_restart/' + o + '_' + str(field) + '.png')

        if chosen_type == 1:
            if chosen_type != current_type:
                ## This is when the NS_binary system merges and becomes a dark matter particle
                index = np.argwhere(types == chosen_type)[0][0]
                star_type = stars[index]

                NSM_dict[star_type]['output'].append(str(ds))
                NSM_dict[star_type]['time'].append(ds.current_time.to('Myr').v.tolist())
                current_type = chosen_type
            
                break
        else:
            if chosen_type != current_type:
                ## This will update current type
                index = np.argwhere(types == chosen_type)[0][0]
                star_type = stars[index]
                
                NSM_dict[star_type]['output'].append(str(ds))
                NSM_dict[star_type]['time'].append(ds.current_time.to('Myr').v.tolist())
                current_type = chosen_type

    except:
        print('PopIII_NSMParticleID not present in ' + str(ds))
with open('NSM_dict.json', 'w') as file:
    json.dump(NSM_dict, file)

for field in fields:
    image_list = []
    for i, j in enumerate(data_dumps):
        o = str(j).split('output_')[1]
        img_name =  o + '_' + str(field) + '.png'
        for img in os.listdir(image_output):
            if img == img_name:
                image_list.append(image_output + img)
    
    movie_images = " ".join(image_list)
    movie_name = field + '_' + str(width[0]) + '.mov'

    execute = movie + ' ' + movie_output + movie_name + ' 5 ' + movie_images
    os.system(execute)