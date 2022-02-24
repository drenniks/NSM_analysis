import yt
from yt.units import Msun
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import json
from mpl_toolkits.axes_grid1 import AxesGrid

yt.enable_plugins()

fraction = .2
fraction_increase = .1
iterations = 5

data_dir = "/storage/home/hhive1/dskinner6/data/SG64-2020/"

with open('DD_data.json') as f:
    DD_data = json.load(f)

outputs = np.array(list(DD_data.keys()))
final_output = outputs[-1]
final_string = 'DD' + final_output + '/output_' + final_output

ds = yt.load(data_dir + final_string)

stars = ['p3', 'p2', 'p3_living']
for s in stars:
    ds.add_particle_filter(s)

ad = ds.all_data()

## Replace with the central halo 
sp = ds.sphere('c', (5, 'kpc'))

## Iterate over larger radii if there isn't a PopIII star that is within 20-30 Msun:
min_mass = 20*Msun
max_mass = 30*Msun

potential_p3_mass = []
potential_p3_index = []

for i in range(iterations):
    
    sp_central = ds.sphere('c', sp.radius*fraction)
    
    p3_indices = sp_central['p3', 'particle_index']
    p3_test = sp_central['p3', 'particle_mass'].to('Msun')
    huge_number = 1e20

    p3_fixed = []
    for m in p3_test:
        if m < 1e-10:
            p3_fixed.append(m*huge_number)
        else:
            p3_fixed.append(m)
    p3_fixed = np.array(p3_fixed)

    ## Check to see if there are any P3 stars in this radius that fall in the mass range.
    any_stars = np.argwhere((min_mass < p3_fixed) & (p3_fixed < max_mass))

    if len(any_stars) == 0:
        fraction += fraction_increase
        continue
    else:
        for n in any_stars:
            potential_p3_mass.append(p3_test[n])
            potential_p3_index.append(p3_indices[n])
        
        print(i)
        break

p2_pos = ad['p2', 'particle_position'].to('kpc')   
p3_pos = ad['p3', 'particle_position'].to('kpc')
p3_ct = ad['p3', 'creation_time'].to('Myr')
p3_id = ad['p3', 'particle_index']
p3_m = ad['p3', 'particle_mass'].to('Msun')

intersect, p3_index, ind_b = np.intersect1d(p3_id, potential_p3_index, return_indices=True)

p3_position = p3_pos[p3_index]
p3_creation = p3_ct[p3_index]
p3_ident = p3_id[p3_index]
p3_masses = p3_m[p3_index]

times = []
outputs = []
for i, j in enumerate(DD_data):
    times.append(DD_data[j]['time'])
    outputs.append(j)
times = np.array(times)
outputs = np.array(outputs)

closest_outputs = []
for i in range(len(p3_position)):
    diff = times - p3_creation[i].v # want difference to be positive
    positives = [num for num in diff if num > 0]
    diff_index = np.argwhere(diff == np.min(positives))[0][0]
    closest_outputs.append(outputs[diff_index])    

f = open("potential_p3_info.txt", "w")
for i in range(len(p3_position)):
    f.write('------------------------------------------------------\n')
    f.write('P3 index = ' + str(p3_ident[i]) + '\n')
    f.write('P3 mass [Msun] = ' + str(p3_m[i]) + '\n')
    f.write('P3 creation time [Myr] = ' + str(p3_creation.v[i]) + '\n')
    f.write('P3 output creation = ' + str(outputs[diff_index][i]) + '\n')
    f.write('------------------------------------------------------\n')
f.close()

potential = {}
for i in range(len(p3_position)):
    potential[str(p3_ident[i].v)] = {}
    potential[str(p3_ident[i].v)]['creation_time'] = p3_creation[i].v.tolist()
    potential[str(p3_ident[i].v)]['mass'] = p3_masses[i].to('Msun').v.tolist()
    potential[str(p3_ident[i].v)]['position'] = p3_position[i].v.tolist()
    
# Saving potential P3 stars

with open('potential_P3_info.json', 'w') as outfile:
    json.dump(potential, outfile)


fields = [
    ("gas", "density"),
    ("gas", "temperature"),
    ("gas", "El_fraction"),
    ("gas", "metallicity")
]


# Grab min and max values of the below fields
width = (3, 'kpc')
direction = 'x'
fields = ['density', 'temperature', 'El_fraction', 'metallicity']
min_max = {}    
for star in range(len(p3_ident)):
    min_max[str(p3_ident[star].v)] = {}
    p3_position = p3_pos[p3_index[star]].to('kpc')

    print(p3_position)

    for field in fields:
        min_max[str(p3_ident[star].v)][field] = {}

        prj = yt.ProjectionPlot(ds, direction, field, center = p3_position.to('kpc'), width = width, weight_field = 'density', origin='native')
        prj.annotate_marker(p3_position, plot_args={'color':'white'})
        #prj.annotate_marker((0, 0, 0), plot_args={'color':'white'})
        prj.annotate_timestamp(corner='upper_left', redshift=True, draw_inset_box=True)
        prj.annotate_scale(corner='upper_right')
        prj.annotate_title(str(int(p3_ident[star].v)))
        prj.show()

        p = prj.plots[field]
        min_max[str(p3_ident[star].v)][field]['min'] = p.cb.vmin
        min_max[str(p3_ident[star].v)][field]['max'] = p.cb.vmax

# Saving output dict

with open('min_max_' + str(width[0]) + '.json', 'w') as outfile:
    json.dump(min_max, outfile)

outs = []
for i in outputs[diff_index::]:
    outs.append(data_dir + 'DD' + i + '/output_' + i)

### Only iterate over the datasets where p3 star exists

ts = yt.load(outs)

for ds in ts.piter():
    
    stars = ['p3', 'p2', 'p3_living']
    for s in stars:
        ds.add_particle_filter(s)

    ad = ds.all_data()

    ## Replace with the central halo 
    sp = ds.sphere('c', (5, 'kpc'))
    
    p3_id = ad['p3', 'particle_index']
    p3_pos = ad['p3', 'particle_position']
    
    for star in range(len(p3_ident)):
        p3_position = p3_pos[p3_index[star]].to('kpc')
        for field in fields:
            min_val = min_max[str(p3_ident[star].v)][field]['min']
            max_val = min_max[str(p3_ident[star].v)][field]['max']

            prj = yt.ProjectionPlot(ds, direction, field, center = p3_position.to('kpc'), width = width, weight_field = 'density', origin='native')
            prj.annotate_marker(p3_position, plot_args={'color':'white'})
            prj.annotate_timestamp(corner='upper_left', redshift=True, draw_inset_box=True)
            prj.annotate_scale(corner='upper_right')
            prj.annotate_title(str(int(p3_ident[star].v)))
            prj.set_zlim(field, zmin=min_val, zmax = max_val)

            prj.save('images/potential_choices/' + str(ds) + '_' + field + '_' + str(width[0]) + '_' + str(int(p3_ident[star].v)) + '.png')
