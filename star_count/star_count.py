import yt
import json

yt.enable_plugins()
#yt.enable_parallelism()

ds = yt.load('DD0282/DD0282')

stars = ['p3', 'p2', 'p3_living', 'p3_bh']
for s in stars:
    ds.add_particle_filter(s)

ad = ds.all_data()

p3_ct = ad['p3', 'creation_time'].to('Myr')
p3_id = ad['p3', 'particle_index']
p3_m = ad['p3', 'particle_mass'].to('Msun')

p3_living_ct = ad['p3_living', 'creation_time'].to('Myr')
p3_living_id = ad['p3_living', 'particle_index']
p3_living_m = ad['p3_living', 'particle_mass'].to('Msun')

p3_bh_ct = ad['p3_bh', 'creation_time'].to('Myr')
p3_bh_id = ad['p3_bh', 'particle_index']
p3_bh_m = ad['p3_bh', 'particle_mass'].to('Msun')

p2_ct = ad['p2', 'creation_time'].to('Myr')
p2_id = ad['p2', 'particle_index']
p2_m = ad['p2', 'particle_mass'].to('Msun')

star_info = {}
star_info['p3'] = {}
star_info['p3_living'] = {}
star_info['p3_bh'] = {}
star_info['p2'] = {}

for i in range(len(p3_id)):
    star_info['p3'][str(p3_id[i].v)] = {}
    star_info['p3'][str(p3_id[i].v)]['mass'] = p3_m[i].v.tolist()
    star_info['p3'][str(p3_id[i].v)]['ct'] = p3_ct[i].v.tolist()

for i in range(len(p3_living_id)):
    star_info['p3_living'][str(p3_living_id[i].v)] = {}
    star_info['p3_living'][str(p3_living_id[i].v)]['mass'] = p3_living_m[i].v.tolist()
    star_info['p3_living'][str(p3_living_id[i].v)]['ct'] = p3_living_ct[i].v.tolist()

for i in range(len(p3_bh_id)):
    star_info['p3_bh'][str(p3_bh_id[i].v)] = {}
    star_info['p3_bh'][str(p3_bh_id[i].v)]['mass'] = p3_bh_m[i].v.tolist()
    star_info['p3_bh'][str(p3_bh_id[i].v)]['ct'] = p3_bh_ct[i].v.tolist()

with open('star_info.json', 'w') as outfile:
    json.dump(star_info, outfile)

