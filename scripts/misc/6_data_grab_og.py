import yt
import ytree
import numpy as np
import matplotlib.pyplot as plt
import json
import os
import os.path
from os.path import exists


yt.enable_plugins()
yt.enable_parallelism()

data_dir = '../'
#restart_dir = data_dir + 'NSM_restart/'
image_output = 'images/NSM_restart/'
data_output = 'run_data/'

tiny_number = 1e-30

stars = ['p3_living', 'p3_binary', 'ns_binary', 'bh', 'p2', 'dm']
types = [5, 12, 13, 6, 7, 1]

## Profile plot fields for all datasets and for post-NSM
profile_fields = [('gas', 'metallicity'), ('gas', 'metallicity3'), ('gas', 'gas_fraction'), ('gas', 'temperature'), ('gas', 'electron_fraction'), ('gas', 'H2_fraction'), ('gas', 'density')]
profile_fields_postNSM = [("gas", "metallicity4")]

## Phase plot fields for all datasets and for post-NSM
phase_plot_fields = [[("gas", "density"), ("gas", "temperature")], [('index', 'radius'), ('gas', 'metallicity3')]]
phase_plot_colors = [[("gas", "cell_mass"), ("gas", "metallicity3")], [("gas", "cell_mass")]]
phase_plot_fields_postNSM = [[("gas", "metallicity3"), ("gas", "metallicity4")], [("gas", "density"), ("gas", "temperature")], [('index', 'radius'), ('gas', 'metallicity4')]]
phase_plot_colors_postNSM = [[("gas", "cell_mass")], [("gas", "metallicity4")], [("gas", "cell_mass")]]

## Projection plot fields for all datasets and for post-NSM
fields = [("gas", "density"), ('gas', 'temperature'), ('gas', 'metallicity'), ('gas', 'metallicity3')]
fields_postNSM = [('gas', 'metallicity4')]

## Halo quantities to calculate and their weights, pre- and post-NSM + some radiation fields
halo_quant_fields = [('gas', 'metallicity'), ('gas', 'metallicity3'), ('gas', 'gas_fraction'), ('gas', 'temperature'), ('gas', 'electron_fraction'), ('gas', 'H2_fraction'), ('gas', 'density')]
halo_quant_weights = [('gas', 'density'), ('gas', 'cell_mass')]
halo_quant_radfields = [('gas', 'J21_LW'), ('gas', 'J_Lyman')]
halo_quant_radweights = [('index', 'ones')]
halo_quant_fields_postNSM = [('gas', 'metallicity4')]

star_fields_p2 = [('p2', 'P2_metallicity_fraction'), ('p2', 'P3_metallicity_fraction'), ('p2', 'metallicity_fraction')]
star_fields_p3 = [('p3', 'P2_metallicity_fraction'), ('p3', 'P3_metallicity_fraction'), ('p3', 'metallicity_fraction')]
star_fields_postNSM_p2 = [('p2', 'NSM_metallicity_fraction')]
star_fields_postNSM_p3 = [('p3', 'NSM_metallicity_fraction')]

## Units of all fields
units = {("gas", "density"):'g/cm**3', 
         ("gas", "temperature"):'K', 
         ("gas", "cell_mass"):'g', 
         ('gas', 'metallicity'):'Zsun',
         ("gas", "metallicity3"):'Zsun', 
         ("gas", "metallicity4"):"Zr",
         ('gas', 'gas_fraction'):"dimensionless",
         ('gas', 'electron_fraction'):"dimensionless",
         ('gas', 'H2_fraction'):"dimensionless",
         ('gas', 'J21_LW'):"dimensionless", 
         ('gas', 'J_Lyman'):"erg/cm**2",
         ('index', 'radius'):"pc"
        }

## Resolution of the fixed resoution buffers and width of projection
res = 1000
width = (10, 'kpc')

with open('run_DD_data.json') as f:
    DD_data = json.load(f)

with open('param_info.json') as f:
    param_info = json.load(f)

with open('massive_prog_info.json') as f:
    massive_prog = json.load(f)

with open('run_halos_incl_vir.json') as f:
    run_halos = json.load(f)

run_dirs = []

for i in os.listdir(data_dir):
        if i.startswith('run_'):
                if os.path.isdir(data_dir + i) == True:
                        run_dirs.append(data_dir + i)
                        run_dirs = np.sort(run_dirs).tolist()

### First load in last dataset to find most massive progenitor for each run
last_dataset = []

for i, run_name in enumerate(DD_data):
    ## We already have the progenitor tree from the third script. Only need to do this for new runs. 
    if run_name != 'run_original':
        keys = list(DD_data[run_name].keys())
        last_dataset.append(data_dir + run_name + '/DD' + keys[-1] + '/output_' + keys[-1])

fns = yt.DatasetSeries(last_dataset)
'''
## Dictionary for the progenitor halo information for each NSM run
temp_run_info = {}

for sto, ds in fns.piter(storage = temp_run_info, dynamic=True):
    print(ds)
    run_name = str(ds.directory).split('/')[-2]
    
    o = str(ds).split('output_')[-1]
    
    ## Initializing a dictionary for the run
    temp_dict = {}
    
    redshifts = []
    for i, j in enumerate(DD_data[run_name]):
        redshifts.append(DD_data[run_name][j]['redshift'])
        temp_dict[j] = {}
        temp_dict[j]['halo_info'] = {}
        temp_dict[j]['halo_info']['halo_id'] = -1
        temp_dict[j]['halo_info']['mass'] = -1
        temp_dict[j]['halo_info']['position'] = -1
        temp_dict[j]['halo_info']['rvir'] = -1
        
    redshifts = np.array(redshifts)
    
    outputs = np.array(list(DD_data[run_name].keys()))
    
    hds = yt.load(data_dir + run_name + '/rockstar_halos/halos_DD' + o + '.0.bin')
    halos = hds.all_data()
    
    ## Biggest halo of the last dataset
    print('Finding the biggest halo')
    big_index = np.argwhere(halos['particle_mass'] == np.max(halos['particle_mass']))[0][0]
    big_mass = halos['particle_mass'][big_index].to('Msun')
    big_pos = halos['particle_position'][big_index].to('kpc')
    big_rvir = halos['virial_radius'][big_index].to('kpc')
    big_ident = halos['particle_identifier'][big_index]

    tree = ytree.load(data_dir + run_name + '/rockstar_halos/trees/tree_0_0_0.dat')

    ## Halo in tree:
    tree_index = np.argwhere(tree['halo_id'] == big_ident)[0][0]
    my_tree = tree[tree_index]
    tree_nodes = list(my_tree['tree'])
    
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

        if m > temp_dict[key]['halo_info']['mass']:
            temp_dict[key]['halo_info']['mass'] = m.v.tolist()
            temp_dict[key]['halo_info']['halo_id'] = identifier.tolist()
            temp_dict[key]['halo_info']['position'] = pos.v.tolist()
            temp_dict[key]['halo_info']['rvir'] = rvir.v.tolist()
        
    sto.result_id = run_name
    sto.result = temp_dict

if yt.is_root():
    with open('temp_run_info.json', 'w') as outfile:
        json.dump(temp_run_info, outfile)

         ### DO I EVEN USE TEMP_RUN?? DANIELLE
'''

## List of all datasets, excluding ones where a halo was not identified
dataset_list = []

for i, run_name in enumerate(run_halos):
    if run_name != 'run_original':
        for j, ds in enumerate(run_halos[run_name]):
            #        if temp_run_info[run_name][ds]['halo_info']['halo_id'] != -1:
            dataset_list.append(data_dir + run_name + '/DD' + ds + '/output_' + ds)
    else:
        for j, ds in enumerate(massive_prog):
            run_name = 'run_original'
            dataset_list.append(data_dir + run_name + '/DD' + ds + '/output_' + ds)

## Dictionary which will hold the halo quantities and star info
quantities = {}

fns = yt.DatasetSeries(dataset_list)

for sto, ds in fns.piter(storage=quantities, dynamic=True):
    run_name = str(ds.directory).split('/')[-2]
    o = str(ds).split('output_')[-1]
    print(run_name)
    print(o)
    
    temp_dict = {}
    temp_dict[o] = {}
    temp_dict[o]['star_info'] = {}
    temp_dict[o]['star_info']['chosen_star'] = {}
    temp_dict[o]['star_info']['P2'] = {}
    temp_dict[o]['star_info']['P3'] = {}
    
    temp_dict[o]['halo_quantities'] = {}
    
    hds = yt.load(data_dir + run_name + '/rockstar_halos/halos_DD' + o + '.0.bin')
    
    for s in stars:
        ds.add_particle_filter(s)
    ds.add_particle_filter('p3')
    
    ad = ds.all_data()
    halos = hds.all_data()
    
    ## Grab halo info from progenitor list found with ytree
    if run_name == 'run_original':
        halo_pos = ds.arr(massive_prog[o]['position'], 'unitary').to('kpc')
        halo_rvir = ds.quan(massive_prog[o]['rvir'], 'unitary').to('kpc')
        halo_mass = ds.quan(massive_prog[o]['mass'], 'Msun')
    else:
        halo_pos = ds.arr(run_halos[run_name][o]['position'], 'unitary').to('kpc')
        halo_rvir = ds.quan(run_halos[run_name][o]['rvir'], 'unitary').to('kpc')
        halo_mass = ds.quan(run_halos[run_name][o]['mass'], 'Msun')

    sp = ds.sphere(halo_pos, halo_rvir)
    reg = ds.sphere(halo_pos, (20, 'kpc'))
    
    chosen_index = np.argwhere(ad['all', 'particle_index'] == int(param_info[run_name]['id']))

    ## It's possible the chosen star doesn't exist yet. Check to see if it does.
    if len(chosen_index) != 0:
        print(ad)
        print(ad[('all', 'particle_type')])
        print('chosen_index = ', chosen_index)
        chosen_type = ad['all', 'particle_type'][chosen_index[0][0]]
        chosen_pos = ad['all', 'particle_position'][chosen_index[0][0]]
        chosen_mass = ad['all', 'particle_mass'][chosen_index[0][0]]
        chosen_ct = ad['all', 'creation_time'][chosen_index[0][0]]
        type_name = stars[np.argwhere(types == chosen_type)[0][0]]
    
        ## Adding chosen star info to dictionary
        temp_dict[o]['star_info']['chosen_star']['id'] = int(param_info[run_name]['id'])
        temp_dict[o]['star_info']['chosen_star']['position'] = chosen_pos.to('kpc').v.tolist()
        temp_dict[o]['star_info']['chosen_star']['mass'] = chosen_mass.to('Msun').v.tolist()
        temp_dict[o]['star_info']['chosen_star']['creation_time'] = chosen_ct.to('Myr').v.tolist()
        temp_dict[o]['star_info']['chosen_star']['type'] = chosen_type.v.tolist()
        
    else:
        temp_dict[o]['star_info']['chosen_star']['id'] = -1
        temp_dict[o]['star_info']['chosen_star']['position'] = -1
        temp_dict[o]['star_info']['chosen_star']['mass'] = -1
        temp_dict[o]['star_info']['chosen_star']['creation_time'] = -1
        temp_dict[o]['star_info']['chosen_star']['type'] = -1
        
    ## Adding P2 and P3 star info to dictionary ## check star_field DANIELLE
    if len(ad['p2', 'particle_index']) > 0:
        temp_dict[o]['star_info']['P2']['id'] = ad['p2', 'particle_index'].v.tolist()
        temp_dict[o]['star_info']['P2']['position'] = ad['p2', 'particle_position'].to('kpc').v.tolist()
        temp_dict[o]['star_info']['P2']['mass'] = ad['p2', 'particle_mass'].to('Msun').v.tolist()
        temp_dict[o]['star_info']['P2']['creation_time'] = ad['p2', 'creation_time'].to('Myr').v.tolist()
        for star_field in star_fields_p2:
            try:
                temp_dict[o]['star_info']['P2'][star_field[1]] = ad[star_field].to('Zsun').v.tolist()
            except Exception as e:
                print(e)

    
    temp_dict[o]['star_info']['P3']['id'] = ad['p3', 'particle_index'].v.tolist()
    temp_dict[o]['star_info']['P3']['position'] = ad['p3', 'particle_position'].to('kpc').v.tolist()
    temp_dict[o]['star_info']['P3']['mass'] = ad['p3', 'particle_mass'].to('Msun').v.tolist()
    temp_dict[o]['star_info']['P3']['creation_time'] = ad['p3', 'creation_time'].to('Myr').v.tolist()
    for star_field in star_fields_p3:
        try:
            temp_dict[o]['star_info']['P3'][star_field[1]] = ad[star_field].to('Zsun').v.tolist()
        except Exception as e:
            print(e)
    
    ## Adding halo quantity info to dictionary
    for hq_field in halo_quant_fields:
        temp_dict[o]['halo_quantities'][hq_field[1]] = {}
        for hq_weight in halo_quant_weights:
            temp_dict[o]['halo_quantities'][hq_field[1]][hq_weight[1]] = sp.quantities.weighted_average_quantity(hq_field, weight=hq_weight).to(units[hq_field]).v.tolist()
    
    for hq_field in halo_quant_radfields:
        temp_dict[o]['halo_quantities'][hq_field[1]] = {}
        for hq_weight in halo_quant_radweights:
            temp_dict[o]['halo_quantities'][hq_field[1]][hq_weight[1]] = sp.quantities.weighted_average_quantity(hq_field, weight=hq_weight).to(units[hq_field]).v.tolist()
    
    ## Profile plots
    print('Making profile plots for ' + run_name + ' output_' + o)

    make_prof_fields = []
    for l, prof_field in enumerate(profile_fields):
        prof_name = 'run_data/profile_' + prof_field[1] + '_' + run_name + '_' + o + '.npy'

        if exists(prof_name) == False:
            make_prof_fields.append(prof_field)
        else:
            print('Profile plot ' + prof_field[1] + ' already made for ' + run_name + ' ' + o)
    
    if len(make_prof_fields) != 0:
        try:
            plot = yt.ProfilePlot(sp, "radius", make_prof_fields, weight_field='cell_mass')
            plot.set_unit("radius", "kpc")
            profile = plot.profiles[0]

            with open('run_data/profile_x_' + run_name + '_' + o + '.npy', 'wb') as f:
                np.save(f, np.array(profile.x))

            for l, prof_field in enumerate(make_prof_fields):
                plot.set_unit(prof_field, units[prof_field])
                with open('run_data/profile_' + prof_field[1] + '_' + run_name + '_' + o + '.npy', 'wb') as f:
                    np.save(f, profile[prof_field])
            del plot
            del profile
        except Exception as e:
            print(run_name + ' DD' + o + ': ' + str(e))
    else:
        print('Profile plots complete for ' + run_name + ' ' + o)

    ## Phase plots
    print('Making phase plots for ' + run_name + ' output_' + o)
    for k, phase_field in enumerate(phase_plot_fields):
        for m, phase_color in enumerate(phase_plot_colors[k]):
            x = phase_field[0]
            y = phase_field[1]
            z = phase_color

            H_name = 'run_data/phase_H_' + x[1] + '_' + y[1] + '_' + z[1] + '_'+ run_name + '_' + o + '.npy'
            yedges_name = 'run_data/phase_yedges_' + x[1] + '_' + y[1] + '_' + z[1] + '_'+ run_name + '_' + o + '.npy'
            xedges_name = 'run_data/phase_xedges_' + x[1] + '_' + y[1] + '_' + z[1] + '_'+ run_name + '_' + o + '.npy'

            if ((exists(H_name) == False) and (exists(yedges_name) == False) and (exists(xedges_name) == False)):
                x_data = sp[x].to(units[x])
                y_data = sp[y].to(units[y])
                z_data = sp[z].to(units[z])

                ## Make sure that there are no zeros in the data. Set any zeros to a tiny number.
                x_zeroes = np.argwhere(np.array(x_data) == 0.0)
                if len(x_zeroes) != 0:
                    for i in x_zeroes:
                        x_data[i] = ds.quan(tiny_number, units[x])

                y_zeroes = np.argwhere(np.array(y_data) == 0.0)
                if len(y_zeroes) != 0:
                    for i in y_zeroes:
                        y_data[i] = ds.quan(tiny_number, units[y])
                        
                z_zeroes = np.argwhere(np.array(z_data) == 0.0)
                if len(z_zeroes) != 0:
                    for i in z_zeroes:
                        z_data[i] = ds.quan(tiny_number, units[z])

                x_bins = np.logspace(np.log10(np.min(x_data)), np.log10(np.max(x_data)), 128)
                y_bins = np.logspace(np.log10(np.min(y_data)), np.log10(np.max(y_data)), 128)

                H, yedges, xedges = np.histogram2d(x_data, y_data, [x_bins, y_bins], weights = z_data)

                with open('run_data/phase_H_' + x[1] + '_' + y[1] + '_' + z[1] + '_'+ run_name + '_' + o + '.npy', 'wb') as f:
                    np.save(f, H)
                with open('run_data/phase_yedges_' + x[1] + '_' + y[1] + '_' + z[1] + '_'+ run_name + '_' + o + '.npy', 'wb') as f:
                    np.save(f, yedges)
                with open('run_data/phase_xedges_' + x[1] + '_' + y[1] + '_' + z[1] + '_'+ run_name + '_' + o + '.npy', 'wb') as f:
                    np.save(f, xedges)

                del x
                del y
                del z
                del x_data
                del y_data
                del z_data
                del x_bins
                del y_bins
                del H
                del yedges
                del xedges
            else:
                print(x[1] + '_' + y[1] + '_' + z[1] + '_'+ run_name + '_' + o + ' created.')

    ## Projections
    print('Making projection plots for ' + run_name + ' output_' + o)
    for field in fields:
        img_name = 'run_data/projection_' + field[1] + '_' + run_name + '_' + o +'.npy'
        if exists(img_name) == False:
            prj = yt.ProjectionPlot(ds, 'x', field, center = halo_pos, weight_field='density', width=width, data_source = reg)
            frb = prj.data_source.to_frb(width, res)
            with open(img_name, 'wb') as f:
                np.save(f, np.array(frb[field]))
            del prj
            del frb
        else:
            print(img_name + ' created.')


    ## Check to see if the NSMRProcess field is present -- this indicates the NSM model has been turned on
    field_present = ("enzo", "NSMRProcess") in ds.derived_field_list
    print(run_name + ' DD' + o + ' field_present = ' + str(field_present))
    if field_present == True:
        def _metallicity4(field, data):
            return data['NSMRProcess'] / data['Density']
        ds.add_field(name=('gas', 'metallicity4'), units='Zr', display_name='NSM R-Process Metallicity', validators=[yt.ValidateDataField(('enzo', 'NSMRProcess'))], sampling_type='cell', function=_metallicity4)
        
        #print(ds.derived_field_list)
        ## Adding in post-NSM fields to dictionary

        for star_field in star_fields_postNSM_p2:
            temp_dict[o]['star_info']['P2'][star_field[1]] = ad[star_field].to('Zr').v.tolist()
        for star_field in star_fields_postNSM_p3:
            temp_dict[o]['star_info']['P3'][star_field[1]] = ad[star_field].to('Zr').v.tolist()

        for hq_field in halo_quant_fields_postNSM:
            temp_dict[o]['halo_quantities'][hq_field[1]] = {}
            for hq_weight in halo_quant_weights:
                temp_dict[o]['halo_quantities'][hq_field[1]][hq_weight[1]] = sp.quantities.weighted_average_quantity(hq_field, weight=hq_weight).to(units[hq_field]).v.tolist()
        
        ## Total R-process mass
        temp_dict[o]['halo_quantities']['total_NSMRProcess'] = {}
        temp_dict[o]['halo_quantities']['total_met4'] = {}
        
        temp_dict[o]['halo_quantities']['total_NSMRProcess']['none'] = sp.quantities.total_quantity('NSMRProcess').to('g/cm**3').v.tolist()
        temp_dict[o]['halo_quantities']['total_met4']['none'] = sp.quantities.total_quantity('metallicity4').v.tolist() ## Zr

        ## Profile plots
        print('Making profile plots for ' + run_name + ' output_' + o)
        make_prof_fields = []
        for l, prof_field in enumerate(profile_fields_postNSM):
            prof_name = 'run_data/profile_' + prof_field[1] + '_' + run_name + '_' + o + '.npy'

            if exists(prof_name) == False:
                make_prof_fields.append(prof_field)
            else:
                print('Profile plot ' + prof_field[1] + ' already made for ' + run_name + ' ' + o)
        
        if len(make_prof_fields) != 0:
            plot = yt.ProfilePlot(sp, "radius", make_prof_fields, weight_field='cell_mass')
            plot.set_unit("radius", "kpc")
            profile = plot.profiles[0]

            with open('run_data/profile_x_' + run_name + '_' + o + '.npy', 'wb') as f:
                np.save(f, np.array(profile.x))

            for l, prof_field in enumerate(make_prof_fields):
                plot.set_unit(prof_field, units[prof_field])
                with open('run_data/profile_' + prof_field[1] + '_' + run_name + '_' + o + '.npy', 'wb') as f:
                    np.save(f, profile[prof_field])
            del plot
            del profile
        else:
            print('Profile plots complete for ' + run_name + ' ' + o)

        ## Phase plots
        
        print('Making phase plots for ' + run_name + ' output_' + o)
        for k, phase_field in enumerate(phase_plot_fields_postNSM):
            for m, phase_color in enumerate(phase_plot_colors_postNSM[k]):
                x = phase_field[0]
                y = phase_field[1]
                z = phase_color

                H_name = 'run_data/phase_H_' + x[1] + '_' + y[1] + '_' + z[1] + '_'+ run_name + '_' + o + '.npy'
                yedges_name = 'run_data/phase_yedges_' + x[1] + '_' + y[1] + '_' + z[1] + '_'+ run_name + '_' + o + '.npy'
                xedges_name = 'run_data/phase_xedges_' + x[1] + '_' + y[1] + '_' + z[1] + '_'+ run_name + '_' + o + '.npy'

                if ((exists(H_name) == False) and (exists(yedges_name) == False) and (exists(xedges_name) == False)):

                    x_data = sp[x].to(units[x])
                    y_data = sp[y].to(units[y])
                    z_data = sp[z].to(units[z])

                    ## Make sure that there are no zeros in the data. Set any zeros to a tiny number.
                    x_zeroes = np.argwhere(np.array(x_data) == 0.0)
                    if len(x_zeroes) != 0:
                        for i in x_zeroes:
                            x_data[i] = ds.quan(tiny_number, units[x])

                    y_zeroes = np.argwhere(np.array(y_data) == 0.0)
                    if len(y_zeroes) != 0:
                        for i in y_zeroes:
                            y_data[i] = ds.quan(tiny_number, units[y])
                            
                    z_zeroes = np.argwhere(np.array(z_data) == 0.0)
                    if len(z_zeroes) != 0:
                        for i in z_zeroes:
                            z_data[i] = ds.quan(tiny_number, units[z])

                    x_bins = np.logspace(np.log10(np.min(x_data)), np.log10(np.max(x_data)), 128)
                    y_bins = np.logspace(np.log10(np.min(y_data)), np.log10(np.max(y_data)), 128)

                    H, yedges, xedges = np.histogram2d(x_data, y_data, [x_bins, y_bins], weights = z_data)

                    with open('run_data/phase_H_' + x[1] + '_' + y[1] + '_' + z[1] + '_'+ run_name + '_' + o + '.npy', 'wb') as f:
                        np.save(f, H)
                    with open('run_data/phase_yedges_' + x[1] + '_' + y[1] + '_' + z[1] + '_'+ run_name + '_' + o + '.npy', 'wb') as f:
                        np.save(f, yedges)
                    with open('run_data/phase_xedges_' + x[1] + '_' + y[1] + '_' + z[1] + '_'+ run_name + '_' + o + '.npy', 'wb') as f:
                        np.save(f, xedges)

                    del x
                    del y
                    del z
                    del x_data
                    del y_data
                    del z_data
                    del x_bins
                    del y_bins
                    del H
                    del yedges
                    del xedges
                else:
                    print(x[1] + '_' + y[1] + '_' + z[1] + '_'+ run_name + '_' + o + ' created.')

        
        ## Projection plots
        
        print('Making projection plots for ' + run_name + ' output_' + o)
        for field in fields_postNSM:
            img_name = 'run_data/projection_' + field[1] + '_' + run_name + '_' + o +'.npy'
            if exists(img_name) == False:
                prj = yt.ProjectionPlot(ds, 'x', field, center = halo_pos, weight_field='density', width=width, data_source = reg)
                frb = prj.data_source.to_frb(width, res)
                with open(img_name, 'wb') as f:
                    np.save(f, np.array(frb[field]))
            
                del prj
                del frb
            else:
                print(img_name + ' created.')

        ## Save the data

        sto.result_id = run_name + '_' + o
        sto.result = temp_dict
        
    else:
        ## No NSM field present - save the data

        sto.result_id = run_name + '_' + o
        sto.result = temp_dict

    ds.index.clear_all_data()


if yt.is_root():

    ## Let's reorganize the dictionaries slightly.

    ## This dictionary will contain only the halo quantities.
    halo_quantities = {}
    for i, j in enumerate(run_dirs):
        run_name = j.split('/')[-1]
        halo_quantities[run_name] = {}

    for i, run_name_ds in enumerate(quantities):
        run_name = run_name_ds.split('_')[0] + '_' + run_name_ds.split('_')[1]
        ds = run_name_ds.split('_')[-1]
        
        halo_quantities[run_name][ds] = {}
        for j, o in enumerate(quantities[run_name_ds]):
            for k, quant in enumerate(quantities[run_name_ds][o]['halo_quantities']):
                halo_quantities[run_name][ds][quant] = {}
                for l, weight in enumerate(quantities[run_name_ds][o]['halo_quantities'][quant]):
                    halo_quantities[run_name][ds][quant][weight] = quantities[run_name_ds][o]['halo_quantities'][quant][weight]

    with open('halo_quantities.json', 'w') as outfile:
        json.dump(halo_quantities, outfile)

    ## This dictionary will contain only information related to stars.
    star_info = {}
    for i, j in enumerate(run_dirs):
        run_name = j.split('/')[-1]
        star_info[run_name] = {}

    for i, run_name_ds in enumerate(quantities):
        run_name = run_name_ds.split('_')[0] + '_' + run_name_ds.split('_')[1]
        for j, ds in enumerate(quantities[run_name_ds]):
            star_info[run_name][ds] = {}
            star_info[run_name][ds]['star_info'] = {}
            for key in quantities[run_name_ds][ds]['star_info'].keys():
                star_info[run_name][ds]['star_info'][key] = quantities[run_name_ds][ds]['star_info'][key]

    with open('star_info.json', 'w') as outfile:
        json.dump(star_info, outfile)
