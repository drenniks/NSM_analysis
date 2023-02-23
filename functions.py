import yt
import numpy as np
import json
import sys
import os
from os.path import exists
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import load_ins as load

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 18}

matplotlib.rc('font', **font)

run_dir = load.run_dir
data_dir = load.data_dir
image_outputs = load.image_outputs

runs = load.runs
stars = load.stars
types = load.types

res = load.res
width = load.width

tiny_number = load.tiny_number

restart_ds = load.restart_ds

DD_data = load.DD_data
param_info = load.param_info
massive_prog = load.massive_prog
run_halos = load.run_halos
unit_convert = load.unit_convert

## These are for 7_make_plots -- created in 6_data_grab
try: 
    halo_quantities = load.halo_quantities
    star_info = load.star_info
except Exception as e:
    print('Star_info and halo_quantities dont exist yet.')

halo_rad_fields = load.halo_rad_fields

##################### GRABBING DATA LIST #####################

def data_list():
    """
    Returns a yt.DatasetSeries for the runs. 
    This is the same structure used for finding the halo positions.
    """
    dataset_list = []

    for i, run_name in enumerate(runs):
        if run_name == 'run_original':
            for j, o in enumerate(massive_prog):
                if int(o) > int(restart_ds):
                    pass
                else:
                    dataset_list.append(f'{run_dir}{run_name}/DD{o}/output_{o}')
            for j, o in enumerate(run_halos[run_name]):
                if int(o) >= int(restart_ds):
                    dataset_list.append(f'{run_dir}{run_name}/DD{o}/output_{o}')
                else:
                    pass
        else:
            for j, o in enumerate(run_halos[run_name]):
                dataset_list.append(f'{run_dir}{run_name}/DD{o}/output_{o}')

    fns = yt.DatasetSeries(dataset_list)
    
    return fns

##################### GRABBING DATA LIST #####################

def data_list_noyt():
    """
    Returns a yt.DatasetSeries for the runs. 
    This is the same structure used for finding the halo positions.
    """
    dataset_list = []

    for i, run_name in enumerate(runs):
        if run_name == 'run_original':
            for j, o in enumerate(massive_prog):
                if int(o) > int(restart_ds):
                    pass
                else:
                    dataset_list.append(f'{run_dir}{run_name}/DD{o}/output_{o}')
            for j, o in enumerate(run_halos[run_name]):
                if int(o) >= int(restart_ds):
                    dataset_list.append(f'{run_dir}{run_name}/DD{o}/output_{o}')
                else:
                    pass
        else:
            for j, o in enumerate(run_halos[run_name]):
                dataset_list.append(f'{run_dir}{run_name}/DD{o}/output_{o}')

    fns = np.array(dataset_list)
    
    return fns
    
##################### GRABBING HALO INFO #####################

def get_halo(ds, run_name, o):
    """
    Inputs: run_name, o
    Returns: halo_pos, halo_rvir, halo_mass
    
    This grabs the halo information for the run and output.
    For the original run, we grab the information from massive_prog. After that we grab from run_halos.
    """
    if run_name == 'run_original':
        if int(o) > int(restart_ds):
            halo_pos = ds.arr(run_halos[run_name][o]['position'], 'unitary').to('kpc')
            halo_rvir = ds.quan(run_halos[run_name][o]['rvir'], 'unitary').to('kpc')
            halo_mass = ds.quan(run_halos[run_name][o]['mass'], 'Msun')
        else:
            halo_pos = ds.arr(massive_prog[o]['position'], 'unitary').to('kpc')
            halo_rvir = ds.quan(massive_prog[o]['rvir'], 'unitary').to('kpc')
            halo_mass = ds.quan(massive_prog[o]['mass'], 'Msun')
    else:
        halo_pos = ds.arr(run_halos[run_name][o]['position'], 'unitary').to('kpc')
        halo_rvir = ds.quan(run_halos[run_name][o]['rvir'], 'unitary').to('kpc')
        halo_mass = ds.quan(run_halos[run_name][o]['mass'], 'Msun')
        
    return halo_pos, halo_rvir, halo_mass

##################### GRABBING HALO INFO WITHOUT LOADING IN DATASET ##################### 

def get_halo_nods(run_name, o):
    """
    Inputs: run_name, o
    Returns: halo_pos, halo_rvir, halo_mass
    
    This grabs the halo information for the run and output.
    For the original run, we grab the information from massive_prog. After that we grab from run_halos.
    """
    if run_name == 'run_original':
        if int(o) > int(restart_ds):
            halo_pos = np.array(run_halos[run_name][o]['position']) * unit_convert[run_name][o]['length']
            halo_rvir = np.array(run_halos[run_name][o]['rvir']) * unit_convert[run_name][o]['length']
            halo_mass = run_halos[run_name][o]['mass']
        else:
            halo_pos = np.array(massive_prog[o]['position']) * unit_convert[run_name][o]['length']
            halo_rvir = np.array(massive_prog[o]['rvir']) * unit_convert[run_name][o]['length']
            halo_mass = massive_prog[o]['mass']
    else:
        halo_pos = np.array(run_halos[run_name][o]['position']) * unit_convert[run_name][o]['length']
        halo_rvir = np.array(run_halos[run_name][o]['rvir']) * unit_convert[run_name][o]['length']
        halo_mass = run_halos[run_name][o]['mass']
        
    return halo_pos, halo_rvir, halo_mass

##################### GETTING HALO QUANTITY INFO #####################

def get_halo_info(ds, sp, fields, weights, units):
    '''
    Input: ds, sp (yt data object), fields, weights (list of weight fields), units (dictionary of units).
    Returns: a dictionary with structure dictionary[field][weight] = value
    '''
    halo_info = {}
    for field in fields:
        ## Check if field exists
        field_lives = field in ds.derived_field_list
        rad_field = field in halo_rad_fields
        
        halo_info[field[1]] = {}
        if field_lives == True:
            if rad_field == True:
                value = sp.quantities.weighted_average_quantity(field, weight=('index', 'ones')).to(units[field]).v.tolist()
                halo_info[field[1]]['ones'] = value
            
            else:
                for weight in weights:
                    value = sp.quantities.weighted_average_quantity(field, weight=weight).to(units[field]).v.tolist()

                    halo_info[field[1]][weight[1]] = value
        else:
            if rad_field == True:
                halo_info[field[1]]['ones'] = -1
            else:
                for weight in weights:
                    halo_info[field[1]][weight[1]] = -1
    
    return halo_info

##################### GETTING TOTAL HALO QUANTITY INFO #####################

def get_total_info(ds, sp, fields, units):
    '''
    Input: ds, sp (yt data object), fields (list of fields), units (dictionary of units).
    Returns: a dictionary with structure dictionary[field][weight] = value
    '''
    halo_info = {}
    for field in fields:
        ## Check if field exists
        field_lives = field in ds.derived_field_list

        if field_lives == True:
            value = sp.quantities.total_quantity(field).to(units[field]).v.tolist()
            halo_info[field[1]] = value
        
        else:
            halo_info[field[1]] = -1
    
    return halo_info

##################### GRABBING PROFILE PLOT DATA #####################

def profile_plot_data(ds, sp, fields, run_name, o, units):
    '''
    Input: ds, sp (yt data object), fields (list of fields), units (dictionary of units).
    Returns: Saves numpy arrays to data_dir
    '''
    for field in fields:
        field_lives = field in ds.derived_field_list
        
        if field_lives == True:
            prof_name = f'{data_dir}profile_{field[1]}_{run_name}_{o}.npy'
            
            if exists(prof_name) == False:
                try:
                    plot = yt.ProfilePlot(sp, "radius", field, weight_field='cell_mass')
                    plot.set_unit("radius", "kpc")
                    profile = plot.profiles[0]

                    with open(f'{data_dir}profile_x_{run_name}_{o}.npy', 'wb') as f:
                        np.save(f, np.array(profile.x))

                    for l, prof_field in enumerate(field):
                        plot.set_unit(prof_field, units[field])
                        with open(f'{data_dir}profile_{field[1]}_{run_name}_{o}.npy', 'wb') as f:
                            np.save(f, profile[field])
                    del plot
                    del profile
                except Exception as e:
                    print(f'{run_name} DD{o}: {str(e)}')
            else:
                print(f'Profile plot {field} already made for {run_name} {o}')
        else:
            print(f'Field {field[1]} doesnt exist for {run_name} {o}.')
    return 

##################### GRABBING PHASE PLOT DATA #####################

def phase_plot_data(ds, sp, fields, units, run_name, o):
    for k, fields in enumerate(fields):
        x = fields[0]
        y = fields[1]
        z = fields[2]

        x_lives = x in ds.derived_field_list
        y_lives = y in ds.derived_field_list
        z_lives = z in ds.derived_field_list

        if (x_lives == True) & (y_lives == True) & (z_lives == True):
            H_name = f'{data_dir}phase_H_{x[1]}_{y[1]}_{z[1]}_{run_name}_{o}.npy'
            yedges_name = f'{data_dir}phase_yedges_{x[1]}_{y[1]}_{z[1]}_{run_name}_{o}.npy'
            xedges_name = f'{data_dir}phase_xedges_{x[1]}_{y[1]}_{z[1]}_{run_name}_{o}.npy'

            if ((exists(H_name) == False) and (exists(yedges_name) == False) and (exists(xedges_name) == False)):
                try:
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

                    with open(f'run_data/phase_H_{x[1]}_{y[1]}_{z[1]}_{run_name}_{o}.npy', 'wb') as f:
                        np.save(f, H)
                    with open(f'run_data/phase_yedges_{x[1]}_{y[1]}_{z[1]}_{run_name}_{o}.npy', 'wb') as f:
                        np.save(f, yedges)
                    with open(f'run_data/phase_xedges_{x[1]}_{y[1]}_{z[1]}_{run_name}_{o}.npy', 'wb') as f:
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
                except Exception as e:
                    print(f'{run_name} {o} phase plot {x[1]} {y[1]} {z[1]} couldnt be made: {e}')
            else:
                print(f'{x[1]}_{y[1]}_{z[1]}_{run_name}_{o} created.')
        else:
            field_truth = np.array([x_lives, y_lives, z_lives])
            
            false = field_truth == False
            for i, f in enumerate(false):
                if f == True:
                    print(f'{run_name} {o} field {fields[i]} doesnt exist')
            print('\n')

##################### GRABBING PROJECTION PLOT DATA #####################

def projection_plot_data(ds, reg, center, fields, units, run_name, o, width = (10, 'kpc'), res = 1000):
    for field in fields:
        field_exists = field in ds.derived_field_list
        if field_exists == True:
            try:
                img_name = f'{data_dir}projection_{field[1]}_{run_name}_{o}.npy'
                if exists(img_name) == False:
                    rad_field = field in load.halo_rad_fields
                    if rad_field == True:
                        weight = 'ones'
                    else:
                        weight = 'density'
                    prj = yt.ProjectionPlot(ds, 'x', field, center = center, weight_field=weight, width=width, data_source = reg)
                    prj.set_unit(field, units[field])
                    frb = prj.data_source.to_frb(width, res)
                    with open(img_name, 'wb') as f:
                        np.save(f, np.array(frb[field]))
                    del prj
                    del frb
                else:
                    print(f'{img_name} created.')
            except:
                print(f'ERROR {run_name} {o} projection plot {field[1]}.')
        else:
            print(f'{run_name} {o} field {field} doesnt exist.')

##################### GRABBING STAR DATA #####################

def star_data(ds, ad, fields, star_types, particular_stars = -1):
    star_dict = {}
    for stype in star_types:
        star_dict[stype] = {}
        ## Check to see if there are particles of this type
        num_stars = len(ad[stype, 'particle_index'])
        if num_stars > 0:
            for i, field in enumerate(fields):
                field_exists = (stype, field) in ds.derived_field_list
                
                if field_exists == True:
                    value = ad[stype, field].to(fields[field]).v.tolist()
                    star_dict[stype][field] = value
                else:
                    star_dict[stype][field] = list(np.zeros(num_stars) - 1)
        else:
            star_dict[stype] = -1
    
    if particular_stars != -1:
        star_dict['chosen_stars'] = {}
        for star in particular_stars:
            star_dict['chosen_stars'][str(star)] = {}
            
            part_index = np.argwhere(ad['all', 'particle_index'] == int(star))
            
            if len(part_index) != 0:
                part_type = ad['all', 'particle_type'][part_index[0][0]]
                part_pos = ad['all', 'particle_position'][part_index[0][0]]
                part_mass = ad['all', 'particle_mass'][part_index[0][0]]
                part_ct = ad['all', 'creation_time'][part_index[0][0]]
                
                star_dict['chosen_stars'][str(star)]['particle_position'] = part_pos.to('kpc').v.tolist()
                star_dict['chosen_stars'][str(star)]['particle_mass'] = part_mass.to('Msun').v.tolist()
                star_dict['chosen_stars'][str(star)]['creation_time'] = part_ct.to('Myr').v.tolist()
                star_dict['chosen_stars'][str(star)]['particle_type'] = part_type.v.tolist()
            else:
                star_dict['chosen_stars'][str(star)]['particle_position'] = -1
                star_dict['chosen_stars'][str(star)]['particle_mass'] = -1
                star_dict['chosen_stars'][str(star)]['creation_time'] = -1
                star_dict['chosen_stars'][str(star)]['particle_type'] = -1
                
    return star_dict

##################### MAKING PROJECTIONS #####################

def make_projections(run_name, o, halo_rvir, halo_pos, time, redshift, fields, field_labels, zlim, colorbars, starcolors, set_limits = False):
    """
    Input:
        run_name: name of the run (e.g. 'run_A')
        o: output (e.g. '0050')
        halo_rvir: halo virial radius (no units -- assumes default from 6_data_grab)
        halo_pos: halo position (no units -- assumes default from 6_data_grab)
        time: time from DD_data
        redshift: redshift from DD_data
        fields: field list
        field_labels: colorbar label
        zlim: colorbar limits
        colorbars: colormap of colorbar
        starcolors: colors of particle types you want to plut
        set_limits: True or False
    """
    for field in fields:
        proj_name = f'{data_dir}projection_{field[1]}_{run_name}_{o}.npy'
        
        if exists(proj_name) == True: ## Data exists
            if set_limits == True:
                proj_img_name = f'{image_outputs}{run_name}/limits/projection_{field[1]}_{o}.png'
            else:
                proj_img_name = f'{image_outputs}{run_name}/no_limits/projection_{field[1]}_{o}.png'
        
            if exists(proj_img_name) == False: ## Image doesn't exist
                try:
                    cbar_label =  field_labels[field[1]]
                    cmap = colorbars[field[1]]
                    try:
                        with open(proj_name, 'rb') as f:
                            frb = np.load(f)
                    except Exception as e:
                        print(f'ERROR {run_name} {o}: {str(e)}')
                        continue

                    plot_time_label = 't = ' + str(round(time, 2)) + ' Myr \nz = ' + str(round(redshift, 2))
                    bbox = dict(boxstyle='square', fc='white', ec='white', alpha=0.5)
                    radius = halo_rvir * 100

                    fig = plt.figure(figsize=(10,10))
                    ax = plt.subplot()

                    if set_limits == True:
                        vmin = zlim[field[1]][0]
                        vmax = zlim[field[1]][1]
                        img = ax.imshow(frb, origin="lower",extent=[-5, 5, -5, 5], cmap = cmap, norm=LogNorm(vmin=vmin, vmax=vmax))
                    else:
                        img = ax.imshow(frb, origin="lower",extent=[-5, 5, -5, 5], cmap = cmap, norm=LogNorm())

                    ## See if there are stars to plot

                    star_types = list(star_info[run_name][o].keys())

                    for star in star_types:
                        star_dict = star_info[run_name][o][star]

                        if star_dict != -1:
                            if star != 'chosen_stars':
                                position = np.array(star_info[run_name][o][star]['particle_position']) - halo_pos
                                pdr = np.sqrt(((position)**2).sum(1))

                                new_pos = []
                                for dr in pdr:
                                    if dr < 5:
                                        index = np.argwhere(pdr == dr)[0][0]
                                        new_pos.append(position[index])

                                list_x = [i[1] for i in new_pos]
                                list_y = [i[2] for i in new_pos]

                                #posx = position[1]
                                #posy = position[2]

                                c = starcolors[star]

                                plt.scatter(list_x, list_y, marker = '*', color = c, alpha=0.5)
                            else:
                                for l, part in enumerate(star_info[run_name][o][star]):
                                    position = star_info[run_name][o][star][part]['particle_position'] - halo_pos
                                    ptype = star_info[run_name][o][star][part]['particle_type']

                                    posx = position[1]
                                    posy = position[2]
                                    
                                    c = starcolors[star]

                                    plt.scatter(posx, posy, marker = '*', color = c, alpha=1.0)

                    plt.plot(0, 0, 'o', ms=radius, mec='white', mfc='none')   
                    divider = make_axes_locatable(ax)
                    cax = divider.append_axes("right", size="5%", pad=0)
                    cbar = plt.colorbar(img, cax = cax)
                    cbar.set_label(cbar_label)

                    ax.set_xlabel('y (kpc)')
                    ax.set_ylabel('z (kpc)')
                    ax.text(-4.8, 4.15, plot_time_label, color='black', bbox=bbox)

                    fig.tight_layout()
                    #plt.show()
                    plt.savefig(proj_img_name, dpi=300)

                    plt.close()
                    print(f'{proj_img_name} is done!')
                except Exception as e:
                    print(f'ERROR {run_name} {o} e: {str(e)}')
                    exc_type, exc_obj, exc_tb = sys.exc_info()
                    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                    print(exc_type, fname, exc_tb.tb_lineno)
                    print(
                    type(e).__name__,          
                    __file__,                  
                    e.__traceback__.tb_lineno  
                    )
            else:
                print(f'{proj_img_name} has been created.')

##################### MAKING PHASE PLOTS #####################

def make_phase(run_name, o, time, redshift, fields, field_labels, xlim, ylim, zlim, set_limits = False, chosen_star = -1):
    """
    Input:
        run_name: name of the run (e.g. 'run_A')
        o: output (e.g. '0050')
        time: time from DD_data
        redshift: redshift from DD_data
        fields: field list
        field_labels: colorbar label
        xlim: xlimits
        ylim: ylimits
        zlim: colorbar limits
        set_limits: True or False
        chosen_star = -1 (If set to the particle ID of the desired star, plot title will be the particle type name)
    """
    for k, field in enumerate(fields):
            x_field = field[0]
            y_field = field[1]
            z_field = field[2]

            phase_name = f'{data_dir}phase_xedges_{x_field[1]}_{y_field[1]}_{z_field[1]}_{run_name}_{o}.npy'
            
            if set_limits == True:
                phase_img_name = f'{image_outputs}{run_name}/limits/phase_{x_field[1]}_{y_field[1]}_{z_field[1]}_{run_name}_{o}.png'
            else:
                phase_img_name = f'{image_outputs}{run_name}/no_limits/phase_{x_field[1]}_{y_field[1]}_{z_field[1]}_{run_name}_{o}.png'

            if exists(phase_name) == True: ## See if the data exists
                if exists(phase_img_name) == False: ## See if the image exists
                    try:
                        with open(phase_name, 'rb') as f:
                            xedges = np.load(f)

                        with open(f'{data_dir}phase_yedges_{x_field[1]}_{y_field[1]}_{z_field[1]}_{run_name}_{o}.npy', 'rb') as f:
                            yedges = np.load(f)

                        with open(f'{data_dir}phase_H_{x_field[1]}_{y_field[1]}_{z_field[1]}_{run_name}_{o}.npy', 'rb') as f:
                            H = np.load(f)

                        vmin = zlim[z_field[1]][0]
                        vmax = zlim[z_field[1]][1]

                        X, Y = np.meshgrid(xedges, yedges)

                        plot_time_label = 't = ' + str(round(time, 2)) + ' Myr \nz = ' + str(round(redshift, 2))
                        bbox = dict(boxstyle='square', fc='black', ec='black', alpha=0.5)

                        fig = plt.figure(figsize=(10,10))
                        ax = plt.subplot()

                        if set_limits == True:
                            img = ax.pcolormesh(Y, X, H, norm=LogNorm(vmin = vmin, vmax = vmax))
                            ax.set_xlim(xlim[x_field[1]][0], xlim[x_field[1]][1])
                            ax.set_ylim(ylim[y_field[1]][0], ylim[y_field[1]][1])
                        else:
                            img = ax.pcolormesh(Y, X, H, norm=LogNorm())
                        ax.set_yscale('log')
                        ax.set_xscale('log')
                        divider = make_axes_locatable(ax)
                        cax = divider.append_axes("right", size="5%", pad=0)
                        cbar = plt.colorbar(img, cax = cax)

                        plt.text(0.02, 0.92, plot_time_label, color='white', bbox=bbox, transform=ax.transAxes)
                        ax.set_xlabel(field_labels[x_field[1]])
                        ax.set_ylabel(field_labels[y_field[1]])

                        if chosen_star != -1:
                            star_type = star_info[run_name][o]['chosen_stars'][chosen_star]['particle_type']
                            type_name = np.array(stars)[np.argwhere(np.array(types) == star_type)[0][0]]

                            ax.set_title(type_name)
                        cbar.set_label(field_labels[z_field[1]])

                        fig.tight_layout()

                        plt.savefig(phase_img_name, dpi=300)
                        plt.close()
                        print(f'{phase_img_name} is done!')
                    except Exception as e:
                        print(f'ERROR {run_name} {o} e: {str(e)}')
                        print(f'x: {x_field}, y {y_field}, z {z_field}')
                        exc_type, exc_obj, exc_tb = sys.exc_info()
                        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                        print(exc_type, fname, exc_tb.tb_lineno)
                        print(
                        type(e).__name__,          
                        __file__,                  
                        e.__traceback__.tb_lineno  
                    )
                else:
                    print(f'{phase_img_name} has been created.')

##################### MAKING PROFILE PLOTS #####################
xlim_profile = load.xlim_profile
ylim_profile = load.ylim_profile

def make_profile(run_name, o, time, redshift, fields, field_labels, xlim = xlim_profile, ylim = ylim_profile, set_limits = False, chosen_star = -1, halo_pos = -1):
    try:
        profile_x = f'{data_dir}profile_x_{run_name}_{o}.npy'
        with open(profile_x, 'rb') as f:
            xx = np.load(f)

        plot_time_label = 't = ' + str(round(time, 2)) + ' Myr \nz = ' + str(round(redshift, 2))
        bbox = dict(boxstyle='square', fc='grey', ec='black', alpha=0.5)
        for l, field in enumerate(fields):
            profile_y = f'{data_dir}profile_{field[1]}_{run_name}_{o}.npy'

            if set_limits == True: 
                profile_img_name = f'{image_outputs}{run_name}/limits/profile_prof_{field[1]}_{run_name}_{o}.png'
            else:
                profile_img_name = f'{image_outputs}{run_name}/no_limits/profile_prof_{field[1]}_{run_name}_{o}.png'

            if exists(profile_y) == True: ## See if the data is there
                if exists(profile_img_name) == False: ## See if the image is there

                    with open(profile_y, 'rb') as f:
                        y = np.load(f)

                    fig = plt.figure(figsize=(10, 8))
                    ax = plt.subplot()

                    plt.loglog(xx, y)

                    if chosen_star != -1:
                        chosen_pos = np.array(star_info[run_name][o]['chosen_stars'][chosen_star]['particle_position']) - halo_pos
                        chosen_dr = np.sqrt(((chosen_pos)**2).sum())
                        
                        star_type = star_info[run_name][o]['chosen_stars'][chosen_star]['particle_type']
                        type_name = np.array(stars)[np.argwhere(np.array(types) == star_type)[0][0]]

                        
                        plt.axvline(x = chosen_dr, color='black')
                        plt.title(type_name)
                    plt.xlabel(r'Radius (kpc)')
                    plt.ylabel(field_labels[field[1]])
                    plt.text(0.02, 0.90, plot_time_label, color='black', bbox=bbox, transform=ax.transAxes)

                    if set_limits == True:
                        plt.xlim(xlim[0], xlim[1])
                        plt.ylim(ylim[field[1]][0], ylim[field[1]][1])

                    fig.tight_layout()
                    plt.savefig(profile_img_name, dpi=300)
                    
                    plt.close()
                    print(f'{profile_img_name} is done!')
                else:
                    print(f'{profile_img_name} has been created.')
    except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(f'ERROR {run_name} {o} e: {str(e)}')
        print(exc_type, fname, exc_tb.tb_lineno)
        print(
        type(e).__name__,          
        __file__,                  
        e.__traceback__.tb_lineno  
    )
