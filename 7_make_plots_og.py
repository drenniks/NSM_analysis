import yt
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
import json
import numpy as np
import numpy.ma as ma
import os
import os.path
from os.path import exists

yt.enable_plugins()
yt.enable_parallelism()

set_limits = True

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 18}

matplotlib.rc('font', **font)

run_dir = '../'
data_dir = 'run_data/'
image_outputs = 'images/'

stars = ['p3_living', 'p3_binary', 'ns_binary', 'bh', 'p2', 'dm']
types = [5, 12, 13, 6, 7, 1]

with open('run_DD_data.json') as f:
    DD_data = json.load(f)
    
with open('temp_run_info.json') as f:
    temp_run_info = json.load(f)
    
with open('halo_quantities.json') as f:
    halo_quantities = json.load(f)
    
with open('star_info.json') as f:
    star_info = json.load(f)

with open('massive_prog_info.json') as f:
    massive_prog = json.load(f)

with open('run_halos.json') as f:
    run_halos = json.load(f)

profile_fields = [("gas", "metallicity4"), ('gas', 'metallicity'), ('gas', 'metallicity3'), ('gas', 'gas_fraction'), ('gas', 'temperature'), ('gas', 'electron_fraction'), ('gas', 'H2_fraction'), ('gas', 'density')]

phase_plot_fields = [[("gas", "density"), ("gas", "temperature")], [("gas", "metallicity3"), ("gas", "metallicity4")], [('index', 'radius'), ('gas', 'metallicity3')], [('index', 'radius'), ('gas', 'metallicity4')]]
phase_plot_colors = [[("gas", "cell_mass"), ("gas", "metallicity3"), ("gas", "metallicity4")], [("gas", "cell_mass")], [("gas", "cell_mass")], [("gas", "cell_mass")]]
fields = [("gas", "density"), ('gas', 'temperature'), ('gas', 'metallicity'), ('gas', 'metallicity3'), ('gas', 'metallicity4')]


halo_quant_fields = [('gas', 'metallicity'), ('gas', 'metallicity3'), ('gas', 'gas_fraction'), ('gas', 'temperature'), ('gas', 'electron_fraction'), ('gas', 'H2_fraction'), ('gas', 'density')]
halo_quant_radfields = [('gas', 'J21_LW'), ('gas', 'J_Lyman')]
halo_quant_weights = [('gas', 'density'), ('gas', 'cell_mass')]
halo_quant_radweights = [('index', 'ones')]
halo_quant_fields_postNSM = [('gas', 'metallicity4')]


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

field_labels = {"density":r'Density $\left( \frac{g}{cm^{3}} \right)$', 
         "temperature":'Temperature (K)', 
         "cell_mass":'Mass [g]', 
         "metallicity3":r'Pop III Metallicity (Z$_{\odot}$)', 
         "metallicity4":r"NSM R-Process Metallicity (Z$_{r}$)",
         "metallicity":r"Pop II Metallicity (Z$_{\odot}$)",
        'gas_fraction':"Gas Fraction",
         'electron_fraction':"Electron Fraction",
         'H2_fraction':"H2 Fraction",
         'J21_LW':r"J$_{21}$", 
         'J_Lyman':r"J(Lyman) $\left( \frac{erg}{cm^{2}} \right)$",
        'radius':"r (pc)"
        }

colorbars = {"density":'viridis', 
         "temperature":'inferno', 
         "metallicity3":'magma', 
         "metallicity4":'magma',
         "metallicity":'magma'
        }

## Projection limits
zlim = {}
zlim["density"] = (1e-27, 1e-22)
zlim["temperature"] = (1e2, 1e4)
zlim["metallicity"] = (1e-8, 1e-3)
zlim["metallicity3"] = (1e-10, 1e-2)
zlim["metallicity4"] = (1e-3, 1e3)

## Phase plot limits
zlim_phase = {}
zlim_phase["metallicity3"] = (1e-4, 1e-1)
zlim_phase["metallicity4"] = (1e-6, 10)
zlim_phase["cell_mass"] = (1e32, 1e37)

xlim_phase = {}
xlim_phase["density"] = (1e-28, 1e-20)
xlim_phase["metallicity3"] = (1e-10, 1e-1)
xlim_phase['radius'] = (1, 1e3)

ylim_phase = {}
ylim_phase["temperature"] = (10, 1e6)
ylim_phase["metallicity4"] = (1e-10, 1e4)
ylim_phase["metallicity3"] = (1e-10, 1e-1)

## Profile plot limits
xlim_profile = (1e-3, 1)

ylim_profile = {}
ylim_profile["metallicity4"] = (1e-2, 1e3)
ylim_profile['metallicity'] = (1e-6, 1e-2)
ylim_profile['metallicity3'] = (1e-4, 5e1)
ylim_profile['gas_fraction'] = (1e-3, 1)
ylim_profile['temperature'] = (1e2, 1e5)
ylim_profile['electron_fraction'] = (1e-2, 1)
ylim_profile['H2_fraction'] = (1e-7, 1e-2)
ylim_profile['density'] = (1e-26, 1e-20)


run_names = list(run_halos.keys())

dataset_list = []

for i, run_name in enumerate(run_halos):
    if run_name == 'run_D':
        for j, ds in enumerate(run_halos[run_name]):
#        if temp_run_info[run_name][ds]['halo_info']['halo_id'] != -1:
            dataset_list.append(run_dir + run_name + '/DD' + ds + '/output_' + ds)

#for j, ds in enumerate(massive_prog):
#    run_name = 'run_original'
#    dataset_list.append(run_dir + run_name + '/DD' + ds + '/output_' + ds)

fns = yt.DatasetSeries(dataset_list)
for ds in fns.piter(dynamic=True):
    run_name = str(ds.directory).split('/')[-2]
    o = str(ds).split('output_')[-1]

    time = DD_data[run_name][o]['time']
    redshift = DD_data[run_name][o]['redshift']

    print(o)
    print(run_name)
    
    ## Halo info
    if run_name == 'run_original':
        halo_rvir = ds.quan(massive_prog[o]['rvir'], 'unitary').to('kpc')
        halo_pos = ds.arr(massive_prog[o]['position'], 'unitary').to('kpc')

    else:
        halo_rvir = ds.quan(run_halos[run_name][o]['rvir'], 'unitary').to('kpc')
        halo_pos = ds.arr(run_halos[run_name][o]['position'], 'unitary').to('kpc')
    
    sp = ds.sphere(halo_pos, halo_rvir)

    ## Star info
    ##Chosen star - if available
    plot_chosen_star = False
    chosen_id = star_info[run_name][o]['star_info']['chosen_star']['id']
    if chosen_id != -1:
        plot_chosen_star = True
        chosen_pos = ds.arr(star_info[run_name][o]['star_info']['chosen_star']['position'], 'kpc') - halo_pos
        chosen_dr = np.sqrt(((chosen_pos)**2).sum())
        chosen_type = star_info[run_name][o]['star_info']['chosen_star']['type']
        type_name = stars[np.argwhere(np.array(types) == chosen_type)[0][0]]

        chosen_pos_x = [chosen_pos[1].v]
        chosen_pos_y = [chosen_pos[2].v]

    ##P2 stars
    plot_p2_stars = False

    try:
        num_p2 = len(star_info[run_name][o]['star_info']['P2']['position'])
        plot_p2_stars = True
    except:
        plot_p2_stars == False
    
    if plot_p2_stars == True:
        p2_pos = [ds.arr(i, 'kpc') for i in star_info[run_name][o]['star_info']['P2']['position']]
        p2_plot_pos = (p2_pos - halo_pos)
        p2_dr = np.sqrt(((p2_plot_pos)**2).sum(1))
        new_p2_plot_pos = []
        for dr in p2_dr:
            if dr < 5:
                index = np.argwhere(p2_dr == dr)[0][0]
                new_p2_plot_pos.append(p2_plot_pos[index])

        p2_list_x = [i[1].v for i in new_p2_plot_pos]
        p2_list_y = [i[2].v for i in new_p2_plot_pos]

    ##P3 stars
    plot_p3_stars = False
    p3_pos = [ds.arr(i, 'kpc') for i in star_info[run_name][o]['star_info']['P3']['position']]
    if len(p3_pos) != 0:
        plot_p3_stars = True
        p3_plot_pos = (p3_pos - halo_pos)
        p3_dr = np.sqrt(((p3_plot_pos)**2).sum(1))
        new_p3_plot_pos = []
        for dr in p3_dr:
            if dr < 5:
                index = np.argwhere(p3_dr == dr)[0][0]
                new_p3_plot_pos.append(p3_plot_pos[index])

        p3_list_x = [i[1].v for i in new_p3_plot_pos]
        p3_list_y = [i[2].v for i in new_p3_plot_pos]
    
    ## Projections
    print('Making projections')
    for field in fields:
        proj_name = data_dir + 'projection_' + field[1] + '_' + run_name + '_' + o + '.npy'
        if set_limits == True:
            proj_img_name = image_outputs + run_name + '/limits/projection_' + field[1] + '_' + o + '.png'
        else:
            proj_img_name = image_outputs + run_name + '/no_limits/projection_' + field[1] + '_' + o + '.png'

        if exists(proj_name) == True: ## See if the data exists
            if exists(proj_img_name) == False: ## See if the image hasn't been made

                cbar_label =  field_labels[field[1]]
                vmin = zlim[field[1]][0]
                vmax = zlim[field[1]][1]
                cmap = colorbars[field[1]]
                try:
                    with open(proj_name, 'rb') as f:
                        frb = np.load(f)
                except Exception as e:
                    print(run_name + ' DD' + o + ': ' + str(e))
                    continue

                plot_time_label = 't = ' + str(round(time, 2)) + ' Myr \nz = ' + str(round(redshift, 2))
                bbox = dict(boxstyle='square', fc='white', ec='white', alpha=0.5)
                radius = halo_rvir * 100
                

                fig = plt.figure(figsize=(10,10))
                ax = plt.subplot()

                if set_limits == True:
                    img = ax.imshow(frb, origin="lower",extent=[-5, 5, -5, 5], cmap = cmap, norm=LogNorm(vmin=vmin, vmax=vmax))
                else:
                    img = ax.imshow(frb, origin="lower",extent=[-5, 5, -5, 5], cmap = cmap, norm=LogNorm())
                if plot_p2_stars == True:
                    plt.scatter(p2_list_x, p2_list_y, marker = '*', color = 'black', alpha=0.5)
                if plot_p3_stars == True:
                    plt.scatter(p3_list_x, p3_list_y, marker = '*', color = 'red', alpha=0.5)
                if plot_chosen_star == True:
                    plt.scatter(chosen_pos_x, chosen_pos_y, marker = "*", color = 'white')
                    ax.set_title(type_name)
                    
                plt.plot(0, 0, 'o', ms=radius, mec='white', mfc='none')   
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="5%", pad=0)
                cbar = plt.colorbar(img, cax = cax)
                cbar.set_label(cbar_label)

                ax.set_xlabel('y (kpc)')
                ax.set_ylabel('z (kpc)')
                ax.text(-4.8, 4.15, plot_time_label, color='black', bbox=bbox)
                #plt.show()
                fig.tight_layout()

                plt.savefig(proj_img_name, dpi=300)
            else:
                print(proj_img_name + ' has been created.')

    if plot_p2_stars == True:
        del p2_list_x
        del p2_list_y
    
    if plot_p3_stars == True:
        del p3_list_x
        del p3_list_y
    
    ## Phase Plots
    print('Making phase plots')
    for k, field_list in enumerate(phase_plot_fields):
        x_field = field_list[0]
        y_field = field_list[1]

        for z_field in phase_plot_colors[k]:
            phase_name = data_dir + 'phase_xedges_'+ x_field[1] + '_' + y_field[1] + '_'+ z_field[1] + '_' + run_name + '_' + o + '.npy'
            if set_limits == True:
                phase_img_name = image_outputs + run_name + '/limits/phase_'+ x_field[1] + '_' + y_field[1] + '_'+ z_field[1] + '_' + run_name + '_' + o + '.png'
            else:
                phase_img_name = image_outputs + run_name + '/no_limits/phase_'+ x_field[1] + '_' + y_field[1] + '_'+ z_field[1] + '_' + run_name + '_' + o + '.png'
            
            if exists(phase_name) == True: ## See if the data exists
                if exists(phase_img_name) == False: ## See if the image exists
                    try:
                        with open(phase_name, 'rb') as f:
                            xedges = np.load(f)

                        with open(data_dir + 'phase_yedges_'+ x_field[1] + '_' + y_field[1] + '_'+ z_field[1] + '_' + run_name + '_' + o + '.npy', 'rb') as f:
                            yedges = np.load(f)

                        with open(data_dir + 'phase_H_'+ x_field[1] + '_' + y_field[1] + '_'+ z_field[1] + '_' + run_name + '_' + o + '.npy', 'rb') as f:
                            H = np.load(f)

                        vmin = zlim_phase[z_field[1]][0]
                        vmax = zlim_phase[z_field[1]][1]

                        X, Y = np.meshgrid(xedges, yedges)

                        plot_time_label = 't = ' + str(round(time, 2)) + ' Myr \nz = ' + str(round(redshift, 2))
                        bbox = dict(boxstyle='square', fc='black', ec='black', alpha=0.5)

                        fig = plt.figure(figsize=(10,10))
                        ax = plt.subplot()

                        if set_limits == True:
                            img = ax.pcolormesh(Y, X, H, norm=LogNorm(vmin = vmin, vmax = vmax))
                            ax.set_xlim(xlim_phase[x_field[1]][0], xlim_phase[x_field[1]][1])
                            ax.set_ylim(ylim_phase[y_field[1]][0], ylim_phase[y_field[1]][1])
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
                        if plot_chosen_star == True:
                            ax.set_title(type_name)
                        cbar.set_label(field_labels[z_field[1]])

                        fig.tight_layout()

                        plt.savefig(phase_img_name, dpi=300)
                        
                    except:
                        try:
                            plot = yt.PhasePlot(sp, x_field, y_field, z_field, weight_field = None)
                            if set_limits == True:
                                plot.save(image_outputs + run_name + '/limits/phase_plot_'+ x_field[1] + '_' + y_field[1] + '_'+ z_field[1] + '_' + run_name + '_' + o + '.png')   
                            else:
                                plot.save(image_outputs + run_name + '/no_limits/phase_plot_'+ x_field[1] + '_' + y_field[1] + '_'+ z_field[1] + '_' + run_name + '_' + o + '.png')   
                        except Exception as e:
                            print(run_name + ' DD' + o + ': ' + str(e))
    ## Profile Plots
    print('Making profile plots')
    try:
        profile_x = data_dir + 'profile_x_' + run_name + '_' + o + '.npy'
        with open(profile_x, 'rb') as f:
            xx = np.load(f)
    
        plot_time_label = 't = ' + str(round(time, 2)) + ' Myr \nz = ' + str(round(redshift, 2))
        bbox = dict(boxstyle='square', fc='grey', ec='black', alpha=0.5)
        for l, prof_field in enumerate(profile_fields):
            profile_y = data_dir + 'profile_' + prof_field[1] + '_' + run_name + '_' + o + '.npy'

            if set_limits == True: 
                profile_img_name = image_outputs + run_name + '/limits/profile_' + prof_field[1] + '_' + run_name + '_' + o + '.png'
            else:
                profile_img_name = image_outputs + run_name + '/no_limits/profile_' + prof_field[1] + '_' + run_name + '_' + o + '.png'
        
            if exists(profile_y) == True: ## See if the data is there
                if exists(profile_img_name) == False: ## See if the image is there

                    with open(profile_y, 'rb') as f:
                        y = np.load(f)
        
                    fig = plt.figure(figsize=(10, 8))
                    ax = plt.subplot()
                    
                    plt.loglog(xx, y)
                
                    if plot_chosen_star == True:
                        plt.axvline(x = chosen_dr, color='black')
                        plt.title(type_name)
                    plt.xlabel(r'Radius (kpc)')
                    plt.ylabel(field_labels[prof_field[1]])
                    plt.text(0.02, 0.90, plot_time_label, color='black', bbox=bbox, transform=ax.transAxes)

                    if set_limits == True:
                        plt.xlim(xlim_profile[0], xlim_profile[1])
                        plt.ylim(ylim_profile[prof_field[1]][0], ylim_profile[prof_field[1]][1])

                    fig.tight_layout()
                    #plt.show()
                    plt.savefig(profile_img_name, dpi=300)
    except Exception as e:
        print(run_name + ' DD' + o + ': ' + str(e))
