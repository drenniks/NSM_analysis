import os
import json
from os.path import exists
import load_ins as load

set_limits = load.set_limits

DD_data = load.DD_data

image_outputs = load.image_outputs

profile_fields = load.profile_fields
phase_plot_fields = load.phase_fields
fields = load.proj_fields

movie = './movie.sh'

for i, run_name in enumerate(DD_data):
    if set_limits == True:
        image_dir = run_name + '/limits/'
    else:
        image_dir = run_name + '/no_limits/'
        
    print(run_name)

    ## Projections
    print('Making projections.')
    for field in fields:
        print(field)
        movie_name = f'{image_outputs}{image_dir}movies/projection_{field[1]}.mov'

        if exists(movie_name) == False:
            movie_images = f'{image_outputs}{image_dir}projection_{field[1]}_*.png'
            execute = movie + ' ' + movie_name + ' 5 ' + movie_images
            os.system(execute)
        else:
            print(f'{movie_name} already made.')
        
    ## Phase plots
    print('Making phase.')
    for k, field_list in enumerate(phase_plot_fields):
        print(field_list)
        x_field = field_list[0]
        y_field = field_list[1]
        z_field = field_list[2]
        movie_name = f'{image_outputs}{image_dir}movies/phase_{x_field[1]}_{y_field[1]}_{z_field[1]}.mov'

        if exists(movie_name) == False:
            movie_images = f'{image_outputs}{image_dir}phase_{x_field[1]}_{y_field[1]}_{z_field[1]}_*.png'
            execute = movie + ' ' + movie_name + ' 5 ' + movie_images
            os.system(execute)
        else:
            print(f'{movie_name} already made.')
    
    ## Profile plots
    print('Making profile.')
    for l, prof_field in enumerate(profile_fields):
        print(prof_field)
        movie_name = f'{image_outputs}{image_dir}movies/profile_{prof_field[1]}.mov'   
        
        if exists(movie_name) == False:
            movie_images = f'{image_outputs}{image_dir}profile_{prof_field[1]}_{run_name}_*.png'
            execute = movie + ' ' + movie_name + ' 5 ' + movie_images
            os.system(execute)
        else:
            print(movie_name + ' already made.')
