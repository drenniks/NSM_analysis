import os
import json
import numpy as np

fields = ['density', 'temperature', 'El_fraction', 'metallicity']
width = (10, 'kpc')

with open('DD_data.json') as f:
    DD_data = json.load(f)

with open('potential_p3_info.json') as f:
    potential_p3 = json.load(f)

p3_ids = [] 

for i, j in enumerate(potential_p3):
    try:
        star_identity = potential_p3[j]['index']
        if star_identity[0] not in p3_ids:
            p3_ids.append(star_identity[0])
    
    except:
        print('No stars in dataset ' + j)
p3_ids = np.array(p3_ids)

movie = './movie.sh'
image_output = 'images/potential_choices/'
movie_output = image_output + 'movies/'

outputs = np.array(list(DD_data.keys()))

for i in p3_ids:
    for field in fields:
        image_list = []
        for out in outputs:
            img_name = 'output_' + out + '_' + field + '_' + str(width[0]) + '_' + str(float(i)) +'.png'
            for img in os.listdir(image_output):
                if img == img_name:
                    image_list.append(image_output + img)
        if len(image_list) == 1:
            print('Only one dataset to view for star ID # ' + str(float(i)))
            continue
        else:
            movie_images = " ".join(image_list)
            movie_name = str(int(float(i))) + '_' + field + '_' + str(width[0]) + '.mov'

            execute = movie + ' ' + movie_output + movie_name + ' 5 ' + movie_images
            os.system(execute)
