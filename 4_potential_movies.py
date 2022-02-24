import os
import json


direction = ['x', 'y', 'z']
fields = ['density', 'temperature', 'El_fraction', 'metallicity']
width = (5, 'kpc')

with open('DD_data.json') as f:
    DD_data = json.load(f)

with open('potential_P3_info.json') as f:
    potential_p3 = json.load(f)

p3_ids = np.array(list(potential_p3.keys()))

movie = './movie.sh'
image_output = 'images/potential_choices/'
movie_output = image_output + 'movies/'

outputs = np.array(list(DD_data.keys()))

for i in p3_ids:
    for field in fields:
        for d in direction:
            image_list = []
            for out in outputs:
                img_name = str(int(float(i))) + '_output_' + out + '_' + field + '_' + str(d) + '.png'
                for img in os.listdir(image_output):
                    if img == img_name:
                        image_list.append(image_output + img)
                    
            movie_images = " ".join(image_list)
            movie_name = str(int(float(i))) + '_' + field + '_' + str(d) + '.mov'
        
            execute = movie + ' ' + movie_output + movie_name + ' 5 ' + movie_images
            os.system(execute)
