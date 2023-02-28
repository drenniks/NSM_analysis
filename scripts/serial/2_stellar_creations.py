import time
startTime = time.time()

import yt
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import json

yt.enable_plugins()

### Produce histograms of creation times for P3 and P2 stars. Also print information about number of stars

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 22}

matplotlib.rc('font', **font)

data_dir = "../"

with open('DD_data_serial.json') as f:
    DD_data = json.load(f)

outputs = np.array(list(DD_data.keys()))
final_output = outputs[-1]
final_string = 'DD' + final_output + '/output_' + final_output

ds = yt.load(data_dir + final_string)

stars = ['p3', 'p2', 'p3_living']
for s in stars:
    ds.add_particle_filter(s)

ad = ds.all_data()

p2_id = ad['p2', 'particle_index']
p3_id = ad['p3', 'particle_index']
p3_living_id = ad['p3_living', 'particle_index']

p2_ct = ad['p2', 'creation_time']
p3_ct = ad['p3', 'creation_time']

p2_ct_sortby = np.argsort(p2_ct)
p3_ct_sortby = np.argsort(p3_ct)

f = open("number_of_stars.txt", "w")
f.write('Number of P3_living stars = ' + str(len(p3_living_id)) + '\n')
f.write('Number of P3 stars = ' + str(len(p3_id)) + '\n')
f.write('Number of P2 stars = ' + str(len(p2_id)) + '\n')
f.close()

## Plot P2 Creation Times histogram
fig = plt.figure(figsize=(10, 8))
plt.hist(p2_ct[p2_ct_sortby].to('Myr').v)

plt.title('P2 Creation Times')
plt.ylabel('Number of Stars')
plt.xlabel('Creation Time [Myr]')

#plt.show()
plt.savefig('images/P2_creation_times.png')


## Plot P3 Creation Times histogram
fig = plt.figure(figsize=(10, 8))
plt.hist(p3_ct[p3_ct_sortby].to('Myr').v)

plt.title('P3 Creation Times')
plt.ylabel('Number of Stars')
plt.xlabel('Creation Time [Myr]')

#plt.show()
plt.savefig('images/P3_creation_times.png')

import time
executionTime = (time.time() - startTime)
print('Execution time in seconds: ' + str(executionTime))
