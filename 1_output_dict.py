import yt
import numpy as np
import json

yt.enable_plugins()

### Find outputs, time and z -- creates output.txt

data_dir = "/storage/home/hhive1/dskinner6/data/SG64-2020/"

fns = yt.load(data_dir + "DD????/output_????")

time = []
redshift = []
output = []

for ds in fns:
    time.append(ds.current_time.to('Myr'))
    redshift.append(ds.current_redshift)
    output.append(str(ds).split('output_')[1])

time = np.array(time)
redshift = np.array(redshift)
output = np.array(output)

DD_data = {}

for o in range(len(output)):
    DD_data[output[o]] = {}
    DD_data[output[o]]['time'] = time[o]
    DD_data[output[o]]['redshift'] = redshift[o]

# Saving output dict

with open('DD_data.json', 'w') as outfile:
    json.dump(DD_data, outfile)
