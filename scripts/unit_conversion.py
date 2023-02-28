import os
import yt
import json

## Need conversion dictionary from unitary to kpc

## positions - unitary
## rvir - unitary

unit_conversion = {}

runs = ['run_original', 'run_fiducial', 'run_A', 'run_B', 'run_C', 'run_D', 'run_E']

data_dir = '../'

for i, run_name in enumerate(runs):
    unit_conversion[run_name] = {}
    for dd in os.listdir(data_dir + run_name):
        if dd.startswith('DD'):
            o = dd.split('DD')[1]
            
            unit_conversion[run_name][o] = {}
            try:
                ds = yt.load(data_dir + run_name + '/' + dd + '/output_' + o)

                length = ds.quan(1, 'unitary')

                conv = length.to('kpc').v.tolist()

                unit_conversion[run_name][o]['length'] = conv
            except Exception as e:
                print(e)

with open('unit_conversion.json', 'w') as f:
    json.dump(unit_conversion, f)