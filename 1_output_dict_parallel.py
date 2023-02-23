import time
startTime = time.time()

import yt
import numpy as np
import json

yt.enable_plugins()
yt.enable_parallelism()

### Find outputs, time and z -- creates output.txt

data_dir = "../"

fns = yt.DatasetSeries(data_dir + "DD????/output_????")#, parallel=2)

DD_data = {}

for sto, ds in fns.piter(storage=DD_data, dynamic=True):
    output = str(ds).split('output_')[1]
    time = ds.current_time.to('Myr').v.tolist()
    redshift = ds.current_redshift
    
    
    
    sto.result_id = output
    sto.result = [time, redshift]
    
if yt.is_root():
    ## Need to reorder the data
    keys = DD_data.keys()
    new_DD_data = {}

    new_keys = []
    for i in keys:
        new_keys.append(i)

    new_keys = np.sort(new_keys)

    for i, j in enumerate(new_keys):
        new_DD_data[j] = DD_data[j]

    with open('DD_data.json', 'w') as outfile:
        json.dump(new_DD_data, outfile)
    import time
    executionTime = (time.time() - startTime)
    print('Execution time in seconds: ' + str(executionTime))
