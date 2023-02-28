import json

run_dir = '../'
data_dir = 'run_data/'
image_outputs = 'images/'

runs = ['run_original', 'run_fiducial', 'run_A', 'run_B', 'run_D', 'run_E']
stars = ['p3_living', 'p3_binary', 'ns_binary', 'bh', 'p2', 'dm']
types = [5, 12, 13, 6, 7, 1]

res = 1000
width = (10, 'kpc')
direction = 'x'

tiny_number = 1e-30

restart_ds = '0037'

## For 7_make_plots.py
set_limits = False
plot_stars = True

## Dictionaries
try:
    with open('param_info.json') as f: ## produced in 6a_output_dict.py
        param_info = json.load(f)
except:
    print('param_info.json doesnt exist yet -- run all_output_info.py')

try:
    with open('massive_prog_info_incl_vir.json') as f: ## produced in 3a_massive_prog.py
        massive_prog = json.load(f)
except:
        print('massive_prog_info_incl_vir.json doesnt exist yet -- 6_data_grab.py')

try:
    with open('run_halos_incl_vir.json') as f: ## produced in 6aa_halo_finder.py
        run_halos = json.load(f)
except:
    print('run_halos_incl_vir.json doesnt exist yet -- run 6aa_halo_finder.py')

try:
    with open('unit_conversion.json') as f: ## produced in unit_conversion.ipynb
        unit_convert = json.load(f)
except: 
    print('unit_conversion.json doesnt exist yet -- run unit_conversion.py')

try: 
    with open('run_DD_data.json') as f:
        DD_data = json.load(f)
except:
    print('run_DD_data.json doesnt exist yet -- run 6_data_grab.py')

try: 
    with open('halo_quantities.json') as f:
        halo_quantities = json.load(f)
except:
    print('halo_quantities.json doesnt exist yet -- run 6_data_grab.py')

try:
    with open('star_info.json') as f:
        star_info = json.load(f)
except Exception as e:
    print('star_info.json exist yet -- run 6_data_grab.py')

## Radiation fields

halo_rad_fields = [('gas', 'J21_LW'),
                   ('gas', 'J_Lyman'),
                   ('enzo', 'HI_kph'),
                   ('enzo', 'HeI_kph'),
                   ('enzo', 'HeII_kph')]

## Halo quantities to calculate and their weights, pre- and post-NSM + some radiation fields
halo_fields = [('gas', 'metallicity'), 
               ('gas', 'metallicity3'), 
               ('gas', 'gas_fraction'), 
               ('gas', 'temperature'), 
               ('gas', 'electron_fraction'), 
               ('gas', 'H2_fraction'), 
               ('gas', 'density'), 
               ('gas', 'metallicity4'), 
               ('enzo', 'NSMRProcess'),
               ('gas', 'J21_LW'),
               ('gas', 'J_Lyman'),
               ('enzo', 'HI_kph'), 
               ('enzo', 'HeI_kph'), 
               ('enzo', 'HeII_kph')]

halo_weights = [('gas', 'density'), 
                ('gas', 'cell_mass')]

total_fields = [('gas', 'metallicity4'), 
               ('enzo', 'NSMRProcess')]

## Profile plot fields
profile_fields = [('gas', 'metallicity'),
                  ('gas', 'metallicity3'),
                  ('gas', 'gas_fraction'),
                  ('gas', 'temperature'),
                  ('gas', 'electron_fraction'),
                  ('gas', 'H2_fraction'),
                  ('gas', 'density'),
                  ("gas", "metallicity4"),
                  ('enzo', 'NSMRProcess')]

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

## Phase plot fields
phase_fields = [[("gas", "density"), ("gas", "temperature"), ("gas", "cell_mass")],
                [("gas", "density"), ("gas", "temperature"), ('gas', 'metallicity3')],
                [('index', 'radius'), ('gas', 'metallicity3'), ("gas", "cell_mass")],
                [("gas", "metallicity3"), ("gas", "metallicity4"), ("gas", "cell_mass")], 
                [("gas", "density"), ("gas", "temperature"), ("gas", "metallicity4")], 
                [('index', 'radius'), ('gas', 'metallicity4'), ("gas", "cell_mass")]]

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

## Projection plot fields
proj_fields = [("gas", "density"),
          ('gas', 'temperature'),
          ('gas', 'metallicity'),
          ('gas', 'metallicity3'),
          ('gas', 'metallicity4'), 
          ('enzo', 'HI_kph'),
          ('enzo', 'HeI_kph'),
          ('enzo', 'HeII_kph'), 
          ('gas', 'gas_fraction')]

## Projection limits
zlim = {}
zlim["density"] = (1e-27, 1e-22)
zlim["temperature"] = (1e2, 1e4)
zlim["metallicity"] = (1e-8, 1e-3)
zlim["metallicity3"] = (1e-10, 1e-2)
zlim["metallicity4"] = (1e-3, 1e3)

## Star types and fields to get information for
star_types = ['p2', 'p3']
star_fields = {'particle_index':'dimensionless',
               'particle_position':'kpc',
               'particle_mass':'Msun',
               'creation_time':'Myr',
               'metallicity_fraction':'Zsun',
               'P2_metallicity_fraction':'Zsun',
               'P3_metallicity_fraction':'Zsun',
               'NSM_metallicity_fraction':'Zr'}

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
         ('index', 'radius'):"pc", 
         ('enzo', 'NSMRProcess'):"g/cm**3",
         ('enzo', 'HI_kph'):"1/s", 
         ('enzo', 'HeI_kph'):"1/s", 
         ('enzo', 'HeII_kph'):"1/s"
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
        'radius':"r (pc)",
        "HI_kph":"HI_kph",
        "HeI_kph":"HeI_kph",
        "HeII_kph":"HeII_kph" 
        }

colorbars = {"density":'viridis', 
        "temperature":'inferno', 
        "metallicity3":'magma', 
        "metallicity4":'magma',
        "metallicity":'magma', 
        "HI_kph":'viridis',
        "HeI_kph":'viridis',
        "HeII_kph":'viridis',
        "gas_fraction":'viridis'
        }

starcolors = {'p2':'black', 'p3':'red', 'chosen_stars':'white'}