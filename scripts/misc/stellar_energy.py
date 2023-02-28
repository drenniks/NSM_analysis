import yt
import numpy as np
import matplotlib.pyplot as plt
import json
import os
import os.path
from os.path import exists
from unyt import mass_hydrogen_cgs, eV, erg, Msun, s, Myr
from bisect import bisect
yt.enable_plugins()
yt.enable_parallelism()

run_dir = '../'
data_dir = 'run_data/'
image_outputs = 'images/'

tiny_number = 1e-30

stars = ['p3_age', 'p3_binary', 'ns_binary', 'bh', 'p2', 'dm']
types = [5, 12, 13, 6, 7, 1]

kph_fields = [('enzo', 'HI_kph'), ('enzo', 'HeI_kph'), ('enzo', 'HeII_kph')]
eng = [21.62, 30, 60] * eV
eng_P3 = [28.0, 30.0, 58.0, 12.8] * eV
StarClusterSNEnergy = 1e49 * erg / Msun ## erg / Msun -- releases this energy constantly over 16 Myr
StarClusterIonizingLuminosity = 1.9e+46 ## Ionizing photons / Msun / s
EnergyFractionHeI = 0.2951
EnergyFractionHeII = 2.818e-4

TypeIILowerMass = 11
TypeIIUpperMass = 40.1
PISNLowerMass = 140
PISNUpperMass = 260
PopIIISupernovaRadius = 10
SNExplosionMass = [19.99, 25, 30, 35, 40.01]

## In type II range but less than 20 Msun
SN_type_II = 1e51 * erg
SNExplosionEnergy = SN_type_II

## In type II range but greater than 20 Msun
def hypernova(m):
    index = bisect(SNExplosionMass, m)
    frac = SNExplosionMass[index] - m / (SNExplosionMass[index] - SNExplosionMass[index-1])
    SNEnergy = 1e51 #* (SNExplosionEnergy[index-1] + frac * (SNExplosionEnergy[index] - SNExplosionEnergy[index-1]))

    return SNEnergy*erg

## If PISN
def PISN(m):
    HeliumCoreMass = (13./24.) * (m - 20)
    SNEnergy = (5.0 + 1.304 * (HeliumCoreMass - 64)) * 1e51
    
    return SNEnergy*erg

with open('run_DD_data.json') as f:
    DD_data = json.load(f)
    
    
with open('star_info.json') as f:
    star_info = json.load(f)


stellar_energy = {}

#for sto, ds in fns.piter(storage = stellar_energy, dynamic=True):
for n, run_name in enumerate(star_info):
    dead = [] ## Dead Pop III stars
    outputs = np.array([o for o in star_info[run_name].keys()])

    stellar_energy[run_name] = {}

    for j, o in enumerate(star_info[run_name]):

    #run_name = str(ds.directory).split('/')[-2]
    #o = str(ds).split('output_')[-1]

        time = DD_data[run_name][o]['time']
        redshift = DD_data[run_name][o]['redshift']

        print(o)
        print(run_name)

        stellar_energy[run_name][o] = {} 
        stellar_energy[run_name][o]['P2'] = {}
        stellar_energy[run_name][o]['P3'] = {}
        
        ## Halo info
        #if run_name == 'run_original':
        #    halo_rvir = ds.quan(massive_prog[o]['rvir'], 'unitary').to('kpc')
        #    halo_pos = ds.arr(massive_prog[o]['position'], 'unitary').to('kpc')

        #else:
        #    halo_rvir = ds.quan(run_halos[run_name][o]['rvir'], 'unitary').to('kpc')
        #    halo_pos = ds.arr(run_halos[run_name][o]['position'], 'unitary').to('kpc')
        
        #sp = ds.sphere(halo_pos, halo_rvir)

        ## Calculating the stellar radiative energies directly from the fields (for all radiative particles)

        #N_H = sp['H_mass'] / mass_hydrogen_cgs
        #N_He = sp['He_mass'] / (mass_hydrogen_cgs*4)
        #N_HeI = (sp['HeI_Density'] * sp['cell_volume']).to('g') / (mass_hydrogen_cgs*4)

        #N = [N_H, N_He, N_HeI]

        #for i, field in enumerate(kph_fields):
        #    kph = sp[field].to('1/s') * N[i] 
        #    #total_kph = sp.sum(field).to('1/s')
        #    total_lum = kph * eng[i]
        #    total_lum = np.sum(total_lum.to('erg/s'))
        #    print(field)
        #    print(total_lum)
        #    print('\n')

        ## Calculating stellar radiative energies from P2 stars
        #ds.add_particle_filter('p2')
        #p2_mass = sp[('p2', 'particle_mass')].to('Msun')
        #p2_age = sp[('p2', 'age')].to('Myr')
        P2 = star_info[run_name][o]['star_info']['P2']

        if len(P2) == 0: ## No P2 stars
            pass
        elif len(P2) == 1: ## No P2 stars
            pass
        else:
            p2_mass = P2['mass']
            p2_mass = yt.YTArray(p2_mass, 'Msun')
            time_list = np.zeros(len(p2_mass)) + time
            
            p2_age = (time_list - P2['creation_time']) * Myr


            Q_0 = p2_mass.v * StarClusterIonizingLuminosity
            Q_1 = EnergyFractionHeI * Q_0
            Q_2 = EnergyFractionHeII * Q_0
            Q_0 = Q_0 * (1 - EnergyFractionHeI - EnergyFractionHeII)

            E_0 = (Q_0 * eng[0]).to('erg') * 1/s
            E_1 = (Q_1 * eng[1]).to('erg') * 1/s
            E_2 = (Q_2 * eng[2]).to('erg') * 1/s

            total_E_0 = p2_age * E_0 #erg
            total_E_1 = p2_age * E_1
            total_E_2 = p2_age * E_2

            P2_radiative = [total_E_0.v.tolist(), total_E_1.v.tolist(), total_E_1.v.tolist()]
            stellar_energy[run_name][o]['P2']['radiative'] = P2_radiative

            ## Calculating SN energy (thermal) from P2 stars
            sn_eng = p2_mass * StarClusterSNEnergy
            therm = []
            for i, age in enumerate(p2_age):
                if age.v > 20:
                    therm.append(sn_eng[i].v.tolist())
                elif age.v < 4:
                    pass
                else:
                    time_pass = age - 4*Myr

                    frac = time_pass / (16 * Myr)
                    frac_eng = sn_eng[i]*frac

                    therm.append(frac_eng.v.tolist())
            therm = therm * erg

            stellar_energy[run_name][o]['P2']['thermal'] = therm.v.tolist()

        ## Calculating stellar radiative energies from P3 stars
        P3 = star_info[run_name][o]['star_info']['P3']
        if len(P3) == 0: ## No P3 stars
            pass
        else:
            p3_mass = P3['mass']
            p3_mass = yt.YTArray(p3_mass, 'Msun')
            p3_id = P3['id']

            Q0_p3 = [] # photon/s
            Q1_p3 = []
            Q2_p3 = []
            Q3_p3 = []

            for mass in p3_mass:
                x = np.log10(mass)
                x2 = x**2

                if (9 < mass.v <= 500):
                    Q0 = 10**(43.61 + 4.9 * x - 0.83*x2)
                    Q1 = 10**(42.51 + 5.69 * x - 1.01*x2)
                    Q2 = 10**(26.71 + 18.14 * x - 3.58*x2)
                    Q3 = 10**(44.03 + 4.59 * x - 0.77*x2)

                    Q0_p3.append(Q0)
                    Q1_p3.append(Q1)
                    Q2_p3.append(Q2)
                    Q3_p3.append(Q3)
                elif (5 < mass.v <= 9):
                    Q0 = 10**(39.29 + 8.55 * x)
                    Q1 = 10**(29.24 + 18.49 * x)
                    Q2 = 10**(26.71 + 18.14 * x - 3.58*x2)
                    Q3 = 10**(44.03 + 4.59 * x - 0.77*x2)

                    Q0_p3.append(Q0)
                    Q1_p3.append(Q1)
                    Q2_p3.append(Q2)
                    Q3_p3.append(Q3)
                else:
                    Q0_p3.append(0)
                    Q1_p3.append(0)
                    Q2_p3.append(0)
                    Q3_p3.append(0)

            E0_p3 = (Q0_p3 * eng_P3[0]).to('erg')
            E1_p3 = (Q1_p3 * eng_P3[1]).to('erg')
            E2_p3 = (Q2_p3 * eng_P3[2]).to('erg')
            E3_p3 = (Q3_p3 * eng_P3[3]).to('erg')

            P3_radiative = [E0_p3.v.tolist(), E1_p3.v.tolist(), E2_p3.v.tolist(), E3_p3.v.tolist()]
            stellar_energy[run_name][o]['P3']['radiative'] = P3_radiative

            ## Calculating thermal energy from P3 stars
            P3_SN = []
            for k, m in enumerate(p3_mass):
                if p3_id[k] in dead:
                    pass
                else:
                    if m < 1e-10:
                         ## Check previous output to see if mass was greater than 1

                        previous_output_index = np.argwhere(outputs == o)[0][0] - 1
                        previous_output = outputs[previous_output_index]

                        prev_masses = star_info[run_name][previous_output]['star_info']['P3']['mass']
                        prev_ids = star_info[run_name][previous_output]['star_info']['P3']['id']

                        prev_index = np.argwhere(np.array(prev_ids) == p3_id[k])
                        #print(prev_index)

                        if len(prev_index) == 0: ## Case where a P3 star has fallen into the main halo but went SN outside
                            print('This Pop III star fell in')
                            pass
                        else:
                            prev_mass = prev_masses[prev_index[0][0]]
                            print(prev_mass)
                            if prev_mass > 1:
                                ## Died this output! Calculate SN energy
                                dead.append(p3_id[k])

                                if (TypeIILowerMass <= prev_mass <= TypeIIUpperMass):
                                    if prev_mass < 20:
                                        P3_SN.append(SN_type_II.v.tolist())
                                    else:
                                        P3_SN.append(hypernova(prev_mass).v.tolist())
                                elif (PISNLowerMass <= prev_mass <= PISNUpperMass):
                                    P3_SN.append(PISN(prev_mass).v.tolist())
            stellar_energy[run_name][o]['P3']['thermal'] = P3_SN

with open('stellar_energy.json', 'w') as outfile:
    json.dump(stellar_energy, outfile)



        














