# NSM_analysis
Analysis pipeline for NSMs in the early universe.

Directory structure:

Let's assume you have run a simulation that you now want to restart with a NSM. Pipeline runs in the following order:
1. `output_OG_info.py`
   This creates a filed called `DD_data_OG.json` which contains the time and redshift for each dataset.
2. `stellar_info.py`
   This script grabs the final output of the OG simulation and produces three things:
    1. `number_of_stars.txt`: Lists the number of living P3 stars, all P3 stars, and P2 stars.
    2. `images/P2_creation_times.png`: A histogram of the creation times of all P2 stars.
    3. `images/P3_creation_times.png`: A histogram of the creation times of all P3 stars. 
    The purpose is to use this information to help make a decision about which P3 star to choose to become the binary. The next script will also help with this. 
3. `potential_choices.py`:
    This script looks within the central part of the main halo for P3 stars that fall within the required mass range. It will also create a dictionary called `potential_p3_info.json` with information about these potential P3 stars. 
    Next, it will find the min and max values of four different fields in order to create movies for each potential P3 star. It writes these values to a dictionary called `min_max_WIDTH.json` where WIDTH is the width of the plot. 
    Finally, it iterates over the datasets where P3 stars exist and makes plots for each potential star. 
4. `potential_choices_movies.py`:
    Creates movies for each potential p3 star.

Once you use the previous information to choose your P3 sacrifice, create a file called `PopIII_params.txt` and include these four lines:
    
```
PopIII_NeutronStarMergers       = 1         # Turns NSMs on
PopIII_NSMParticleID            = 761313    # Chosen particle ID to become NSM
PopIII_NSMExplosionEnergy       = 1e+50     # Desired NSM explosion energy in ergs
PopIII_NSMDelayTime             = 100       # Time between P3_binary explosion and NSM in Myr
PopIII_NSMMetalMass             = 0.01      # Amount of metals ejected in solar masses
```
Run the NSM simulations. Once they are done, run the following scripts:
5. `all_output_info.py`
    This does a similar thing as step 1, creating a dictionary called `run_DD_data.json` of the time and redshift for every dataset from each of your runs. This also writes a dictionary called `param_info.json` which has the parameters for each of your runs.

6. `find_halo_tree.py`
    This script will try to match the merger trees between the NSM run and the OG run up until the restart dataset. After which it will just find the largest halo that the NSM particle is in. It will produce a dictionary called `run_halos_incl_vir.json` which contains all the halo information from the tree. 

7. `NSM_dict.py`
   This script determines when the NSM takes place in each simulation. Produces `NSM_dict.json`.

8. `data_grab.py`
    This script grabs a bunch of information, designated in `load_ins.py` to grab data from the runs for easy plotting. 

9. `make_plots.py`
    Makes the profile, phase, and projection plots from the data grabbed in the previous step.

10. `make_movies.py` 
    Makes movies made from images produced in the previous step.

If your original simulation needs to run longer, run `output_OG_info.py`, `OG_most_massive.py`, `all_output_info.py` before you run step 6. 