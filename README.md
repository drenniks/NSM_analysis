# NSM_analysis
Analysis pipeline for NSMs in the early universe.

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
5. `5_NSM_dict.py`
    This script will create a dictionary `NSM_dict.json` to indicate when the chosen P3 star became a P3_binary, a NS_binary, and then a bh. The output and time are given for each particle type.
