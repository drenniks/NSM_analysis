# NSM_analysis
Analysis pipeline for NSMs in the early universe.

Pipeline runs in the following order: 
1. `1_output_dict.py`
    This creates a file called `DD_data.json` which contains information relevant to each dataset (specifically the time and redshift for each output). It iss used in subsequent analysis.
    
2. `2_stellar_creations.py`
    This script takes the final output and produces three things:
    1. `number_of_stars.txt`: Lists the number of living P3 stars, all P3 stars, and P2 stars.
    2. `images/P2_creation_times.png`: A histogram of the creation times of all P2 stars.
    3. `images/P3_creation_times.png`: A histogram of the creation times of all P3 stars. 
    The purpose is to use the creation times of each group of stars to help make a decision about which P3 star to choose to become the binary. The next script will aid in this process.
    
3. `3_potential_choices.py`:
    This script looks within the central part of the main halo for P3 stars that fall within the required mass range. It creates a file called `potential_p3_info.txt` with information about the P3 stars that fall in the required mass range. It will also create a dictionary called `potential_p3_info.json` with the same information. 
    Next, it will find the min and max values of four different fields in order to create movies for each potential P3 star. It writes these values to a dictionary called `min_max_WIDTH.json` where WIDTH is the width of the plot. 
    Finally, it iterates over the datasets where P3 stars exist and makes plots for each potential star. 
4. `4_potential_movies.py`:
    Creates movies for each potential p3 star.

Once you use the previous information to choose your P3 sacrifice, create a file called `PopIII_params.txt` and include these four lines:
    ```
    PopIII_NeutronStarMergers       = 1         # Turns NSMs on
    PopIII_NSMParticleID            = 761313    # Chosen particle ID to become NSM
    PopIII_NSMExplosionEnergy       = 1e+50     # Desired NSM explosion energy in ergs
    PopIII_NSMDelayTime             = 100       # Time between P3_binary explosion and NSM in Myr
    PopIII_NSMMetalMass             = 0.01      # Amount of metals ejected in solar masses
    ```
