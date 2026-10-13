README for Crystal Clear: A shortcut to peak picking and indexing in powder x-ray diffraction (PXRD)

Dependencies

This program requires Python 3 and imported modules/libraries, including NumPy and Matplotlib.
These can be installed by using:
"conda install numpy matplotlib" OR "pip install numpy matplotlib"

The program also requires additional files within the same working directory to function, including:
crystal_clear.py (the main code/script)
peak_picking_functions.py (python file containing peak picking related functions)
experimental data file (e.g HA_experimental.csv)
reference data file (e.g HA_standard.csv)
ICDD database file (e.g HA_ICDD.txt)

User code preparation

Before running the code, ensure that the required test files are in their required format, before calling on the command line of the terminal.
Also ensure that the chemical symbol for the incident radiation is present in the experimental file (not case sensitive).
The advanced user settings should only be edited if the user requires more than what the default settings offer, which include:
intensity_threshold_value = 0.5 (peak detection sensitivity)
angle_tolerance = 0.5 (peak matching tolerance)
intensity_tolerance = 0.1 (background filtering)

Running the code

The script should be ran via the terminal, calling the files in the following order:
python crystal_clear.py experimental_file.csv reference_file.csv ICDD_file.txt 

Output 

If the program is supplied with correct input information it should:
Print the peak positions and associated d spacings above the intensity_threshold_value
Print miller indices for peaks within the experimental range
Show an overlayed PXRD plot and associated RMSD plot
Print the RMSD value and associated matched peak pairs





