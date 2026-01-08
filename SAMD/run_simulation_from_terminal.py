from subprocess import PIPE, run
from pathlib import Path

# Set paths
script_path = Path(__file__).resolve() # get path of current script
script_dir = script_path.parent # get directory
file_name = 'SAMD_simulator.py' # file name
full_path = script_dir / file_name # build path 

# Run script
output = run("abaqus cae noGUI="+str(full_path), stdout=PIPE, stderr=PIPE, universal_newlines=True, shell=True)
print(output)
