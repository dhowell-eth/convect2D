## Quick script to run some test cases of the finite Prandtl number code
## by varying the inputs file

import subprocess
import os
import shutil
import time
#Define Input files
params_file = 'convect_inputs.txt'
out_files=['T_field.dat','phi_field.dat','w_field.dat']
fortran_program='./a.out'

# Define test cases
test_cases = [
    {'ra':'1e5','total_time':0.1,'pr':10},
    {'ra':'1e5','total_time':0.1,'pr':1},
    {'ra':'1e5','total_time':0.1,'pr':0.1},
    {'ra':'1e5','total_time':0.1,'pr':0.01},
    {'ra':'1e5','total_time':1.0,'pr':0.01},
    {'ra': '1e7', 'total_time':0.1,'pr': 0.1}]

# Read in current parameters file
with open(params_file) as f:
    param_lines = [x.rstrip('\r\n') for x in f.readlines()]
original_params = param_lines

# Loop through test cases
for j,test in enumerate(test_cases):
    # Update parameters text
    for c in test.items():
        for i,x in enumerate(param_lines):
            if c[0].lower() == x.strip().split('=')[0].lower():
                param_lines[i] = '\t\t'+str(c[0])+'='+str(c[1])
    # Write new parameters file
    with open(params_file,'w') as o:
        for line in param_lines:
            print(line,file=o)

    # Call fortran program
    subprocess.call(('./a.out', ''))
    # Move output files to subdirectories
    if not os.path.exists(str(j)):
        os.mkdir(str(j))
    for x in out_files:
        shutil.move(x, str(j))

