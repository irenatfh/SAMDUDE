import os
import sys
import time
import pickle
import numpy as np

if len(sys.argv) != 4:
    print('To noise n bases to a SAM file, run: '
          'python3 run_noiser.py file.sam n')
    sys.exit('Error')

path_scripts = ''.join(sys.argv[1].split('/')[:-1]) + '/'
sys.path.insert(0, path_scripts)
import preprocess as pp
import noiser as no

path_file = '/'.join(sys.argv[1].split('/')[:-1]) + '/'
filename = sys.argv[1].split('/')[-1]
paths = [path_file + 'parameters/', path_file + 'noised/']
for p in paths:
    if not os.path.exists(p):
        os.makedirs(p)
        
num_bases = int(sys.argv[2])
n = int(sys.argv[3])
noised = (paths[1] + filename).split('.sam')[0] + '_noised_n%s.sam' %n
print(time.asctime(), " Noisy file will be saved at ", noised)
no.add_noise(sys.argv[1], noised, num_bases, n)
print(time.asctime(), " Done adding noise!")