import numpy as np
import bisect
np.set_printoptions(threshold=np.inf)
import sys
sys.path.insert(0, '/mnt/Irena/denoising/scripts/')
import subfunctions as sf
import random
import time
from collections import deque


################################################################################
def add_noise(file_o, file_n, num_bases, n,
              ALPH='ATGC'):
    """
    Script adds noise at a rate of n bases out of num_bases bases
    randomly. Note that selected bases will be changed to one of the other
    three possible bases with equal probability.
    
    Args:
        file_o (SAM file): Original file
        file_n (SAM file): Noisy file
        num_bases(int): Number of bases in file
        n (int): Number of bases to change.
        ALPH (str): The genomic alphabet (default 'ATGC').
       
    Returns:
        A noisy file under path/filename file_n
        
    """
    print(time.asctime(), " Determining noising locations")
    noise_locs = random.sample(range(num_bases), n)
    noise_locs.sort()
    noise_locs = deque(noise_locs)
    counter = 0
    print(time.asctime(), " Adding noise")
    with open(file_o, 'r') as fo, \
         open(file_n, 'w') as fn:
            for line in fo:
                if line[0] != '@':
                    data = line.split('\t')
                    read = data[9]
                    noised_read = list(read)
                    for ind in range(len(read)):
                        if len(noise_locs) > 0:
                            if counter == noise_locs[0]:
                                base = read[ind]
                                if base in ALPH:
                                    base_ind = ALPH.index(read[ind])
                                    base_list = list(
                                        {0, 1, 2, 3}.difference({base_ind})
                                    )
                                    noisy_base = np.random.choice(
                                        base_list, 1)[0]
                                else:
                                    noisy_base = np.random.choice(
                                        [0, 1, 2, 3], 1)[0]
                                noised_read[ind] = ALPH[noisy_base]
                                noise_locs.popleft()
                        counter += 1
                    data[9] = ''.join(noised_read)
                    fn.write('\t'.join(data))
                else:
                    fn.write(line)
    return