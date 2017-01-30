import sys
import time

if len(sys.argv) != 3:
    print('Wrong number of arguments')
    sys.exit('Error')

import preprocess as pp
import denoise as de

noisy_file_split = sys.argv[1].split('/')
path_noisy_file = ''.join(noisy_file_split[:-1]) + '/'
file_name = noisy_file_split[-1]
denoised_file = sys.argv[2]
        

################################################################################
print(time.asctime(), ' Getting parameters')
num_reads, pos_min, pos_max = pp.get_params(sys.argv[1])
print(time.asctime(), ' Estimating sequence')
seq_est, len_seq = pp.estimate_sequence(
    sys.argv[1], num_reads, pos_min, pos_max)
print(time.asctime(), ' Estimating quality channel')
channel = pp.estimate_quality_channel(sys.argv[1], pos_min, seq_est)
print(time.asctime(), ' Counting contexts')
contexts = pp.count_quality_contexts(sys.argv[1], pos_min, seq_est, len_seq)
print(time.asctime(), ' Denoising')
de.quality_dude(
    sys.argv[1], denoised_file, pos_min, seq_est, len_seq, channel, contexts)
print(time.asctime(), ' Done denoising!')