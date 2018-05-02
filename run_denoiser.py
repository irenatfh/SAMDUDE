import sys
import time
import timeit
import resource

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
start = timeit.default_timer()
print(time.asctime(), ' Getting parameters')
num_reads, pos_min, pos_max = pp.get_params(sys.argv[1])
print(time.asctime(), ' Estimating sequence')
seq_est, len_seq = pp.estimate_sequence(
    sys.argv[1], num_reads, pos_min, pos_max)
print(time.asctime(), ' Estimating quality channel')
channel = pp.estimate_quality_channel(sys.argv[1], pos_min, seq_est)
print('Quality channel:\n', channel)
print(time.asctime(), ' Counting contexts')
contexts = pp.count_quality_contexts(sys.argv[1], pos_min, seq_est, len_seq)
print(time.asctime(), ' Denoising')
de.quality_dude(
    sys.argv[1], denoised_file, pos_min, seq_est, len_seq, channel, contexts)
stop = timeit.default_timer()
print(time.asctime(), ' Done denoising!')
print('Maximum memory used: %.3f GB' % (resource.getrusage(
        resource.RUSAGE_SELF).ru_maxrss / 1000000))

total_time = stop - start
print('Total running time: %.3f' % total_time)


################################################################################
import argparse
parser=argparse.ArgumentParser(
    description='SAMDUDE is a genomic sequence denoiser that operates on aligned SAM files')
parser.add_argument('noisy_file', type=str, help='Filename (including path) of the SAM file to be denoised.')
parser.add_argument('denoised_file', type=str, help='Filename (including path) of the denoised SAM file.')
