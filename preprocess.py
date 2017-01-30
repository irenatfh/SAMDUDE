import numpy as np
import bisect
import subfunctions as sf


################################################################################
def get_params(
    file):
    """Get file parameters that will be used for other functions.
    
    Args:
        file (str): Path to the SAM file.
    
    Returns:
        num_reads (int): Total number of reads.
        pos_min (int): Minimum mapping position in the SAM file.
        pos_max (int): Maximum mapping position, after adjustment for clipping, 
        insertions and deletions.
                       
    """
    
    counter = 0
    with open(file, 'r') as f:
        for line in f:
            if line[0] != '@':
                data = line.split('\t')
                pos_map = int(data[3])
                cigar = data[5]
                read = data[9]
                if (cigar != '*'):
                    len_read = len(read)
                    if cigar != str(len_read) + 'M':
                        list_pos = sf.map_pos(pos_map, cigar, len_read)
                    else:
                        list_pos = list(range(pos_map, pos_map + len_read))
                    if counter == 0:
                        pos_min = pos_map
                        pos_max = list_pos[-1]
                    pos_max = max(list_pos[-1], pos_max)
                    counter += 1
    num_reads = counter - 1
    return num_reads, pos_min, pos_max


################################################################################
def estimate_sequence(
    file, num_reads, pos_min, pos_max,
    THRESHOLD=0.9, ALPH='ATGC'):
    """Estimate the sequence by taking the majority (defined by THRESHOLD) base
    over a pileup, for all pileups produced by the reads.
       
    Args:
        file (str): Path to the SAM file.
        num_reads (int): Total number of reads.
        pos_min (int): Minimum mapping position in the SAM file.
        pos_max (int): Maximum mapping position, after adjustment for clipping, 
        insertions and deletions.
        THRESHOLD (float): Threshold for majority voting (default 0.9)
        ALPH (str): The genomic alphabet (default 'ATGC').
           
    Returns:
        seq_est (numpy array): Sequence estimate encoded as integers.
        len_seq (int): Length of the sequence estimate.
        
     """
           
    seq_est = np.zeros(pos_max - pos_min + 1, dtype=np.int)
    counter = 0
    with open(file, 'r') as f:
        for line in f:
            if line[0] != '@':
                data = line.split('\t')
                pos_map = int(data[3])
                cigar = data[5]
                read = data[9]
                len_read = len(read)
                if (cigar != '*'):
                    if cigar != str(len_read) + 'M':
                        read_mapped, list_pos = sf.map_read(
                            pos_map, cigar, read, len_read)
                    else:
                        read_mapped = read
                        list_pos = list(range(pos_map, pos_map + len_read))
                    if counter == 0:
                        pos_prev = pos_map
                        list_prev = list_pos
                        len_pileups = list_pos[-1] - pos_map + 1
                        pileups = np.zeros((len_pileups, 4), dtype=np.int)
                    extra = list_pos[-1] - (pos_prev + len(pileups) - 1)
                    if extra > 0:
                        pileups = np.concatenate(
                            (pileups, np.zeros(
                                    (extra, 4), dtype=np.int)
                            ), axis=0)
                    for base, pos in zip(read_mapped, list_pos):
                        if base in ALPH:
                            pos_base = ALPH.find(base)
                            pileups[pos - pos_prev][pos_base] += 1
                    diff = pos_map - pos_prev
                    if (counter == num_reads) or (diff > 0):
                        pos_start = pos_prev - pos_min
                        if counter == num_reads:
                            diff = len(pileups)
                            pos_stop = pos_start + len(pileups)
                        else:
                            pos_stop = pos_map - pos_min
                        a = np.argmax(pileups[:diff], axis=1)
                        b = np.max(pileups[:diff], axis=1)
                        c = THRESHOLD * np.sum(pileups[:diff], axis=1)
                        a[b < c] = -1
                        seq_est[pos_start:pos_stop] = a
                        pos_last = pos_prev + len(pileups)
                        if pos_map > pos_last:  # no info about these positions
                            seq_est[pos_last:pos_map] = -1
                        pileups = pileups[diff:]
                    pos_prev = pos_map
                    list_prev = list_pos
                    counter += 1
        len_seq = len(seq_est)
    return seq_est, len_seq


################################################################################
def estimate_quality_channel(
    file, pos_min, seq_est, 
    ALPH='ATGC', BIN_LIMITS = [2, 10, 20, 25, 30, 35, 40]):
    """Estimate a channel that takes as inputs nucleotides and outputs
    (nucleotide, quality score bin) pairs. 
    
    Args:
        file (str): Path to the SAM file.
        pos_min (int): Minimum mapping position in the SAM file.
        seq_est (numpy array): Sequence estimate encoded as integers.
        ALPH (str): The genomic alphabet (default 'ATGC').
        BIN_LIMITS (list): List of quality score bin limits.
       
    Returns:
        channel (numpy array): (4 x 4 * number of bins) Probability transition
        matrices for each read position. Format for xth matrix row is: 
        [p(A, bin = 0|x) p(T, bin = 0|x) p(G, bin = 0|x) p(C, bin = 0|x)...
         p(A, bin = n|x) ... p(C, bin = n|x) ]    
    
    """
    num_bins = len(BIN_LIMITS) + 1
    channel = np.zeros((4, 4 * num_bins))
    with open(file, 'r') as f:
        for line in f:
            if line[0] != '@':
                data = line.split('\t')
                pos_map = int(data[3])
                cigar = data[5]
                read = data[9]
                len_read = len(read)
                scores = data[10]
                if (cigar != '*'):
                    if cigar != str(len_read) + 'M':
                        read_mapped, scores_mapped, list_pos = \
                        sf.map_read_scores(
                            pos_map, cigar, read, scores, len_read)
                    else:
                        read_mapped = read
                        scores_mapped = scores
                        list_pos = list(range(pos_map, pos_map + len_read))
                    scores_mapped = np.array([ord(i) 
                                              for i in scores_mapped]) - 33
                    for base, pos_seq, score in zip(
                        read_mapped, list_pos, scores_mapped):
                        pos_seq -= pos_min
                        if base in ALPH:
                            if seq_est[pos_seq] >= 0:
                                pos_base = ALPH.find(base)
                                base_maj = seq_est[pos_seq]
                                index_bin = bisect.bisect(BIN_LIMITS, score)
                                index_col = pos_base + 4 * index_bin
                                channel[base_maj][index_col] += 1
    channel = channel / channel.sum(axis=1, keepdims=True)
    return channel


################################################################################
def count_quality_contexts(
    file, pos_min, seq_est, len_seq, 
    k=10, ALPH='ATGC', BIN_LIMITS = [2, 10, 20, 25, 30, 35, 40]):
    """Create a dictionary of context counts, where the key is the string of
    length 2k surrounding the central symbol, and the values are counts of the
    central symbol as a (nucleotide, quality score bin) pair.
    
    Args:
        file (str): Path to the SAM file.
        pos_min (int): Minimum mapping position in the SAM file.
        seq_est (numpy array): Sequence estimate encoded as integers.
        len_seq (int): Length of the sequence estimate.
        k (int): One-sided context length, default 10.
        ALPH (str): The genomic alphabet (default 'ATGC').
        BIN_LIMITS (list): List of quality score bin limits.
       
    Returns:
        contexts (dict): Dictionary of contexts, where keys are double-sided 
        contexts (strings) and values are a lists (of length 4 * number of bins)
        of counts.
    
    """
    num_bins = len(BIN_LIMITS) + 1
    contexts = {}
    with open(file, 'r') as f:
        for line in f:
            if line[0] != '@':
                data = line.split('\t')
                pos_map = int(data[3])
                cigar = data[5]
                read = data[9]
                len_read = len(read)
                scores = data[10]
                if (cigar != '*'):
                    read_contextualized, scores_contextualized = \
                    sf.contextualize_read(
                        pos_map, cigar, read, scores, pos_min, len_read,
                        seq_est, len_seq, k)
                    invalid_context = 2*k + 1
                    for ind in range(len(read_contextualized)):
                        base_incoming = read_contextualized[ind]
                        if base_incoming not in ALPH:
                            invalid_context = 2*k + 1
                        else:
                            invalid_context -= 1
                        if invalid_context <= 0:
                            lhs = read_contextualized[(ind - 2*k):(ind - k)]
                            rhs = read_contextualized[(ind - k + 1):(ind + 1)]
                            context = lhs + rhs
                            base = read_contextualized[ind - k]
                            score = scores_contextualized[ind - k]
                            pos_base = ALPH.find(base)
                            index_bin = bisect.bisect(BIN_LIMITS, score)
                            index_col = pos_base + 4 * index_bin
                            if context not in contexts:
                                contexts[context] = [0] * 4 * num_bins
                            contexts[context][index_col] += 1
    return contexts