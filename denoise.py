import numpy as np
import bisect
import subfunctions as sf
                    

################################################################################
def quality_dude(
    file_noisy, file_denoised, pos_min, seq_est, len_seq, channel, contexts, 
    k=7, ALPH='ATGC', BIN_LIMITS = [2, 10, 20, 25, 30, 35, 40]):
    """
    DUDE-based denoising algorithm. Input data is assumed to be mapped reads
    (and their accompanying quality score strings) taken from sorted SAM files.
    The algorithm both denoises the reads and updates the quality score string
    when necessary, i.e.: if the algorithm chooses a base that is different from
    the original base, it replaces the quality score with one derived from its
    confidence calculation of the new base. If the algorithm chooses the same 
    base, it replaces the quality score with one derived from the average of the
    probabilities of the original basecall and its own confidence calculation.
    In order to perform denoising, the algorithm assumes that the sequencer
    behaves as a memoryless noise channel which is assumed to take as inputs 
    bases A, T, G and C, but to output quality (nucleotide, quality score) 
    tuples.
    
    Quality scores are binned as per Illumina's quality score 8-level mapping
    (2014). The channel is a 4 x (4 x 8) matrix, where columns correspond to
    (nucleotide, bin) pairs. The channel was populated with counts from
    positions with a clear majority. Contexts were recorded from all positions
    (except for clipped regions).
    
    Bases are denoised if the following criteria are met:
        1. The read has an alignment position (denoising is attempted on all
        reads with a left-most mapping position, even if the cigar string = *).
        2. The context and base comprise only A, T, G and C bases.
        3. The probability of the basecall, calculated from the quality score,
        is < 0.9.
        
    During denoising, bases are denoised based on the double-sided context of 
    length 2*k flanking the base. The context may contain clipped regions and 
    bases in the read marked as insertions relative to the reference genome. 
    When the read ends in mappable bases, the read is further augmented with 
    bases from the sequence estimate in order to aid context acquisition.
    
    Args:
        file_noisy (SAM file): Noisy file
        file_denoised (SAM file): Denoised file where quality scores are
        replaced wherever DUDE chooses a different base, and are replaced with a
        quality score based on the average of the probability corresponding to
        the original quality score and the DUDE conditional probability
        estimate.
        pos_min (int): Minimum mapping position in the SAM file.
        seq_est (numpy array): Sequence estimate encoded as integers.
        len_seq (int): Length of sequence estimate
        channel (numpy array): Channel transition matrix
        k (int): One-sided context length, default 10.
        ALPH (str): The genomic alphabet (default 'ATGC').
        BIN_LIMITS (list): List of quality score bin limits.
       
    Returns:
        A denoised file under path/filename file_denoised
        
    """
    pi_transinv = np.dot(
        np.linalg.inv(
            np.dot(channel, 
                   np.transpose(channel)
                  ),
            ),
        channel
    )
    hist_original = np.zeros(50, dtype=int)
    hist_qsa = np.zeros(50, dtype=int)
    counter = 0
    with open(file_noisy, 'r') as fn, open(file_denoised, 'w') as fd:
        for line in fn:
            if line[0] != '@':
                data = line.split('\t')
                pos_map = int(data[3])
                cigar = data[5]
                read = data[9]
                len_read = len(read)
                scores = data[10]
                read_denoised = list(read)
                scores_denoised = list(scores)
                read_prepared, scores_prepared, list_ind = sf.prepare_read(
                    pos_map, cigar, read, scores, len_read, 
                    pos_min, seq_est, len_seq, k)
                probs = 1 - 10 ** (-scores_prepared / 10)
                for ind in range(k, len(read_prepared) - k):
                    base = read_prepared[ind]
                    if base in ALPH:
                        score = scores_prepared[ind]
                        hist_original[score] += 1
                        prob = probs[ind]
                        lhs = read_prepared[ind - k:ind]
                        rhs = read_prepared[ind + 1:ind + k + 1]
                        context = lhs + rhs
                        if ((context in contexts) and (prob < 0.9)):
                            pos_base = ALPH.find(base)
                            index_bin = bisect.bisect(BIN_LIMITS, score)
                            index_col = pos_base + 4 * index_bin
                            counts_est = np.multiply(
                                channel[:, index_col],
                                np.dot(pi_transinv, contexts[context])
                            )
                            # check validity as a probability vector
                            if np.sum(counts_est) > 1e-5:
                                prob_est = counts_est / np.sum(counts_est)
                                base_opt = np.argmax(prob_est)
                                prob_opt = np.min((0.9999369, np.max(prob_est)))
                                if base_opt == ALPH.index(base):
                                    prob_opt = np.mean((prob, prob_opt))
                                else:
                                    counter += 1
                                score_opt = int(round(-10 * np.log10(
                                            1 - prob_opt)))
                                hist_qsa[score_opt] += 1
                                score_opt = chr(score_opt + 33)
                                ind_write = list_ind[ind]
                                read_denoised[ind_write] = ALPH[base_opt]
                                scores_denoised[ind_write] = score_opt
                data[9] = ''.join(read_denoised)
                data[10] = ''.join(scores_denoised)
                fd.write('\t'.join(data))
            else:
                fd.write(line)
    print('Bases changed: %s' %counter)
    print('Original quality scores:\n%s' %hist_original)
    print('QSA scores:\n%s' %hist_qsa)
    return
