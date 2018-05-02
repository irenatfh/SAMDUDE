import numpy as np
import re


################################################################################
def map_pos(
    pos_map, cigar, len_read):
    """Use the cigar string to get a list of mapping positions of the read to
    the reference genome. The positions will be processed as follows: inserted
    positions are removed from the read while deleted bases result in mapping   
    position adjustment, and soft-clipped bases are discarded.
    
    Args:
        pos_map (int): Mapping position of the read.
        cigar (str): SAM file CIGAR string.
        len_read (int): Read length.
       
    Returns:
        list_pos (int): List of mapping positions for each base in the processed
        read.
        
    """
    
    ns_all = re.findall('\d+', cigar)
    ops_all = re.findall('\D+', cigar)
    keep_indices = [i for i in range(len(ops_all)) if ops_all[i] is not 'H']
    ns = [int(ns_all[i]) for i in keep_indices]
    ops = [ops_all[i] for i in keep_indices]
    index_m = ops.index('M')
    padding = sum(ns[:index_m])
    pos_map -= padding
    list_pos = list(range(pos_map, pos_map + len_read))
    pointer = 0
    for (i, j) in enumerate(zip(ns, ops)):
        n = j[0]
        op = j[1]
        if op == 'S':
            if i == 0:
                list_pos = list_pos[n:]
            else:
                list_pos = list_pos[:-n]
        elif (op == 'M') or (op == '=') or (op == 'X'):
            pointer += n
        elif (op == 'D') or (op == 'N'):
            list_pos = [l + n if m >= pointer else l for (l, m) in 
                        zip(list_pos, range(len(list_pos)))]
        elif (op == 'I') or (op == 'P'):
            list_pos = list_pos[:-n]
    return list_pos


################################################################################
def map_read(
    pos_map, cigar, read, len_read):
    """Process the read to be mappable to the reference genome based on the
    CIGAR string. Inserted bases will be removed from the read while deleted
    bases will result in mapping position adjustment, and soft-clipped bases are discarded.
    
    Args:
        pos_map (int): Mapping position of the read.
        cigar (str): SAM file CIGAR string.
        read (str): Read comprising bases in AGCTN.
        len_read (int): Read length.
       
    Returns:
        read (str): Processed read string.
        list_pos (int): List of mapping positions for each base in the processed
        read.
        
    """
    
    ns_all = re.findall('\d+', cigar)
    ops_all = re.findall('\D+', cigar)
    keep_indices = [i for i in range(len(ops_all)) if ops_all[i] is not 'H']
    ns = [int(ns_all[i]) for i in keep_indices]
    ops = [ops_all[i] for i in keep_indices]
    index_m = ops.index('M')
    padding = sum(ns[:index_m])
    pos_map -= padding
    list_pos = list(range(pos_map, pos_map + len_read))
    pointer = 0
    for (i, j) in enumerate(zip(ns, ops)):
        n = j[0]
        op = j[1]
        if op == 'S':
            if i == 0:
                read = read[n:]
                list_pos = list_pos[n:]
            else:
                read = read[:-n]
                list_pos = list_pos[:-n]
        elif (op == 'M') or (op == '=') or (op == 'X'):
            pointer += n
        elif (op == 'D') or (op == 'N'):
            list_pos = [l + n if m >= pointer else l for (l, m) in 
                        zip(list_pos, range(len(list_pos)))]
        elif (op == 'I') or (op == 'P'):
            read = read[:pointer] + read[pointer + n:]
            list_pos = list_pos[:-n]
    return read, list_pos


################################################################################
def map_read_scores(
    pos_map, cigar, read, scores, len_read):
    """Process the read and quality score string to be mappable to the
    reference genome based on the SAM file CIGAR string. Inserted bases will be 
    removed from the read while deleted bases will result in mapping position 
    adjustment, and soft-clipped bases are discarded.
    
    Args:
        pos_map (int): Mapping position of the read.
        cigar (str): SAM file CIGAR string.
        read (str): Read comprising bases in AGCTN.
        scores (str): Quality score string.
        len_read (int): Read length.
       
    Returns:
        read (str): Processed read string.
        scores (str): Processed scores string.
        list_pos (int): List of mapping positions for each base in the processed
        read.
        
    """
    
    ns_all = re.findall('\d+', cigar)
    ops_all = re.findall('\D+', cigar)
    keep_indices = [i for i in range(len(ops_all)) if ops_all[i] is not 'H']
    ns = [int(ns_all[i]) for i in keep_indices]
    ops = [ops_all[i] for i in keep_indices]
    index_m = ops.index('M')
    padding = sum(ns[:index_m])
    pos_map -= padding
    list_pos = list(range(pos_map, pos_map + len_read))
    pointer = 0
    for (i, j) in enumerate(zip(ns, ops)):
        n = j[0]
        op = j[1]
        if op == 'S':
            if i == 0:
                read = read[n:]
                scores = scores[n:]
                list_pos = list_pos[n:]
            else:
                read = read[:-n]
                scores = scores[:-n]
                list_pos = list_pos[:-n]
        elif (op == 'M') or (op == '=') or (op == 'X'):
            pointer += n
        elif (op == 'D') or (op == 'N'):
            list_pos = [l + n if m >= pointer else l for (l, m) in 
                        zip(list_pos, range(len(list_pos)))]
        elif (op == 'I') or (op == 'P'):
            read = read[:pointer] + read[pointer + n:]
            scores = scores[:pointer] + scores[pointer + n:]
            list_pos = list_pos[:-n]
    return read, scores, list_pos


################################################################################
def contextualize_read(
    pos_map, cigar, read, scores, pos_min, len_read, seq_est, len_seq, k, 
    ALPH='ATGC'):
    """Process the read and quality score string in anticipation of obtaining
    contexts. Keep inserted portions of the read, adjust list position of
    deleted bases, and discard soft-clipped regions. Pad the read with a header 
    and footer from the sequence estimate if the read begins and ends,
    respectively, with matches.
    
    Args:
        pos_map (int): Mapping position of the read.
        cigar (str): SAM file CIGAR string.
        read (str): Read comprising bases in AGCTN.
        scores (str): Quality score string.
        pos_min (int): Minimum mapping position in the SAM file.
        len_read (int): Read length.
        seq_est (numpy array): Sequence estimate encoded as integers.
        len_seq (int): Length of the sequence estimate.
        k (int): One-sided context length.
       
    Returns:
        read (str): Processed read string.
        scores (np array): Processed quality score string converted to integers.
        
    """
    
    ns_all = re.findall('\d+', cigar)
    ops_all = re.findall('\D+', cigar)
    keep_indices = [i for i in range(len(ops_all)) if ops_all[i] is not 'H']
    ns = [int(ns_all[i]) for i in keep_indices]
    ops = [ops_all[i] for i in keep_indices]
    index_m = ops.index('M')
    padding = sum(ns[:index_m])
    pos_map -= padding
    list_pos = list(range(pos_map, pos_map + len_read))
    pointer = 0
    for (i, j) in enumerate(zip(ns, ops)):
        n = j[0]
        op = j[1]
        if op == 'S':
            if i == 0:
                read = read[n:]
                scores = scores[n:]
                list_pos = list_pos[n:]
            elif i == len(ns) - 1:
                read = read[:-n]
                scores = scores[:-n]
                list_pos = list_pos[:-n]
        elif (op == 'M') or (op == '=') or (op == 'X'):
            pointer += n
        elif (op == 'D') or (op == 'N'):
            list_pos = [l + n if m >= pointer else l for (l, m) in 
                        zip(list_pos, range(len(list_pos)))]
        elif (op == 'I') or (op == 'P'):
            list_pos = list_pos[:pointer] + [''] * n + list_pos[pointer:-n]
            pointer += n
    if type(list_pos[0]) is int:
        start_pos = list_pos[0] - k - pos_min
        end_pos = start_pos + k
        head = seq_est[max(0, start_pos):max(0, end_pos)]
        head = [ALPH[i] if i in range(4) else 'N' for i in head]
        head = ''.join(head)
        read = head + read
        scores = 'z' * len(head) + scores
    if type(list_pos[-1]) is int:
        start_pos = list_pos[-1] + 1 - pos_min
        end_pos = start_pos + k
        tail = seq_est[min(len_seq, start_pos):min(len_seq, end_pos)]
        tail = [ALPH[i] if i in range(4) else 'N' for i in tail]
        tail = ''.join(tail)
        read += tail
        scores += 'z' * len(tail)
    scores = np.array([ord(i) for i in scores]) - 33
    return read, scores


################################################################################
def prepare_read(
    pos_map, cigar, read, scores, len_read, pos_min, seq_est, len_seq, k, 
    ALPH='ATGC'):
    """Process the read and quality score string in anticipation of
    denoising. Keep all bases in the read, pad the read with a header only if it
    starts with mappable bases (i.e., 'nM' for any integer n), and pad the read 
    with a footer only if it ends with mappable bases.
    
    Args:
        pos_map (int): Mapping position of the read.
        cigar (str): SAM file CIGAR string.
        read (str): Read comprising bases in AGCTN.
        scores (str): Quality score string.
        len_read: Read length.
        pos_min (int): Minimum mapping position in the SAM file.
        seq_est (numpy array): Sequence estimate encoded as integers.
        len_seq (int): Length of the sequence estimate.
        k (int): One-sided context length.
        ALPH (str): The genomic alphabet (default 'ATGC').
       
    Returns:
        read_prepared (str): Processed read string.
        scores_prepared (np array): Processed quality score string converted to 
        integers.
        list_ind (list): Processed list of positions within the read.
        
    """

    list_ind = list(range(len_read))
    if (cigar != '*') and ('M' in cigar):
        ns_all = re.findall('\d+', cigar)
        ops_all = re.findall('\D+', cigar)
        keep_indices = [i for i in range(len(ops_all)) if ops_all[i] is not 'H']
        ns = [int(ns_all[i]) for i in keep_indices]
        ops = [ops_all[i] for i in keep_indices]
        index_m = ops.index('M')
        padding = sum(ns[:index_m])
        pos_map -= padding
        list_pos = list(range(pos_map, pos_map + len_read))
        pointer = 0
        for (i, j) in enumerate(zip(ns, ops)):
            n = j[0]
            op = j[1]
            if (op == 'M') or (op == '=') or (op == 'X') or (op == 'S'):
                pointer += n
            elif (op == 'D') or (op == 'N'):
                list_pos = [l + n if m >= pointer else l for (l, m) in 
                            zip(list_pos, range(len(list_pos)))]
            elif (op == 'I') or (op == 'P'):
                list_pos = list_pos[:pointer] + [''] * n + list_pos[pointer:-n]
                pointer += n
        if ops[0] == 'M':
            start_pos = list_pos[0] - k - pos_min
            end_pos = start_pos + k
            head = seq_est[max(0, start_pos):max(0, end_pos)]
            head = [ALPH[i] if i in range(4) else 'N' for i in head]
            head = ''.join(head)
            read = head + read
            scores = 'z' * len(head) + scores
            list_ind = [''] * len(head) + list_ind
        if ops[-1] == 'M':
            start_pos = list_pos[-1] + 1 - pos_min
            end_pos = start_pos + k
            tail = seq_est[min(len_seq, start_pos):min(len_seq, end_pos)]
            tail = [ALPH[i] if i in range(4) else 'N' for i in tail]
            tail = ''.join(tail)
            read += tail
            scores += 'z' * len(tail)
            list_ind += [''] * len(tail)
    scores = np.array([ord(i) for i in scores]) - 33
    return read, scores, list_ind
