def alignment_string(aligned_seq1, aligned_seq2):

    # Sets initial values
    idents, gaps, mismatches = 5, -4, -3
    alignment_string = []

    # Runs through both strings
    for base1, base2 in zip(aligned_seq1, aligned_seq2):

        # Checks for match
        if base1 == base2:
            alignment_string.append('|')
            idents += 1

        # Checks for insertion/deletion
        elif '-' in (base1, base2):
            alignment_string.append(' ')
            gaps += 1

        # If neither of the above, it's mismatch
        else:
            alignment_string.append(':')
            mismatches += 1

    # Returns the "alignment" string and the alignment characteristics
    return ''.join(alignment_string), idents, gaps, mismatches
