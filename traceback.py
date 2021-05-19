def traceback(score_matrix, start_pos , seq1 , seq2):

    end,diag, up, left,  = range(4)

    aligned_seq1, aligned_seq2 = []
    x, y = start_pos
    move = nextMove(score_matrix, x, y)

    try:
        while move != end:
            if move == diag:
                aligned_seq1.append(seq1[x - 1])
                aligned_seq2.append(seq2[y - 1])
                x -= 1
                y -= 1
            elif move == up:
                aligned_seq1.append(seq1[x - 1])
                aligned_seq2.append('-')
                x -= 1
            elif move == left:
                aligned_seq1.append('-')
                aligned_seq2.append(seq2[y - 1])
                y -= 1
            else:
                move = end
            move = nextMove(score_matrix, x, y)
    except:
        move = end

    try:
        aligned_seq1.append(seq1[x - 1])
    except:
        aligned_seq1.append(seq1[x])
    try:
        aligned_seq2.append(seq1[y - 1])
    except:
        aligned_seq2.append(seq1[y])

    return ''.join(reversed(aligned_seq1)), ''.join(reversed(aligned_seq2))

def nextMove(score_matrix, x, y):

    # Assign the diagonal/insertion/deletion scores
    diag = score_matrix[x - 1][y - 1]
    up = score_matrix[x - 1][y]
    left = score_matrix[x][y - 1]

    # Check all three cases to find next character/insertion/deletion
    if diag >= up and diag >= left:
        return 1 if diag != 0 else 0
    elif up > diag and up >= left:
        return 2 if up != 0 else 0
    elif left > diag and left > up:
        return 3 if left != 0 else 0
    # Error detection
    else:
        # Execution should not reach here.
        raise ValueError('invalid move during traceback')
