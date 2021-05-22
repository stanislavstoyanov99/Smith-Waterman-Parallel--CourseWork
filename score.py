def score_function(row_index, col_index, score_matrix, other, seq1, seq2, match):
    maxScore = 0
    # print('\nWorker: %s' % num)
    
    similarity = match if seq1[row_index - 1] == seq2[col_index - 1] else other
    diag_score = score_matrix[row_index - 1][col_index - 1] + similarity
    up_score = score_matrix[row_index - 1][col_index] + other
    left_score = score_matrix[row_index][col_index - 1] + other
    curMax = max(0, diag_score, up_score, left_score)

    if curMax > maxScore:
        maxScore = curMax
        maxPosition = (row_index, col_index)

    score_matrix[row_index][col_index] = curMax

    return