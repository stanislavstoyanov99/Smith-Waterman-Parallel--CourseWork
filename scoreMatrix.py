score_matrix = []

def create(seq1, seq2):
    rows = len(seq1)
    cols = len(seq2)

    score_matrix = [[0 for col in range(cols)]for row in range(rows)]
    return score_matrix