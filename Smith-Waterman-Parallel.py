import threading
import time
#import timeit
# from array import array

# Test Case 1
#seq1 = "ATAGACGACATGGGGACAGCATATAGACGACATGGGGACAGCATATAGACGACATGGGGACAGCATATAGACGACATGGGGACAGCATATAGACGACATGGGGACAGCATATAGACGACATGGGGACAGCATATAGACGACATGGGGACAGCATATAGACGACATGGGGACAGCAT"
#seq2 = "TTTAGCATGCGCATATCAGCAATTTAGCATGCGCATATCAGCAATTTAGCATGCGCATATCAGCAATTTAGCATGCGCATATCAGCAATTTAGCATGCGCATATCAGCAATTTAGCATGCGCATATCAGCAATTTAGCATGCGCATATCAGCAATTTAGCATGCGCATATCAGCAA"

# Test Case 2
# seq1 = "ATAGACGACATGGGGACAGCATATAGACGACATGGGGACAGCATATAGACGACATGGGGACAGCATATAGACGACATGGGGACAGCAT"
# seq2 = "TTTAGCATGCGCATATCAGCAATTTAGCATGCGCATATCAGCAATTTAGCATGCGCATATCAGCAATTTAGCATGCGCATATCAGCAA"

# Test Case 3
# seq1 = "ATAGACGACATGGGGACAGCATATAGACGACATGGGGACAGCAT"
# seq2 = "TTTAGCATGCGCATATCAGCAATTTAGCATGCGCATATCAGCAA"

# Test Case 4
seq1 = "ATAGACGACATGGGGACAGCAT"
seq2 = "TTTAGCATGCGCATATCAGCAA"

# Test Case 5
#seq1 = "ATAGACGACAT"
#seq2 = "TTTAGCATGCG"

# Test Case 6
#seq1 = "ATAGAC"
#seq2 = "TTTAGC"

match = 2
mismatch = -1
gap = -2
maxScore = 0
maxPosition = (0, 0)
rowcounter = 0
colcounter = 0
rows = len(seq1) 
cols = len(seq2)
max_d_count = 0
score_matrix = []
score_index_list = []
temp_list = []
threads = []


def score_function(row_index, col_index):
    global score_matrix
    global maxScore
    global maxPosition
    # print('\nWorker: %s' % num)
    similarity = match if seq1[row_index - 1] == seq2[col_index - 1] else mismatch
    diag_score = score_matrix[row_index - 1][col_index - 1] + similarity
    up_score = score_matrix[row_index - 1][col_index] + gap
    left_score = score_matrix[row_index][col_index - 1] + gap
    curMax = max(0, diag_score, up_score, left_score)
    if curMax > maxScore:
        maxScore = curMax
        maxPosition = (row_index, col_index)
    score_matrix[row_index][col_index] = curMax
    return

def createScoreMatrix(rows, cols):
    score_matrix = [[0 for col in range(cols+1)]for row in range(rows+1)]
    return score_matrix

def antidiagonals_list_generator(L):

    h, w = len(L), len(L[0])
    return [[L[p - q][q]
             for q in range(max(p - h + 1, 0), min(p + 1, w))]
            for p in range(h + w - 1)]

def antidiagonals_indices(rows, cols):

    for i in range(cols):
        p = 1
        q = i

        while(p < rows + 1 and q >= 1):
            #print(p, q)
            temp_list.append([p, q])
            p = p + 1
            q = q - 1

    for i in range(rows+1):
        p = i + 1
        q = cols
        while(p < rows+1 and q >= 1):
            #print(p, q)
            temp_list.append([p, q])
            p = p + 1
            q = q - 1

    # print(temp_list)


def traceback(score_matrix, start_pos):

    END, DIAG, UP, LEFT = range(4)
    aligned_seq1 = []
    aligned_seq2 = []
    x, y = start_pos
    move = nextMove(score_matrix, x, y)
    try:
        while move != END:
            if move == DIAG:
                aligned_seq1.append(seq1[x - 1])
                aligned_seq2.append(seq2[y - 1])
                x -= 1
                y -= 1
            elif move == UP:
                aligned_seq1.append(seq1[x - 1])
                aligned_seq2.append('-')
                x -= 1
            elif move == LEFT:
                aligned_seq1.append('-')
                aligned_seq2.append(seq2[y - 1])
                y -= 1
            else:
                move = END
            move = nextMove(score_matrix, x, y)
    except:
        move = END
    
    return ''.join(reversed(aligned_seq1)), ''.join(reversed(aligned_seq2))


def nextMove(score_matrix, x, y):

    # Assign the diagonal score
    diag = score_matrix[x - 1][y - 1]
    # print(diag)
    curr = score_matrix[x][y]
    # Assign insertion/deletion scores
    up = score_matrix[x - 1][y]
    left = score_matrix[x][y - 1]
    if (curr == diag + match or curr == diag + mismatch) and curr != up + gap and curr != left + gap : return 1
    elif up > left and up>=diag: return 2
    elif left > up and left>=diag: return 3
    elif curr == 0 : return 0
###############################################
## Creates the alignment string for printing ##
###############################################


def alignment_string(aligned_seq1, aligned_seq2):

    # Sets initial values
    idents, gaps, mismatches = 0, 0, 0
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


execution_start = time.time()
score_matrix = createScoreMatrix(rows, cols)
antidiagonals_indices(rows, cols)

antidiagonals = antidiagonals_list_generator(score_matrix)
# print(antidiagonals)
num_diag = len(antidiagonals)
offset_counter = 0
thread_start = time.time()
for i in range(num_diag - 2 ):         
    #print(i + 1)
    q = len(antidiagonals[i]) 
    x = 0
    threads.clear()
    # print(q)
    while(x < q):
        #print(temp_list[x + offset_counter])
        r_index = temp_list[x + offset_counter][0]
        c_index = temp_list[x + offset_counter][1]
        #print(r_index, c_index)
        x = x + 1

        t = threading.Thread(target=score_function, args=(r_index, c_index))
        threads.append(t)

    for x in range(len(threads)):
        threads[x].start()

    offset_counter = offset_counter + q - 1 

    # print(x)
thread_end = time.time()
print(seq1)
print(seq2)
for i in score_matrix:
    print(i)
# for y in range(len(score_matrix)):
# print(score_matrix[y])
'''
for i in range(5):
    t = threading.Thread(target=score_function, args=(i,))
    threads.append(t)
    t.start()
'''
print(maxPosition)
print(maxScore)
seq1_aligned, seq2_aligned = traceback(score_matrix, maxPosition)
assert len(seq1_aligned) == len(seq2_aligned)

execution_end = time.time()
print("Overall time to execute: {} seconds".format(
    execution_end - execution_start))
print("Time to run threads: {} seconds".format(thread_end - thread_start))
# Pretty print the results. The printing follows the format of BLAST results
# as closely as possible.
alignment_str, idents, gaps, mismatches = alignment_string(
    seq1_aligned, seq2_aligned)
alength = len(seq1_aligned)
print('\n')
print('Identities = {0}/{1} ({2:.1%}), Gaps = {3}/{4} ({5:.1%})'.format(idents,alength, idents / alength, gaps, alength, gaps / alength))
print('\n')
for i in range(0, alength, 60):
    seq1_slice = seq1_aligned[i:i + 60]
    print('Query  {0:<4}  {1}  {2:<4}'.format(
        i + 1, seq1_slice, i + len(seq1_slice)))
    print('             {0}'.format(alignment_str[i:i + 60]))
    seq2_slice = seq2_aligned[i:i + 60]
    print('Sbjct  {0:<4}  {1}  {2:<4}'.format(
        i + 1, seq2_slice, i + len(seq2_slice)))
    print('\n')
