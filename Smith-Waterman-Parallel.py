import threading
import time
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

#import timeit
# from array import array

# Test Case 1
# seq1 = "ATAGACGACATGGGGACAGCATATAGACGACATGGGGACAGCATATAGACGACATGGGGACAGCATATAGACGACATGGGGACAGCATATAGACGACATGGGGACAGCATATAGACGACATGGGGACAGCATATAGACGACATGGGGACAGCATATAGACGACATGGGGACAGCAT"
# seq2 = "TTTAGCATGCGCATATCAGCAATTTAGCATGCGCATATCAGCAATTTAGCATGCGCATATCAGCAATTTAGCATGCGCATATCAGCAATTTAGCATGCGCATATCAGCAATTTAGCATGCGCATATCAGCAATTTAGCATGCGCATATCAGCAATTTAGCATGCGCATATCAGCAA"

# Test Case 2
# seq1 = "ATAGACGACATGGGGACAGCATATAGACGACATGGGGACAGCATATAGACGACATGGGGACAGCATATAGACGACATGGGGACAGCAT"
# seq2 = "TTTAGCATGCGCATATCAGCAATTTAGCATGCGCATATCAGCAATTTAGCATGCGCATATCAGCAATTTAGCATGCGCATATCAGCAA"

# Test Case 3
# seq1 = "ATAGACGACATGGGGACAGCATATAGACGACATGGGGACAGCAT"
# seq2 = "TTTAGCATGCGCATATCAGCAATTTAGCATGCGCATATCAGCAA"

# Test Case 4
seq1 = "ATAGACGACATGGGGACAGCATCCC"
seq2 = "TTTAGCATGCGCATATCAGCAAAAA"

# Test Case 5
# seq1 = "ATAGACGACAT"
# seq2 = "TTTAGCATGCG"

# Test Case 6
# seq1 = "AAAAGCATTTTTTGCA"
# seq2 = "CCCCGCAGGGGGGGCA"

match = 1
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
            temp_list.append([p, q])
            p = p + 1
            q = q - 1

    for i in range(rows+1):
        p = i + 1
        q = cols
        while(p < rows+1 and q >= 1):
            temp_list.append([p, q])
            p = p + 1
            q = q - 1


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
                aligned_seq2.append('_')
                x -= 1
            elif move == LEFT:
                aligned_seq1.append('_')
                aligned_seq2.append(seq2[y - 1])
                y -= 1
            else:
                move = END
            move = nextMove(score_matrix, x, y)
    except:
        move = END
    
    return ''.join(reversed(aligned_seq1)), ''.join(reversed(aligned_seq2))


def nextMove(score_matrix, x, y):
    diag = score_matrix[x - 1][y - 1]
    curr = score_matrix[x][y]
    up = score_matrix[x - 1][y]
    left = score_matrix[x][y - 1]
    if (curr == diag + match or curr == diag + mismatch) and curr != up + gap and curr != left + gap : return 1
    elif up > left and up>=diag: return 2
    elif left > up and left>=diag: return 3
    elif curr == 0 : return 0


def alignment_string(aligned_seq1, aligned_seq2):
    match, gaps, mismatches = 0, 0, 0
    alignment_string = []

    for base1, base2 in zip(aligned_seq1, aligned_seq2):

        if base1 == base2:
            alignment_string.append('| ')
            match += 1

        elif '-' in (base1, base2):
            alignment_string.append(' ')
            gaps += 1

        else:
            alignment_string.append(': ')
            mismatches += 1

    return ''.join(alignment_string), match, gaps, mismatches


execution_start = time.time()
score_matrix = createScoreMatrix(rows, cols)
antidiagonals_indices(rows, cols)

antidiagonals = antidiagonals_list_generator(score_matrix)
num_diag = len(antidiagonals)
offset_counter = 0
thread_start = time.time()
for i in range(num_diag - 2 ):         
    q = len(antidiagonals[i]) 
    x = 0
    threads.clear()
    while(x < q):
        r_index = temp_list[x + offset_counter][0]
        c_index = temp_list[x + offset_counter][1]
        x = x + 1

        t = threading.Thread(target=score_function, args=(r_index, c_index))
        threads.append(t)

    for x in range(len(threads)):
        threads[x].start()

    offset_counter = offset_counter + q - 1 

thread_end = time.time()
seq1_aligned, seq2_aligned = traceback(score_matrix, maxPosition)
assert len(seq1_aligned) == len(seq2_aligned)

execution_end = time.time()
alignment_str, match, gaps, mismatches = alignment_string(
    seq1_aligned, seq2_aligned)
alength = len(seq1_aligned)

for i in range(0, alength, 60):
    seq1_slice = seq1_aligned[i:i + 60]
    seq2_slice = seq2_aligned[i:i + 60]

match_index_seq1 = seq1.index(seq1_slice)
match_index_seq2 = seq2.index(seq2_slice)

start1 = match_index_seq1 - 2 if match_index_seq1 >= 2 else 0
end1 = match_index_seq1 + len(seq1_slice) + 2 if match_index_seq1 + len(seq1_slice) + 2 < len(seq1) else match_index_seq1 + len(seq1_slice)

start2 = match_index_seq2 - 2 if match_index_seq2 >= 2 else 0
end2 = match_index_seq2 + len(seq1_slice) + 2 if match_index_seq2 + len(seq1_slice) + 2 < len(seq1) else match_index_seq2 + len(seq1_slice)

result_matrix = score_matrix[start1:end1]
result_matrix = [row[start1:end1] for row in result_matrix]

fig, (ax1, ax2) = plt.subplots(1,2)
result_matrix = np.flip(result_matrix, axis=0)
im = ax1.imshow(result_matrix)

ax1.set_xticks(np.arange(len(seq2[start2:end2])))
ax1.set_yticks(np.arange(len(seq1[start1:end1])))
ax1.set_xticklabels(seq2[start2:end2])
ax1.set_yticklabels(seq1[start1:end1][::-1])
ax1.set_xlabel(seq2[start2:end2])
ax1.set_ylabel(seq1[start1:end1])
ax1.set_title('Seq {0:<4}  {1}\n               {2}\nSeq {4:<4}  {5}\nScore: {3:<4}'.format(
            i + 1,
            seq1_slice,
            alignment_str[i:i + 60],
            i+len(seq1_slice),
            i + 2,
            seq2_slice
        )
    )

ax2.set_title(f"Overall time to execute: {'%1.4f' % (execution_end - execution_start)} seconds        \nTime to run threads: {'%1.4f' % (thread_end - thread_start)} seconds        ")

for i in range(len(seq1[start1:end1])):
    for j in range(len(seq2[start1:end1])):
        text = ax1.text(j, i, result_matrix[i][j], ha="center", va="center", color="w")

pie_labels = [f'Matches({match}/{alength})', f'Gaps({gaps}/{alength})', f'Missmatches({mismatches}/{alength})']

# ax2.axis('equal')

percents = [match/alength, gaps/alength, mismatches/alength]

def filter_zeroes(pct):
    return ('%1.1f' % pct) if pct > 0 else ''

def get_new_labels(sizes, labels):
    new_labels = [label if size > 0 else '' for size, label in zip(sizes, labels)]
    return new_labels

ax2.pie(percents, labels=get_new_labels(percents, pie_labels), autopct=filter_zeroes)

figManager = plt.get_current_fig_manager()
figManager.window.state('zoomed')
plt.subplots_adjust(wspace=1)
plt.show()