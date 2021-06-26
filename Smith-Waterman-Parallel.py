import threading
import time
import numpy as np
from matplotlib import pyplot as plt

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
# seq1 = "ATAGACGACATGGGGACAGCAT"
# seq2 = "TTTAGCATGCGCATATCAGCAA"

# Test Case 5
#seq1 = "ATAGACGACAT"
#seq2 = "TTTAGCATGCG"

# Test Case 6
seq1 = "ATAGAC"
seq2 = "TTTAGC"

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
    aligned_seq_og1 = []
    aligned_seq_og2 = []
    x, y = start_pos
    move = nextMove(score_matrix, x, y)
    try:
        while move != END:
            if move == DIAG:
                aligned_seq1.append(seq1[x - 1])
                aligned_seq2.append(seq2[y - 1])
                aligned_seq_og1.append(seq1[x - 1])
                aligned_seq_og2.append(seq2[y - 1])
                x -= 1
                y -= 1
            elif move == UP:
                aligned_seq1.append(seq1[x - 1])
                aligned_seq2.append('-')
                aligned_seq_og1.append(seq1[x - 1])
                x -= 1
            elif move == LEFT:
                aligned_seq1.append('-')
                aligned_seq2.append(seq2[y - 1])
                aligned_seq_og2.append(seq2[y - 1])
                y -= 1
            else:
                move = END
            move = nextMove(score_matrix, x, y)
    except:
        move = END

    return ''.join(reversed(aligned_seq1)), ''.join(reversed(aligned_seq2)), ''.join(reversed(aligned_seq_og1)), ''.join(reversed(aligned_seq_og2))

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

    idents, gaps, mismatches = 0, 0, 0
    alignment_string = []

    for base1, base2 in zip(aligned_seq1, aligned_seq2):
        if base1 == base2:
            alignment_string.append('| ')
            idents += 1
        elif '-' in (base1, base2):
            alignment_string.append(' ')
            gaps += 1
        else:
            alignment_string.append(': ')
            mismatches += 1

    return ''.join(alignment_string), idents, gaps, mismatches

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

print(seq1)
print(seq2)

for i in score_matrix:
    print(i)

print(maxPosition)
print(maxScore)
seq1_aligned, seq2_aligned, seq1_aligned_og, seq2_aligned_og = traceback(score_matrix, maxPosition)
print(seq1_aligned_og)
print(seq2_aligned_og)

assert len(seq1_aligned) == len(seq2_aligned)

execution_end = time.time()

print("Overall time to execute: {} seconds".format(
    execution_end - execution_start))
print("Time to run threads: {} seconds".format(thread_end - thread_start))

alignment_str, idents, gaps, mismatches = alignment_string(
    seq1_aligned, seq2_aligned)

alength = len(seq1_aligned)

seq1_slice = seq1_aligned_og
seq2_slice = seq2_aligned_og

match_index_seq1 = seq1.index(seq1_slice)
match_index_seq2 = seq2.index(seq2_slice)

start1 = match_index_seq1
end1 = match_index_seq1 + len(seq1_slice)

start2 = match_index_seq2 
end2 = match_index_seq2 + len(seq1_slice)

result_matrix = score_matrix[start1+1:end1+1]
result_matrix = [row[start2+1:end2+1] for row in result_matrix]

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
            1,
            seq1_slice,
            alignment_str,
            idents*match + mismatches*mismatch + gaps*gap,
            2,
            seq2_aligned
        )
    )

ax2.set_title(f"Overall time to execute: {'%1.4f' % (execution_end - execution_start)} seconds        \nTime to run threads: {'%1.4f' % (thread_end - thread_start)} seconds        ")


for i in range(len(seq2[start2:end2])):
    for j in range(len(seq1[start1:end1])):
        text = ax1.text(j, i, result_matrix[i][j], ha="center", va="center", color="w")

pie_labels = [f'Matches({idents}/{alength})', f'Gaps({gaps}/{alength})', f'Missmatches({mismatches}/{alength})']

percents = [idents/alength, gaps/alength, mismatches/alength]

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