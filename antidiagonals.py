def antidiagonals_list_generator(score_matrix):
    h, w = len(score_matrix), len(score_matrix[0])

    return [[score_matrix[p - q][q]
             for q in range(max(p - h + 1, 0), min(p + 1, w))]
            for p in range(h + w - 1)]

def antidiagonals_indices(rows, cols):
    temp_list = []
    
    for i in range(cols):
        p = 0
        q = i

        while(p < rows and q >= 0):
            #print(p, q)
            temp_list.append([p, q])
            p = p + 1
            q = q - 1

    for i in range(rows):
        p = i + 1
        q = cols - 1
        while(p < rows and q >= 1):
            #print(p, q)
            temp_list.append([p, q])
            p = p + 1
            q = q - 1

    # print(temp_list)
    return temp_list;