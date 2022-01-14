import math
#since the old version has some bugs , this is a new version based on some changes by Wang june 2020

def compute_similarity_score4alignment(matrix1, matrix2,return_index=False):
    ''' this script is for computing similairy score and position of max score when allign two matrix ,
        where matrix1 is matched one (e.g., pwms) and matrix2 is the known/target/initial one (e.g., rep) and be fixed!!
         Most part of this function is similar to computer_similarity_score but some minor changes for position purpose-
         Input matrix1 is pwm/predicted matrix, matrix2 is known/initial/rep matrix,
         Output  matched max_score, strand, position in matrix1 (small_index), matrix2 (large_index) Wang'''
    #add and modified by wang June 2020
    #here we fix small_index for matched/clustered ones
    #large_index for known/target/initial one
    small_index=(0,0)
    large_index=(0,0)
    temp_max_strand=''

    #removed by wang
    #assign short matrix to small
    if matrix1.shape < matrix2.shape:
        #add wang
        is_first_smaller=True
        small = matrix1
        large = matrix2
    else:
        #add wang
        is_first_smaller=False
        small = matrix2
        large = matrix1

    small_length = small.shape[0]
    large_length1 = large.shape[0]
    initial_pos = 0
    temp_max_score = 0
    #add wang use target matrix length
    minimum_allign_length = small_length + 1
    #add wang
    if minimum_allign_length > 7 and minimum_allign_length < 10:
        minimum_allign_length = int(minimum_allign_length * 0.95)
    if minimum_allign_length > 9 and minimum_allign_length < 14:
        minimum_allign_length = int(minimum_allign_length * 0.88)
    if minimum_allign_length > 13 and minimum_allign_length < 18:
        minimum_allign_length = int(minimum_allign_length * 0.80)
    if minimum_allign_length > 17:
        minimum_allign_length = int(minimum_allign_length * 0.75)

    minimum_allign_length -= 1

    #add wang
    #forward strand scan move small matrix against large one from left to righ
    #initial position for small matrix is always 0
    while initial_pos <= large_length1 - minimum_allign_length:
        temp_score = 0
        temp_len = 0
        for row in range(0, small_length):
            #add wang
            if row + initial_pos >= large_length1 :
                break
            temp_score += ((1 / math.sqrt(2)) * math.sqrt(
                math.pow(small[row, 0] - large[row + initial_pos, 0], 2) +
                math.pow(small[row, 1] - large[row + initial_pos, 1], 2) +
                math.pow(small[row, 2] - large[row + initial_pos, 2], 2) +
                math.pow(small[row, 3] - large[row + initial_pos, 3], 2)
            ))
            temp_len += 1
        temp_score = 1 - ((1 / temp_len) * temp_score)
        if temp_score > temp_max_score:
            #add wang
            #record aligned position (0, aligned_index_small) for small_matrix, (aligned_index_large, length) for large_matrix
            #if *_index[0]> *_index[1] then it is forward match
            small_index=(0, row if row + initial_pos < large_length1 else row -1)
            large_index=(initial_pos, min(row + initial_pos, large_length1 -1))
            temp_max_strand= 'F'
            temp_max_score = temp_score
            temp_max_direction='right'
        initial_pos += 1

    # forward strand to left side
    #add wang
    #forward strand scan , move small matrix against large one from right to left
    #initial position for small matrix is temp_max_pos
    initial_pos = 0
    while initial_pos <= small_length - minimum_allign_length:
        temp_score = 0
        temp_len = 0
        for row in range(0, small_length):
            if row + initial_pos >= small_length:
                break
            temp_score += ((1 / math.sqrt(2)) * math.sqrt(
                math.pow(large[row, 0] - small[row + initial_pos, 0], 2) +
                math.pow(large[row, 1] - small[row + initial_pos, 1], 2) +
                math.pow(large[row, 2] - small[row + initial_pos, 2], 2) +
                math.pow(large[row, 3] - small[row + initial_pos, 3], 2)
            ))
            temp_len += 1
        temp_score = 1 - ((1 / temp_len) * temp_score)
        if temp_score > temp_max_score:
            #add wang
            #record (aligned_index_small, length) for small_matrix, (0, aligned_index_large) for large matrix
            #small_index[0]<small_index[1] is forward match
            #large_index[0]<large_index[1] is forward match
            small_index = (initial_pos, min(row + initial_pos, small_length-1))
            large_index = (0, row if row + initial_pos < small_length else row -1)
            temp_max_strand='F'
            temp_max_score = temp_score
            temp_max_direction='left'
        initial_pos += 1


    initial_pos = 0
    reverse_pos = list(range(small_length - 1, -1, -1))
    #add wang
    #rever strand scan, move reversed samll matrix from left to right against large one
    #initial position to small matrix is always reverse_pos[0]
    while initial_pos <= large_length1 - minimum_allign_length:
        temp_score = 0
        temp_len = 0
        for row in range(0, small_length):
            if row + initial_pos >= large_length1:
                break
            temp_score += ((1 / math.sqrt(2)) * math.sqrt(
                math.pow(small[reverse_pos[row], 3] - large[row + initial_pos, 0], 2) +
                math.pow(small[reverse_pos[row], 2] - large[row + initial_pos, 1], 2) +
                math.pow(small[reverse_pos[row], 1] - large[row + initial_pos, 2], 2) +
                math.pow(small[reverse_pos[row], 0] - large[row + initial_pos, 3], 2)
            ))
            temp_len += 1
        temp_score = 1 - ((1 / temp_len) * temp_score)
        if temp_score > temp_max_score:
            #add wang
            #record (reverse_pos[0], reverse_pos[alighed_index_small]) for small_reversed_matrix,
            #small_index[0]> small_index_[1] for reversed small matrix
            #(aligned_index_large, length) for large_matrix
            #large_index[0]< large_index[1] forward large matrix
            small_index = (reverse_pos[0], reverse_pos[row if row + initial_pos < large_length1 else row -1])
            large_index = (initial_pos, min(large_length1-1,row + initial_pos ))
            temp_max_score = temp_score
            temp_max_strand = 'R'
            temp_max_direction='right'

        initial_pos += 1

    #add wang
    #recerse strand alighment, move reversed small matrix from left to right against large one
    #the initial position of small matrix
    initial_pos = 0
    while initial_pos <= small_length - minimum_allign_length:
        temp_score = 0
        temp_len = 0
        for row in range(0, small_length):
            #add wang
            if row + initial_pos >= small_length :
                break
            temp_score += ((1 / math.sqrt(2)) * math.sqrt(
                math.pow(small[reverse_pos[row + initial_pos], 3] - large[row, 0], 2) +
                math.pow(small[reverse_pos[row + initial_pos], 2] - large[row, 1], 2) +
                math.pow(small[reverse_pos[row + initial_pos], 1] - large[row, 2], 2) +
                math.pow(small[reverse_pos[row + initial_pos], 0] - large[row, 3], 2)
            ))
            temp_len += 1
        temp_score = 1 - ((1 / temp_len) * temp_score)
        if temp_score > temp_max_score:
            #add wang
            if row + initial_pos >= small_length:
                end = reverse_pos[small_length -1]
            else:
                end = reverse_pos[row + initial_pos]
            #add wang
            #record (reverse_pos[aligned_index_small], length) for small_reversed_matrix
            #(0, aligned_index_large) for large_matrix
            #small_index[0]>small_index[1] for reverse small matrix
            #large_index[0]<large_index[1] for forward large matrix
            #This one may be needed double check latter ?
            small_index= (reverse_pos[initial_pos], end)
            large_index= (0, row if row + initial_pos < small_length else row -1 )
            temp_max_score = temp_score
            temp_max_strand = 'R'
            temp_max_direction='left'

        initial_pos += 1


    if return_index is False:
        return temp_max_score
    else:
        if is_first_smaller:
            return temp_max_score, temp_max_strand, small_index, large_index
        else:
            return temp_max_score, temp_max_strand, large_index, small_index

