import math

def compute_similarity_score4alignment(matrix1, matrix2,return_index=False):
    ''' this script is for computing similairy score and position of max score when allign two matrix ,
        where matrix1 is matched one (e.g., pwms) and matrix2 is the known/target/initial one (e.g., rep) and be fixed!!
         Most part of this function is similar to computer_similarity_score but some minor changes for position purpose-
         Input matrix1 is pwm/predicted matrix, matrix2 is known/initial/rep matrix,
         Output  matched max_score, strand, position in matrix1 (small_index), matrix2 (large_index) Wang'''

    small_index=(0,0)
    large_index=(0,0)
    temp_max_strand=''

    if matrix1.shape < matrix2.shape:
        is_first_smaller=True
        small = matrix1
        large = matrix2
    else:
        is_first_smaller=False
        small = matrix2
        large = matrix1

    small_length = small.shape[0]
    large_length1 = large.shape[0]
    initial_pos = 0
    temp_max_score = 0
    minimum_allign_length = small_length + 1
    if minimum_allign_length > 7 and minimum_allign_length < 10:
        minimum_allign_length = int(minimum_allign_length * 0.95)
    if minimum_allign_length > 9 and minimum_allign_length < 14:
        minimum_allign_length = int(minimum_allign_length * 0.88)
    if minimum_allign_length > 13 and minimum_allign_length < 18:
        minimum_allign_length = int(minimum_allign_length * 0.80)
    if minimum_allign_length > 17:
        minimum_allign_length = int(minimum_allign_length * 0.75)

    minimum_allign_length -= 1


    while initial_pos <= large_length1 - minimum_allign_length:
        temp_score = 0
        temp_len = 0
        for row in range(0, small_length):
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
            small_index=(0, row if row + initial_pos < large_length1 else row -1)
            large_index=(initial_pos, min(row + initial_pos, large_length1 -1))
            temp_max_strand= 'F'
            temp_max_score = temp_score
            temp_max_direction='right'
        initial_pos += 1

    # forward strand to left side

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

            small_index = (initial_pos, min(row + initial_pos, small_length-1))
            large_index = (0, row if row + initial_pos < small_length else row -1)
            temp_max_strand='F'
            temp_max_score = temp_score
            temp_max_direction='left'
        initial_pos += 1


    initial_pos = 0
    reverse_pos = list(range(small_length - 1, -1, -1))

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

            small_index = (reverse_pos[0], reverse_pos[row if row + initial_pos < large_length1 else row -1])
            large_index = (initial_pos, min(large_length1-1,row + initial_pos ))
            temp_max_score = temp_score
            temp_max_strand = 'R'
            temp_max_direction='right'

        initial_pos += 1


    initial_pos = 0
    while initial_pos <= small_length - minimum_allign_length:
        temp_score = 0
        temp_len = 0
        for row in range(0, small_length):
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
            if row + initial_pos >= small_length:
                end = reverse_pos[small_length -1]
            else:
                end = reverse_pos[row + initial_pos]

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

