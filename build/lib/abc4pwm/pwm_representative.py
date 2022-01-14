import numpy as np
from abc4pwm.convert_count_to_pwm import motif_weight2p
import os
#add wang May 27 2020
#rewrite by Wang June 2020
from abc4pwm.similarity_score import compute_similarity_score4alignment

#Since the old version does not work.
#here is a new implementation of pwm  alignment by Wang June 2020


class PWMRepresentative(object):

    def __init__(self):
        pass

    @staticmethod
    def __calculate_representative_length_and_pwm(pwms,pwms_name):
        ''' Input pwm and pwm_name list,
            OutPut return both representative motif length and selected index of initial rep pwm
        added wang'''
        len_list = []
        for pwm in pwms:
            len_list.append(pwm.shape[0])
        allign_length = int(np.median(len_list))
        allign_length2= np.min([i for i in len_list if i >= allign_length ])
        sl_pwm_idx=len_list.index(allign_length2)
        return allign_length2, sl_pwm_idx

    def __make_reverse_strand(self, pwm):
        ''' make reverse complementary sequence but assume input matrix is 4 columns -wang'''
        tmp = pwm.copy()
        tmp[:, 0] = pwm[:, 3]
        tmp[:, 1] = pwm[:, 2]
        tmp[:, 2] = pwm[:, 1]
        tmp[:, 3] = pwm[:, 0]
        return tmp[::-1]

  #  def __make_reverse_strand4position(self, pwm):
  #      ''' make reverse complementary sequence but assume input matrix is 4 columns and no flip of row orders -wang'''
  #      tmp = pwm.copy()
  #      tmp[:, 0] = pwm[:, 3]
  #      tmp[:, 1] = pwm[:, 2]
  #      tmp[:, 2] = pwm[:, 1]
  #      tmp[:, 3] = pwm[:, 0]
  #      return tmp

    def __normalize_by_count(self, pwm, count):
        ''' mean of input matrix based on number counts- wang'''
        tmp = pwm.copy()
        tmp_count = count.copy()
        tmp_count[tmp_count == 0] = 1
        for i in range(tmp.shape[1]):
            tmp[:, i] = tmp[:, i] / tmp_count
        return tmp

    def __add_to_rep4position(self, rep, pwm, rep_size, count_full):
            '''add aligned pwm to representative motif based on olde perl code
             Input rep initial/known matrix, pwm predicted matrixm, but count_full is not used at here
              Output aligned formatch rep and reverse match rev_rep, and their corresponding count and rev_count - wang'''
            #rep_copy = np.zeros(rep.shape)
            rep_copy = np.zeros((rep_size,4))
            rev_rep_copy = np.zeros((rep_size,4))
            #here copy initial rep matrix which is fixed and does not need to normalize it based on count
            rep_norm = rep.copy()
            count = np.zeros((rep_size))
            rev_count = np.zeros((rep_size))
            #find matched position for pwm
            scores, strands, index_pwm, index_rep = self.__get_sub_index4position(pwm, rep_norm)

            #find matched index for pwm and rep
            #check whether is pwm forward or pwm reverse matches to rep initial forard or reverse matrix
            is_pwm_reverse=False
            if index_pwm[0]< index_pwm[1]:
                #forward matched
                pwm_st_idx=index_pwm[0]
                pwm_ed_idx=index_pwm[1]+1
            else:
                #reverse matched
                is_pwm_reverse=True
                pwm_st_idx= index_pwm[1]
                pwm_ed_idx= index_pwm[0] + 1

            is_rep_reverse =False
            if index_rep[0] < index_rep[1] :
                #forward matched
                rep_st_idx= index_rep[0]
                rep_ed_idx= index_rep[1] +1
            else :
                rep_st_idx= index_rep[0]
                rep_ed_idx= index_rep[1]
                rep_len=rep.shape[0]
                idx2rep = np.array(range(rep_len))
                reverse_idx2rep = idx2rep[::-1]
                is_rep_reverse=True

            if not is_pwm_reverse and not is_rep_reverse :
                #pwm forward matches rep forward
                rep_copy[rep_st_idx:rep_ed_idx ] = pwm[pwm_st_idx: pwm_ed_idx ]
                count[rep_st_idx:rep_ed_idx ] += 1
            elif is_pwm_reverse and not is_rep_reverse:
                #pwm reverse matches to rep forward
                tmp_rep_copy= np.zeros((pwm.shape[0] ,4))
                tmp_rep_copy[0:rep_ed_idx-rep_st_idx ] = pwm[pwm_st_idx:pwm_ed_idx ]
                tmp_rep_copy = self.__make_reverse_strand(tmp_rep_copy)
                rep_copy[rep_st_idx: rep_ed_idx]= tmp_rep_copy[0:rep_ed_idx-rep_st_idx]
                count[rep_st_idx:rep_ed_idx ] += 1
            elif not is_pwm_reverse and is_rep_reverse:
                #pwm forward matches to rep reverse and keep it in a separate matrix
                rev_rep_copy[reverse_idx2rep[rep_st_idx] : reverse_idx2rep[rep_ed_idx]+1] = pwm[pwm_st_idx: pwm_ed_idx]
                rev_count[reverse_idx2rep[rep_st_idx] :reverse_idx2rep[rep_ed_idx] +1 ] +=1
            else:
                print("Something wrong in pwm_repreesntative.py! Please check !")

            return rep_copy, count, rev_rep_copy, rev_count

    def __get_sub_index4position(self, pwm, rep):
        '''find alighned position index of PWMs -wang '''
        pwm_p = motif_weight2p(pwm)
        rep_p = motif_weight2p(rep)
        #convert matrix row is motif length , column is nucleotides
        pwm_p=pwm_p.transpose()
        rep_p=rep_p.transpose()
        score, strand, index_pwm , index_rep = compute_similarity_score4alignment(pwm_p, rep_p, True)
        return score, strand, index_pwm, index_rep

    def calculate_representative4position(self, pwms, initial_value=None,pwms_name=None):
        #changed wang based on old perl code
        if initial_value is None:
            #add wang
            rep_size, sl_pwm_idx= self.__calculate_representative_length_and_pwm(pwms,pwms_name)
            initial_value = pwms[sl_pwm_idx].copy()
            # print(' Select ' +pwms[sl_pwm_idx] + ' as initial PWM')
            #print(pwms[sl_pwm_idx].shape)
            #print(rep_size)
            #remove selected initial pwm from the pwm list because it will be added in the first round of alignment
            del pwms[sl_pwm_idx]
            del pwms_name[sl_pwm_idx]
        else:
            rep_size=initial_value.shape[0]
            rep_size2, sl_pwm_idx2 = self.__calculate_representative_length_and_pwm(pwms,pwms_name)
            if (rep_size2 - rep_size) > 2:
                rep_size = rep_size2
            if (rep_size - rep_size2) > 2:
                rep_size = rep_size2


        rep = np.zeros((rep_size, pwms[0].shape[1]))
        rev_rep = np.zeros((rep_size, pwms[0].shape[1]))
        sum_rep_size=rep_size
        sum_rep = np.zeros((sum_rep_size, pwms[0].shape[1]))
        rev_sum_rep = np.zeros((sum_rep_size, pwms[0].shape[1]))
        sum_count= np.zeros((sum_rep_size))
        rev_sum_count= np.zeros((sum_rep_size))

        count = np.zeros((rep_size))
        shape = initial_value.shape[0]
        diff = abs(shape - rep_size)
        if shape >= rep_size:
            rep = initial_value[int(diff / 2):rep_size + int(diff / 2)]
            count += 1
        else:
            #rep[int(diff / 2):shape + int(diff / 2)] = initial_value
            #count[int(diff / 2):shape + int(diff / 2)] += 1
            #force rep start from 0 position
            rep[0:shape ] = initial_value
            count[0:shape] +=1

        for idx, pwm in enumerate(pwms):
            #add wang
            # print(pwms_name[idx])
            r, c, rev_r, rev_c = self.__add_to_rep4position(rep, pwm, rep_size, count)
            #r is pwm matches forward rep
            #rev_r is pwm matches reverse rep
            #here may have to check wheter r and rev_r are zero or not before adding them ?? wang
            if np.max(r) >0:
              #forward intitial matrix matched
              if idx ==0 :
                #in the first round alignment, add intial matrix to aligned new matrix
                sum_rep = r + rep
                sum_count = 1 + c
              else:
                sum_rep += r
                sum_count += c

            if np.max(rev_c) >0:
              #reverse intitial matrix matched
              if idx ==0:
                #in the first round alignment, add initial matrix to aligned new matrix 
                rev_sum_rep = rev_r + self.__make_reverse_strand(rep)
                rev_sum_count = 1+ rev_c
              else:
                rev_sum_rep += rev_r
                rev_sum_count += rev_c
        #remove rows with all columns equals zero
        #count = sum_count[~np.all(sum_rep ==0, axis =1)]
        #sum_rep = sum_rep[~np.all(sum_rep == 0, axis=1)]
        rep = self.__normalize_by_count(sum_rep, sum_count)
        rev_rep = self.__normalize_by_count(rev_sum_rep,rev_sum_count)

        return rep, rev_rep, sum_count, rev_sum_count

    def calculate_representative(self, pwms, initial_value=None, pwms_name=None):
        rep, rev_rep, rep_count, rev_rep_count = self.calculate_representative4position(pwms, initial_value, pwms_name)
        # added wang
        if np.max(rep_count) > 0 and np.max(rev_rep_count) > 0:
            # if both forward and reversed matched exists then do the alignment again betwen these two matrix
            new_pwms = []
            new_pwms.append(rep)
            new_pwms.append(rev_rep)
            potential_rep2 = None
            pwm_name_list2 = ['rep', 'rev_rep']
            # align rep forward matched matrix to rep reverse matched matrix
            rep2, rev_rep2, rep_count2, rev_rep_count2 = self.calculate_representative4position(new_pwms, potential_rep2, pwm_name_list2)

            if np.max(rep_count2) > 1:
                # only consider results with count greater 1
                out_rep = rep2
                out_count = rep_count2
            else:
                out_rep = rev_rep2
                out_count = rev_rep_count2

        elif np.max(rep_count) > 0 and np.max(rev_rep_count) == 0:
            # only forward match to initial matrix exist
            out_rep = rep
            out_count = rep_count

        elif np.max(rev_rep_count) > 0 and np.max(rep_count) == 0:
            # only reverse match to inital matrix exit
            out_rep = rev_rep
            out_count = rev_rep_count

        else:
            print("Errors in alighment !")

        return out_rep, out_count
