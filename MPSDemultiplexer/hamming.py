#!/usr/local/bin/python

##################################################################
#Created using code from Hamady et al. (2008) supplemental.#
#################################################################

from numpy import *

# current encoding scheme
INT_TO_BS = {0:"00", 1:"01", 2:"10", 3:"11"}
CUR_ENC_FO = {'A': 3, 'C': 2, 'T': 0, 'G': 1}
CUR_REV_ENC_SI = { "11":"A", "10":"C", "00":"T", "01":"G"}
    

def calc_parity_vector(parity_vector):
    """ Returns even or odd parit for parity vector """
    return reduce(lambda x, y: x^y, parity_vector[1:])

def calc_syndrome(codeword, n):
    """ Calculate syndrome and correct codeword if possible """
    sym = 0
    for i in range(1,n):
        if codeword[i]:
            sym ^= i
    extra_parity = calc_parity_vector(codeword)
    if extra_parity == codeword[0]:
        if sym == 0:
            return 0, sym
        else:
            return 2, sym
    else:
        if sym >= n:
            pass
        else:
            codeword[sym] ^= 1
        return 1, sym
    
    
def nt_to_cw(cur_enc, cur_nt):
    """ Convert nt sequence to codeword """
    return array(map(int, ''.join([INT_TO_BS[cur_enc[x]] for x in \
    cur_nt])))
    
    
def unpack_bitstr(rev_cur_bit, bitstr):
    """ Unpack bistring into nt sequence """
    bstr_len = len(bitstr)
    return ''.join([rev_cur_bit[bitstr[i:i+2]] for i in range(0, bstr_len, \
    2)])
    
def decode_barcode_8(nt_barcode):
    """ Decode length 8 barcode (16 bits) """
    # check proper length
    if len(nt_barcode) != 8:
        raise ValueError, "barcode must be 8 nt long."
    
    # check valid characters
    if set(list(nt_barcode)).difference(CUR_ENC_FO.keys()):
        raise ValueError, "Only A,T,C,G valid chars."
    
    # decode
    decoded = nt_to_cw(CUR_ENC_FO, nt_barcode)
    num_errors, sym = calc_syndrome(decoded, 16)
    
    # check errors
    if num_errors > 1:
#        raise ValueError, "2 bp error detected."
        pass
    
    # convert corrected codeword back to nt sequence
    if num_errors == 1:
        nt_barcode = unpack_bitstr(CUR_REV_ENC_SI, ''.join(map(str, \
    decoded)))
    
    return nt_barcode