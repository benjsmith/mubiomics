#!/usr/local/bin/python

#demultiplex.py
#Class and function definitions providing functionality in the mubiomics
#package.

#    Copyright (C) <2012>  <Benjamin C. Smith>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data import IUPACData
from patricia import *
from hamming import *
from numpy import *
import re, sys

class MultiIndexDict:
    """Thanks to Brad Chapman for posting this answer on Stack Overflow
    
    @usage: indata = SeqIO.index("f001", "fasta")
            pairdata = SeqIO.index("f002", "fasta")
            combo = MultiIndexDict(indata, pairdata)
            
            print combo['gi|3318709|pdb|1A91|'].description
            print combo['gi|1348917|gb|G26685|G26685'].description
            print combo["key_failure"]
    """
    def __init__(self, *indexes):
        self._indexes = indexes
    def __getitem__(self, key):
        for idx in self._indexes:
            try:
                return idx[key]
            except KeyError:
                pass
        raise KeyError("{0} not found".format(key))
    def __len__(self):
        length=0
        for idx in self._indexes:
            length=len(idx)+length
        return length

        
    
def quality_control(entry, av_score=40, min_score=10, \
    n_trim=0, min_len=80, win=10):
    """Generator function to trim FASTQ files according to quality criteria.

    Reads in a SeqRecord object from a FASTQ format file, if its average
    quality score is greater than av_score it scans through the nucleotides
    and returns all nucleotides until a window of size win has an average
    quality score of less than min_score. It then trims off the first n_trim
    nucleotides (which may be unwanted primer padding) and returns a trimmed 
    SeqRecord object if the length of the remaining sequence is greater than 
    min_length.

    >>>for record in quality_control(record_iterator, av_score=25, \
       min_score=15, n_trim=3, min_len=100, win=50]):
           print record.format("fasta").rstrip("\n")
    """
    for i, rec in enumerate(entry) :
        #get the list of basecall scores
        qual_list = rec[n_trim:].letter_annotations["phred_quality"]
        #scan through the list so that the end of the final window corresponds
        #to the last nucleotide in the read
        for j in arange(len(qual_list) - win):
            #test whether average score of window falls below minimum allowed
            #score
            if mean(qual_list[j:j+win]) < min_score :
                break
        #Test if window made it to end of list (has to be -1 because counting
        #in the for loop that sets j starts at zero)
        if j != len(qual_list)-win-1:
            # if not, return only the sequence before the failed window
            #providing it's longer than the min. length and has an average
            #score greater than the min. average score.
            if len(rec[n_trim:j]) >= min_len and \
            mean(rec[n_trim:j].letter_annotations["phred_quality"]) >= \
            av_score :
                yield rec[n_trim:j]
            else: 
                yield None
        else : #if the window made it to the end of the list without failing
            # return the entire sequence
            if len(rec[n_trim:]) >= min_len and \
            mean(rec[n_trim:].letter_annotations["phred_quality"]) >= \
            av_score :
                yield rec[n_trim:]
            else:
                yield None

def mmDNAtrie(seqs, mismatches):
    """Creates a dictionaried patricia tree of seqs with mismatches.
    
    Requires the patricia class. 
    
    seqs should be a list of strings, mismatches should be an integer of the
    number of acceptable mismatches. Returns a patricia trie in the form of
    a dictionary. The full tree can be viewed with triename._d, where triename
    is the name you assigned to the DNAtrie function call.
    
    This function was designed as a way to 
    quickly allow searches with mismatches for short (<12 bp) DNA barcodes.

    Gets slow i.e. takes a few seconds to build the trie, when doing 2 
    mismatches on barcodes of >15 bp. Not sure about memory requirements. 
    """
    trie = patricia()

    all_words = []
    words = seqs
    trie_buffer = []
    count = 0
    while count <= mismatches:
#        print count
#        print len(words)
        for word in words:
            if not trie.isWord(word):
                trie.addWord(word)
                all_words.append(word)
            for i, nt in enumerate(word):
                for alph in ['A', 'C', 'G', 'T']:
                    new_word = ''.join([word[0:i],alph,word[i+1:]])
                    if new_word not in trie_buffer:
                        trie_buffer.append(new_word)
                    
        words = trie_buffer
        trie_buffer=[]
        count += 1
#    for case in all_words:
#        print case
    return trie


def ptrie_search(ptrie, rec, q_len=8, st_pos=-1, max_st_pos=0):
    "Searches a patricia trie for a window, starting at st_pos, in rec."
    
    s_len = len(rec)
    if st_pos + 1 + q_len >= s_len:
        return None 
    pos = st_pos    
    m = False
    while m is False and pos+1 <= max_st_pos:
        pos+=1
        try:#ensure that if end of sequence is reached, loop exits properly
            s = rec[pos:pos+q_len]
        except (IndexError, UnboundLocalError):
            return None
        m = ptrie.isWord(s)
    if m: #Test to be certain loop exited with m=True
        return pos, s#if so, return position and located query
    else:
        return None

def ambiguous_search(rec, query, max_dist=0):
    """Search rec for query and return position of closest match.
    
    query can contain ambiguous nucleotides. Maximum acceptable edit distance
    between the query and the best match can be set with max_dist.
    @returns: (start_position, end_position)
    """
    win = len(query)
    dist_dict = {}
    for j in arange(len(rec) - win):#scan through sequence 
        test = rec[j:j+win]#select window to examine
        dist = ambiguous_seq_dist(test, query)#calculate edit distance
        if dist <= max_dist:
            dist_dict[j] = dist #add to a dictionary with (key, value) = (start_position, distance)
    #sort dist_dict according to value with lowest first and select first entry
    try:
        best_pos = sorted(dist_dict.iteritems(), key=lambda (k,v): (v,k))[0][0]
        return (best_pos, best_pos+win)
    except IndexError:
        return None

def ambiguous_seq_dist(t_seq, q_seq): 
    """Calculate Levenshtein distance between DNA sequences,
    t_seq and q_seq. 
        
    Can use ambiguous values in q_seq
    (like N = A or T or C or G, R = A or G etc.)  
    """ 
    dist = 0 
    for t_nt, q_nt in zip(t_seq, q_seq): 
        q_value = IUPACData.ambiguous_dna_values[q_nt] 
        if len(q_value) == 1: 
            if t_nt!=q_value:
                dist += 1
        else: 
            pattern = '[%s]' % q_value
            if not re.match(pattern, t_nt):
                dist += 1
    return dist
        
def identify_read(ptrie, rec, barcodes, ids, primers, bc_len=8, rpad=4, \
    max_pos=0, max_p_mismatch=2, bc_type='hamming_8') :
    """Finds the position of the barcode.
    
    Looks for the barcode in the mismatch DNA radix trie. If it finds it,
    decodes it using the Hamming decoding function then calculates the 
    Hamming distance between the input primer and the expected
    primer (accepting ambiguous DNA codes in the expected primer). It also
    simultaneously pulls out the id corresponding to the decoded barcode. After
    finding all possible barcodes in the input read, it takes the one with
    the smallest number of mismatches in the primer and returns the start 
    position of the read after trimming the barcode and primer, the id, the
    barcode from the read, it's decoded or true counterpart and the number of
    mismatches in the primer.
    """

    # stores data associated with possible barcode matches at varying window
    #  starts, so as to find which start position minimizes edit distance
    dist_dict = {}
    id_dict = {}
    old_bc_dict = {}
    new_bc_dict = {}
    seq_start_dict = {}
    pos = -1
    while len(str(rec.seq)) > bc_len :
        bc_match = ptrie_search(ptrie, str(rec.seq), q_len=bc_len, st_pos=pos \
                , max_st_pos=max_pos)
        if bc_match:
            pos = bc_match[0]
            if bc_type == 'hamming_8':
                true_bc = decode_barcode_8(bc_match[1]) # hamming, returns correct barcode.
                #if it fails to decode, then the try statement below fails
#             print "Decoded to " + true_bc
            else:
                #if not a hamming barcode, just calculate the distance between the
                #potential bc at this position against all barcodes in the list
                #of known barcodes and return the one with the smallest distance.
                bc_dist_dict = {}
                for bc in barcodes:
                    bc_dist_dict[bc] = ambiguous_seq_dist(bc, bc_match[1])
                bc_dist_list = sorted(bc_dist_dict.iteritems(), \
                                 key=lambda (k,v): (v,k))
                true_bc = bc_dist_list[0][0]
            try:
                i = barcodes[true_bc]
                new_id = ids[i]
                #if there is known primer sequence, minimize the edit distance
                #between it and the sequence where the primer should be, relative
                #to where the potential barcode is.
                if len(primers) != 0:
                    prim = primers[i]
                    len_prim = len(prim)
                    prim_start = pos + bc_len + rpad
                    prim_end = prim_start + len_prim
                    dist = \
                        ambiguous_seq_dist(str(rec.seq)[prim_start: \
                        prim_end], prim)
                    seq_start_dict[pos] = prim_end
                #if no primers, i.e. fragmented DNA, minimize the edit distance
                #of the potential barcodes.
                else:
                    len_prim = 0
                    dist = bc_dist_list[0][1]
                    seq_start_dict[pos] = pos + bc_len
                dist_dict[pos] = dist
                id_dict[pos] = new_id
                old_bc_dict[pos] = bc_match[1]
                new_bc_dict[pos] = true_bc
            except KeyError:
#                 print "Barcode not in list of barcodes."
                pass               
        else:
#             print "No barcode found at position."
            break
    try:
        best_pos = sorted(dist_dict.iteritems(), key=lambda (k,v): (v,k))[0][0]
        #print "======"
        #print "Min. primer dist = " + str(dist_dict[best_pos]) 
        if dist_dict[best_pos] <= max_p_mismatch :
            # return best starting position in sequence, best ID, old barcode,
            #  new barcode, edit distance, and length of primer
            return seq_start_dict[best_pos], id_dict[best_pos], \
               old_bc_dict[best_pos], new_bc_dict[best_pos], \
               dist_dict[best_pos], len_prim
        else:
            return 0, "Unassigned", "Unknown", "Unknown", dist_dict[best_pos] \
            , 0
    except IndexError:
        return 0, "Unassigned", "Unknown", "Unknown", None, 0

