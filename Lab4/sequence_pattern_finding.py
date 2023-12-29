# -*- coding: utf-8 -*-

####################################################################################
### implementing pattern finding in sequences
####################################################################################
def search_first_occ(seq, pattern):
    found = False
    i = 0
    while i <= len(seq)-len(pattern) and not found:
        j = 0
        while j < len(pattern) and pattern[j]==seq[i+j]: 
            j = j + 1
        if j== len(pattern): found = True
        else: i += 1
    if found: return i
    else: return -1
    
def search_all_occurrences(seq, pattern):
    ''' returns a list of the indices where the pattern occurrs within seq'''
    res = []
    i=0
    while i < len(seq):
        r = search_first_occ(seq[i:],pattern)
        if r == -1:
            return res
        else:
            res = res + [r]
            i = i+r+1
    return res

def search_all_occurrences_mismatch(seq, pattern, m):
    '''Find all occurrences of pattern in seq with maximum hamming distance of m'''
    res = []
    i=0
    l= len(pattern)
    while i < len(seq):
        if hamming_distance(pattern,seq[i:i+l]) < m:
            res = res + [i]
        i = i+ 1
    return res

def hamming_distance(s, p):
    dist = 0
    if len(s) != len(p):
        return -1
    else:
        for i in range(0, len(s)):
            if s[i] != p[i]:
                dist = dist + 1
    return dist




def test_pat_search():
    ''' test the implemented function search_all_occurrences here ''' 
    seq = input("Input sequence: ")
    pat = input("Input pattern: ")
    print(search_all_occurrences(seq, pat))
    print(search_all_occurrences_mismatch(seq, pat, 2))
    
    
#test_pat_search()    



### simpler implementation of Booyer-Moore
def simpleBooyerMoore(seq, pattern):
    """Very simplified version of Booyer-Moore with the BCR rule"""
    alphabet = "".join(set(list(seq)))
    # process Bad-character rule
    occ = {}
    for symb in alphabet:
        occ[symb] = -1
    for j in range(len(pattern)):
        c = pattern[j]
        occ[c] = j 
    # search pattern    
    res = []
    i = 0
    while i <= len(seq) - len(pattern):
        j = len(pattern) - 1
        while j >= 0 and pattern[j] == seq[j+i]:
            j -= 1
        if (j < 0):
            res.append(i)
            i += 1
        else:
            c = seq[j + i]
            i += max(1, j - occ[c]) 
    return res


#seqDNA = "ATAGAATAGATAATAGTC"
#print( search_first_occ(seqDNA, "GAAT") )
#print( search_first_occ(seqDNA, "TATA") )
#print( search_all_occurrences(seqDNA, "AAT") )


def repeated_subsequences_frequency(dna_seq, k = 10):
    '''Write a function that, given a DNA sequence, allows to detect if there are repeated sequences of size k
    The result should be a dictionary with sub-sequences as keys, and their frequency as values.'''
    dic ={}
    i = 0
    while i<len(dna_seq)-k:
        seq= dna_seq[i:i+k]
        if seq in dic.keys():
            i = i+1
        else:
            res = search_all_occurrences(dna_seq[i:],seq)
            if len(res) > 1:
                dic[seq] = len(res)
            i = i+1
    return dic


def test():
    seqDNA = "ATTATACACAATCFCATACATATT"
    print( simpleBooyerMoore(seqDNA, "CAAT") )
    import read_fasta
    seq = read_fasta.read_Fasta("test_files/HBA1.DNA.fasta")[0]
    print(seq)
    pattern="ACTCTT"
    print(search_all_occurrences(seq,pattern))
    print(search_all_occurrences_mismatch(seq,pattern,2))
    print(repeated_subsequences_frequency(seq,10))

test()

####################################################################################
### implementing pattern finding in sequences
####################################################################################

def validate_dna_re (seq):
    from re import search
    if search("[^ACTGactg]", seq) != None:
        return False
    else:
        return True
    

def translate_codon_re (cod):
    import re
    if re.search("GC.", cod): aa = "A"
    elif re.search("TG[TC]", cod): aa = "C" 
    elif re.search("GA[TC]", cod): aa = "D"
    elif re.search("GA[AG]", cod): aa = "E"
    elif re.search("TT[TC]", cod): aa = "F"
    elif re.search("GG.", cod): aa = "G"
    elif re.search("CA[TC]", cod): aa = "H"
    elif re.search("AT[TCA]", cod): aa = "I"
    elif re.search("AA[AG]", cod): aa = "K"
    elif re.search("TT[AG]|CT.", cod): aa = "L"
    elif re.search("ATG", cod): aa = "M"
    elif re.search("AA[TC]", cod): aa = "N"
    elif re.search("CC.", cod): aa = "P"
    elif re.search("CA[AG]", cod): aa = "Q"
    elif re.search("CG.|AG[AG]", cod): aa = "R"
    elif re.search("TC.|AG[TC]", cod): aa = "S"
    elif re.search("AC.", cod): aa = "T"
    elif re.search("GT.", cod): aa = "V"
    elif re.search("TGG", cod): aa = "W"
    elif re.search("TA[TC]", cod): aa = "Y"
    elif re.search("TA[AG]|TGA", cod): aa = "_";
    else: aa = ""     
    return aa

def find_pattern_re (seq, pat):
    from re import search
    mo = search(pat, seq)
    if (mo != None):
        return mo.span()[0]
    else:
        return -1

def find_all_occurrences_re (seq, pat):
    from re import finditer
    mos = finditer(pat, seq)
    res = []
    for x in mos:
        res.append(x.span()[0])
    return res

def find_all_overlap(seq, pat):
    return find_all_occurrences_re(seq, "(?="+pat+")")

def finds_patterns_frequency(seqA, seqB, low, high):
    i = 0
    dic = {}
    while i < len(seqA)-low:
        for j in range(low,high+1):
            if i+j < len(seqA):
                pat = seqA[i:i+j]
                if pat not in dic.keys():
                    res = search_all_occurrences(seqB,pat)
                    if len(res) > 0:
                        dic[pat] = len(res)
        i = i+1
    return dic
def test_find_frequency():
    seqA = input("Input sequenceA:")
    seqB = input("Input sequenceB:")
    low = input("Input low:")
    high = input("Input high:")
    print(finds_patterns_frequency(seqA,seqB,int(low),int(high)))

test_find_frequency()



def test_RE():
    seq = input("Input sequence:")
    pat = input("Input pattern (as a regular expression):")
    
    res = find_pattern_re(seq, pat)
    if res >= 0:
        print("Pattern found in position: ", res)
    else:  print("Pattern not found")
    
    all_res = find_all_occurrences_re(seq, pat)
    if len(all_res) > 0:
        print("Pattern found in positions: ", all_res)
    else:  print("Pattern not found")
        
    all_ov = find_all_overlap(seq, pat)
    if len(all_ov) > 0:
        print("Pattern found in positions (overlap): ", all_ov)
    else: 
        print("Pattern not found")


#test_RE()