# -*- coding: utf-8 -*-

def read_fasta_2dictionary(filename):
    fh = open(filename, "r")
    lines = fh.readlines()
    dic={}
    seq=""
    id=""
    for i in lines:
        if ">" in i:
            if seq != "":
                dic[id]=seq
                seq=""
                id=i.replace(">","").split(" ")[0]
            else:
                id=i.replace(">","").split(" ")[0]
        else:
            seq=seq+i.replace("\n","")
    dic[id] = seq
    return dic
#print(read_fasta_2dictionary("../Files/PS00727.fasta"))


def read_seq_from_file(filename):
    """ Reads a sequence from a multi-line text file. """
    fh = open(filename, "r")
    lines = fh.readlines()
    seq = ""
    for l in lines:
        seq += l.replace("\n","")
    fh.close()
    return seq

def write_seq_to_file(seq, filename):
    """ Writes a sequence to file. """
    fh = open(filename,"w")
    fh.write(seq)
    fh.close()
    return None


def read_genetic_code_from_file(filename):
    """ Reads the genetic code to a dictionary from a multi-line text file. """
    f = open(filename,"r")
    lines = f.readlines()
    dic={}
    for l in lines:
        h=l.replace("\n","")
        h=h.replace('"','')
        h.split(" ")
        dic[h[0]]=h[1]
    return dic

def validate_dna (dna_seq):
    """ Checks if DNA sequence is valid. Returns True is sequence is valid, or False otherwise. """
    seqm = dna_seq.upper()
    valid = seqm.count("A") + seqm.count("C") + seqm.count("G") + seqm.count("T")
    if valid == len(seqm): return True
    else: return False

def frequency (seq):
    """ Calculates the frequency of each symbol in the sequence. Returns a dictionary. """
    dic = {}
    for s in seq.upper():
        if s in dic: dic[s] += 1
        else: dic[s] = 1
    return dic

def frequencyPerc (seq):
    """ Calculates the frequency of each symbol in the sequence. Returns a dictionary. """
    dic = {}
    for s in seq.upper():
        if s in dic: dic[s] += 1
        else: dic[s] = 1
    for i in dic.keys():
        dic[i]=str(round(dic[i]/len(seq),3)*100)+"%"
    return dic

def gc_content (dna_seq):
    """ Returns the percentage of G and C nucleotides in a DNA sequence. """
    gc_count = 0
    for s in dna_seq:
        if s in "GCgc": gc_count += 1
    return gc_count / len(dna_seq)

def gc_content_subseq (dna_seq, k=100):
    """ Returns GC content of non-overlapping sub-sequences of size k. """
    i=0
    l=[]
    while i<len(dna_seq):
        gc_c=gc_content(dna_seq[i:i+k])
        l.append(gc_c)
        i=i+k
    return l


def transcription (dna_seq):
    """ Function that computes the RNA corresponding to the transcription of the DNA sequence provided. """
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    return dna_seq.upper().replace("T","U")


def reverse_complement (dna_seq):
    """ Computes the reverse complement of the inputted DNA sequence. """
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    comp_rev = ""
    dic={"A":"T","T":"A","C":"G","G":"C"}
    for i in range(len(dna_seq)-1,-1,-1):
        comp_rev=comp_rev+dic[dna_seq[i]]
    return comp_rev
    # A > T
    # T > A
    # C > G
    # G > C


def translate_codon (cod):
    """Translates a codon into an aminoacid using an internal dictionary with the standard genetic code."""
    tc = {"GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
      "TGT":"C", "TGC":"C",
      "GAT":"D", "GAC":"D",
      "GAA":"E", "GAG":"E",
      "TTT":"F", "TTC":"F",
      "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
      "CAT":"H", "CAC":"H",
      "ATA":"I", "ATT":"I", "ATC":"I",
      "AAA":"K", "AAG":"K",
      "TTA":"L", "TTG":"L", "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
      "ATG":"M", "AAT":"N", "AAC":"N",
      "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
      "CAA":"Q", "CAG":"Q",
      "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R", "AGA":"R", "AGG":"R",
      "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S", "AGT":"S", "AGC":"S",
      "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
      "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
      "TGG":"W",
      "TAT":"Y", "TAC":"Y",
      "TAA":"_", "TAG":"_", "TGA":"_"}
    if cod in tc: return tc[cod]
    else: return None




def translate_seq (dna_seq, ini_pos = 0):
    """ Translates a DNA sequence into an aminoacid sequence. """
    # ini_pos = 0 > frame 1
    # ini_pos = 1 > frame 2
    # ....
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    seqm = dna_seq.upper()
    seq_aa = ""
    i=ini_pos
    while (i<len(dna_seq)):
        s=seqm[i:i+3]
        seq_aa=seq_aa+str(translate_codon(s))
        i=i+3
    return seq_aa


def codon_usage(dna_seq, aa, ini_pos = 0):
    """Provides the frequency of each codon encoding a given aminoacid, in a DNA sequence ."""
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    seqm = dna_seq.upper()
    dic = {}
    total = 0
    for i in range(ini_pos, len(seqm)-2, 3):
        cod = seqm[i:i+3]
        if translate_codon(cod) == aa:
            if cod in dic:
                dic[cod] += 1
            else: dic[cod] = 1
            total += 1
    if total >0:
        for k in dic:
            dic[k] /= total
    return dic

def codon_usage2(dna_seq, aa, ini_pos = 0):
    """Provides the frequency of each codon encoding a given aminoacid, in a DNA sequence ."""
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    seqm = dna_seq.upper()
    dic = {}
    for i in range(ini_pos, len(seqm)-2, 3):
        cod = seqm[i:i+3]
        if translate_codon(cod) == aa:
            if cod in dic:
                dic[cod] += 1
            else: dic[cod] = 1
    return dic


def reading_frames (dna_seq):
    """Computes the six reading frames of a DNA sequence (including the reverse complement."""
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    res = []
    res.append(translate_seq(dna_seq,0))
    res.append(translate_seq(dna_seq,1))
    res.append(translate_seq(dna_seq,2))
    rc = reverse_complement(dna_seq)
    res.append(translate_seq(rc,0))
    res.append(translate_seq(rc,1))
    res.append(translate_seq(rc,2))
    return res

def reading_frames2 (dna_seq):
    """Computes the six reading frames of a DNA sequence (including the reverse complement."""
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    res = []
    print(dna_seq)
    res.append((translate_seq(dna_seq,0),1))
    res.append((translate_seq(dna_seq,1),2))
    res.append((translate_seq(dna_seq,2),3))
    rc = reverse_complement(dna_seq)
    res.append((translate_seq(rc,0),-1))
    res.append((translate_seq(rc,1),-2))
    res.append((translate_seq(rc,2),-3))
    return res


def all_proteins_rf (aa_seq):
    """Computes all posible proteins in an aminoacid sequence."""
    aa_seq = aa_seq.upper()
    current_prot = []
    proteins = []
    for aa in aa_seq:
        if aa == "_":
            if current_prot:
                for p in current_prot:
                    proteins.append(p)
                current_prot = []
        else:
            if aa == "M":
                current_prot.append("")
            for i in range(len(current_prot)):
                current_prot[i] += aa
    return proteins

def all_proteins_rf2 (aa_seq):
    """Computes all posible proteins in an aminoacid sequence."""
    aa_seq = aa_seq.upper()
    current_prot = []
    proteins = []
    for aa in aa_seq:
        if aa == "_":
            if current_prot:
                proteins.append(current_prot[0])
                current_prot = []
        else:
            if aa == "M":
                current_prot.append("")
            for i in range(len(current_prot)):
                current_prot[i] += aa
    return proteins

def all_orfs (dna_seq):
    """Computes all possible proteins for all open reading frames."""
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    rfs = reading_frames(dna_seq)
    res = []
    for i in rfs:
        for j in all_proteins_rf(i):
            res.append(j)
    return res

def all_orfs2 (dna_seq):
    """Computes all possible proteins for all open reading frames."""
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    rfs = reading_frames2(dna_seq)
    res = []
    for i in rfs:
        for j in all_proteins_rf(i[0]):
            res.append((j,i[1]))
    return res

def all_orfs_ord (dna_seq, minsize = 0):
    """Computes all possible proteins for all open reading frames. Returns ordered list of proteins with minimum size."""
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    rfs = reading_frames(dna_seq)
    res = []
    for i in rfs:
        for j in all_proteins_rf(i):
            if len(j)>minsize:
                insert_prot_ord(j,res)
    return res

def all_orfs_ord2 (dna_seq, minsize = 0):
    """Computes all possible proteins for all open reading frames. Returns ordered list of proteins with minimum size."""
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    rfs = reading_frames2(dna_seq)
    res = []
    for i in rfs:
        for j in all_proteins_rf(i[0]):
            if len(j)>minsize:
                insert_prot_ord2((j,i[1]),res)
    return res

def all_orfs_ord3 (dna_seq, minsize = 0):
    """Computes all possible proteins for all open reading frames. Returns ordered list of proteins with minimum size."""
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    rfs = reading_frames2(dna_seq)
    res = []
    for i in rfs:
        for j in all_proteins_rf2(i[0]):
            if len(j) > minsize:
                insert_prot_ord2((j, i[1]), res)
    return res

def insert_prot_ord (prot, list_prots):
    ''' inserts prot in list_prots in a sorted way '''
    i = 0
    while i < len(list_prots) and len(prot) < len(list_prots[i]):
        i += 1
    list_prots.insert(i, prot)

def insert_prot_ord2(prot, list_prots):
    ''' inserts prot in list_prots in a sorted way '''
    i = 0
    while i < len(list_prots) and len(prot[0]) < len(list_prots[i][0]):
        i += 1
    list_prots.insert(i, prot)


def test_sequence_from_file():
    dna_seq = input("DNA sequence file:")
    seq = read_seq_from_file(dna_seq)
    if validate_dna(seq):
        print("translates:", translate_seq(seq))
        print("reverse complemente:", reverse_complement(dna_seq))
        print("Gc content:", gc_content(dna_seq))
    else: print("DNA sequence is not valid")


def test_frequency():
    seq_aa = input("Protein sequence:")
    freq_aa = frequency(seq_aa)
    list_f = sorted(freq_aa.items(), key=lambda x: x[1], reverse = True)
    for (k,v) in list_f:
        print("Aminoacid:", k, ":", v)

def test_all():
    seq = input("Insert DNA sequence:")
    if validate_dna (seq):
        print ("Valid sequence")
        print ("Transcription: ", transcription (seq))
        print("Reverse complement:", reverse_complement(seq))
        print("GC content (global):", gc_content(seq))
        print("Direct translation:" , translate_seq(seq))
        print("All proteins in ORFs (decreasing size): ", all_orfs_ord(seq))
    else: print("DNA sequence is not valid")

def test_files():
    fname = input("Insert input filename:")
    seq = read_seq_from_file(fname)
    if validate_dna (seq):
        print ("Valid sequence")
        print ("Transcription: ", transcription (seq))
        print("Reverse complement:", reverse_complement(seq))
        print("GC content (global):", gc_content(seq))
        print("Direct translation:" , translate_seq(seq))
        print("All proteins in ORFs (decreasing size): ", all_orfs_ord(seq))
        #orfs = all_orfs_ord(seq)
        #i = 1
        #for orf in orfs:
        #    write_seq_to_file(orf, "orf-"+str(i)+".txt")
        #    i += 1
    else: print("DNA sequence is not valid")



def test_aula3():
    #fname = input("Insert input filename:")
    fname= "example_Hinfluenzae.txt"
    seq = read_seq_from_file(fname)
    if validate_dna(seq):
        print("All proteins in ORFs (decreasing size): ", all_orfs_ord2(seq))
        print(longest_Protein(seq))

def longest_Protein(dna_seq):
    prots=all_orfs_ord(dna_seq)
    return prots[0]

if __name__ == "__main__":
    # test here all implemented functions
    # used your own defined sequences or read from example files
    seq = "ATGAGCGAC"
    #print("Reverse:", seq[::-1])
    #print("Reverse complement:", reverse_complement(seq))

    ## uncomment the test function to run
    #test_frequency()
    #test_sequence_from_file()
    #test_all()
    #test_files()
    #test_aula3()

