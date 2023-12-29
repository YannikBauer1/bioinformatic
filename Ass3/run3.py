
def read_seq_from_file(filename):
    """ Reads a sequence from a multi-line text file. """
    fl = open(filename, "r")
    lines = fl.readlines()
    seq = ""
    for l in lines[1:]:
        seq += l.replace("\n", "")
    fl.close()
    return seq


def write_to_file(words, filename):
    fl = open(filename, "a")
    fl.write(str(words))
    fl.write("\n")
    fl.close()

    return None


def read_genetic_code_from_file(filename):
    """ Reads the genetic code to a dictionary from a multi-line text file. """
    dic = {}
    fl = open(filename, "r")
    lines = fl.readlines()
    for line in lines:
        tupulos = line.strip().replace("\n", "").split(" ")
        dic[tupulos[0]] = tupulos[1]

    fl.close()
    return dic


def validate_dna(dna_seq):
    """ Checks if DNA sequence is valid. Returns True is sequence is valid, or False otherwise. """
    seqm = dna_seq.upper()
    valid = seqm.count("A") + seqm.count("C") + \
        seqm.count("G") + seqm.count("T")
    if valid == len(seqm):
        return True
    else:
        return False


def frequency(seq):
    """ Calculates the frequency of each symbol in the sequence. Returns a dictionary. """
    dic = {}
    for s in seq.upper():
        if s in dic:
            dic[s] += 1
        else:
            dic[s] = 1
    return dic


def gc_content(dna_seq):
    """ Returns the percentage of G and C nucleotides in a DNA sequence. """
    gc_count = 0
    for s in dna_seq:
        if s in "GCgc":
            gc_count += 1
    return gc_count / len(dna_seq)


def gc_content_subseq(dna_seq, k=100):
    """ Returns GC content of non-overlapping sub-sequences of size k. """
    res = []

    for i in range(0, len(dna_seq)-k+1, k):
        subseq = dna_seq[i:i+k]
        gc = gc_content(subseq)
        res.append(gc)

    return res


def transcription(dna_seq):
    """ Function that computes the RNA corresponding to the transcription of the DNA sequence provided. """
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    return dna_seq.upper().replace("T", "U")


def reverse_complement(dna_seq):
    """ Computes the reverse complement of the inputted DNA sequence. """
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    comp = ""
    dna_seq.upper()
    comp = dna_seq[::-1]
    rev_seq = ""
    for i in comp:
        if i == 'A':
            rev_seq += 'T'
        elif i == 'T':
            rev_seq += 'A'
        elif i == 'C':
            rev_seq += 'G'
        elif i == 'G':
            rev_seq += 'C'

    return rev_seq


def translate_codon(cod):
    """Translates a codon into an aminoacid using an internal dictionary with the standard genetic code."""
    tc = {"GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
          "TGT": "C", "TGC": "C",
          "GAT": "D", "GAC": "D",
          "GAA": "E", "GAG": "E",
          "TTT": "F", "TTC": "F",
          "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
          "CAT": "H", "CAC": "H",
          "ATA": "I", "ATT": "I", "ATC": "I",
          "AAA": "K", "AAG": "K",
          "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
          "ATG": "M",
          "AAT": "N", "AAC": "N",
          "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
          "CAA": "Q", "CAG": "Q",
          "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
          "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
          "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
          "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
          "TGG": "W",
          "TAT": "Y", "TAC": "Y",
          "TAA": "_", "TAG": "_", "TGA": "_"}
    if cod in tc:
        return tc[cod]
    else:
        return None


def translate_seq(dna_seq, init_pos=0):
    """ Translates a DNA sequence into an aminoacid sequence. """
    # init_pos = 0 > frame 1
    # init_pos = 1 > frame 2
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    seqm = dna_seq.upper()
    seq_aa = ""
    for i in range(init_pos, len(dna_seq)-2, 3):
        codon = seqm[i:i+3]
        seq_aa += translate_codon(codon)

    return seq_aa


def codon_usage(dna_seq, aa):
    """Provides the frequency of each codon encoding a given aminoacid, in a DNA sequence ."""
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    seqm = dna_seq.upper()
    dic = {}
    total = 0
    for i in range(0, len(seqm)-2, 3):
        cod = seqm[i:i+3]
        if translate_codon(cod) == aa:
            if cod in dic:
                dic[cod] += 1
            else:
                dic[cod] = 1
            total += 1
    if total > 0:
        for k in dic:
            dic[k] /= total
    return dic


def reading_frames(dna_seq):
    """Computes the six reading frames of a DNA sequence (including the reverse complement."""
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    res = []
    res.append(translate_seq(dna_seq, 0))
    res.append(translate_seq(dna_seq, 1))
    res.append(translate_seq(dna_seq, 2))
    #rc = reverse_complement(dna_seq)
    # res.append(translate_seq(rc,0))
    # res.append(translate_seq(rc,1))
    # res.append(translate_seq(rc,2))
    return res


def all_proteins_rf(aa_seq):
    """Computes all possible proteins in an aminoacid sequence."""
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


def insert_prot_ord(prot, list_prots):
    ''' inserts prot in list_prots in a sorted way '''
    i = 0
    while i < len(list_prots) and len(prot[0]) < len(list_prots[i][0]):
        i += 1
    list_prots.insert(i, prot)


def all_orfs_ord(dna_seq, minsize=120):
    """Computes all possible proteins for all open reading frames. Returns ordered list of proteins with minimum size."""
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    rfs = reading_frames(dna_seq)
    res = []
    for rf in rfs:
        prots = all_proteins_rf(rf)
        for p in prots:
            if len(p) > minsize:
                insert_prot_ord(p, res)
    return res


# Exercício B
seqToUse = read_seq_from_file('sequence.fasta')

# (1) - Length
toWrite = "The size of the given sequence is " + str(len(seqToUse)) + "."
print(toWrite)
#write_to_file(toWrite, "outputFile.txt")

# (2) - Frequency (in %) of	A, C, G, T
print(toWrite)
#write_to_file(frequency(seqToUse), "outputFile.txt")

# (3) - GC content
print(toWrite)
#write_to_file(gc_content(seqToUse), "outputFile.txt")

# (4) - Number of Start	(AUG) codons found

rfs = reading_frames(seqToUse)

lst = []
for rf in rfs:
    lst.append(frequency(rf))

toWrite = "START CODONS - rf1:" + \
    str(lst[0]['M']) + " | rf2:" + \
    str(lst[1]['M']) + " | rf3:" + str(lst[2]['M'])
#write_to_file(toWrite, "outputFile.txt")
print(toWrite)

# (5) - Number of Stop Codons (UAA, UAG, UGA)
toWrite = "STOP CODONS - rf1:" + \
    str(lst[0]['_']) + " | rf2:" + \
    str(lst[1]['_']) + " | rf3:" + str(lst[2]['_'])
#write_to_file(toWrite, "outputFile.txt")
print(toWrite)

# (6) - Most and less frequentcodon.
#print("lst: ", lst)
lmin = []
lmax = []
for dic in lst:
    minimum = min(dic.values())
    maximum = max(dic.values())
    for codon, freq in dic.items():
        if freq == minimum:
            lmin.append(codon)
        if freq == maximum:
            lmax.append(codon)

toWrite = "rf1: most - " + str(lmax[0]) + " less - " + str(lmin[0]) + \
    " | rf2: most - " + str(lmax[1]) + " less - " + str(lmin[1]) + \
    " | rf3: most - " + str(lmax[2]) + " less - " + str(lmin[2])

#write_to_file(toWrite, "outputFile.txt")
print(toWrite)


# (7) - output a file named "all_potential_proteins.txt"
'''
def longest_protein(dna_seq):
    prots = all_orfs_ord(dna_seq)
    return prots[0]


print("x:", longest_protein(seqToUse))
'''

write_to_file(toWrite,"all_potential_proteins.txt")


