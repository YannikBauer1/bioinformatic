import sys

# Trabalho Grupo D:
#   Paulo Sousa 201805067
#   Sofia Malpique Lopes 201704877
#   Yannik Bauer 201805440


#-------------------------------------------------------------#
#------------Functions of sequence_functions------------------#
#-------------------------------------------------------------#

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

def read_csv(filename):
    fh = open(filename, "r")
    lines = fh.readlines()
    tab = []
    for l in lines:
        tab += [l.replace("\n","").replace('"',"").split(",")]
    fh.close()
    return tab



#-------------------------------------------------------------#
#----------------------Exercicios-----------------------------#
#-------------------------------------------------------------#

inputFile1 = sys.argv[1].strip()
inputFile2 = ""
try:
    inputFile2 = sys.argv[2].strip()
except:
    inputFile2 = "proteins_86693_757732.csv"

### B
dic = read_fasta_2dictionary(inputFile1)
accession = [*dic][0]
seq = dic[accession]
# B_1
print("Length of the sequence: ",len(seq))
# B_2
freq = frequencyPerc(seq)
print("Frequency (in %) of A, C, G, T:")
for i in freq.keys():
    print("       {} : {}".format(i, freq[i]))
# B_3
print("GC-Content: ", round(gc_content(seq),2))
# B_4
print("Number of start (AUG) codons:")
startCodons = [(codon_usage2(seq,"M",i),i+1) for i in [0,1,2]] + \
              [(codon_usage2(reverse_complement(seq),"M",i),-i-1) for i in [0,1,2]]
for i in startCodons:
    print("       ReadingFrame {} : {}".format(i[1], i[0]["ATG"]))
# B_5
print("Number of stop codons (UAA, UAG, UGA) :")
stopCodons = [(codon_usage2(seq,"_",i),i+1) for i in [0,1,2]] + \
              [(codon_usage2(reverse_complement(seq),"_",i),-i-1) for i in [0,1,2]]
for i in stopCodons:
    for j in i[0].keys():
        print("       ReadingFrame {} ({}) : {}".format(i[1],j,i[0][j]))
# B_6
allCodons = [([codon_usage2(seq,j,i) for j in
               ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","X","Y","_"]], i+1) for i in [0,1,2]]+ \
            [([codon_usage2(reverse_complement(seq),j,i) for j in
               ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","X","Y","_"]], -i-1) for i in [0,1,2]]
dic={}
for i in range(len(allCodons)):
    for j in allCodons[i][0]:
        dic=dict(**dic,**j)
    allCodons[i]=[dic,allCodons[i][1]]
    dic={}
for i in range(len(allCodons)):
    allCodons[i][0]=dict(sorted(allCodons[i][0].items(), key=lambda item: item[1]))
print("Less frequent codons:")
for i in allCodons:
    print("       ReadingFrame {} ({}) : {}".format(i[1], [*i[0]][0], i[0][[*i[0]][0]]))
for i in range(len(allCodons)):
    allCodons[i][0]=dict(sorted(allCodons[i][0].items(), key=lambda item: item[1],reverse=True))
print("Most frequent codons:")
for i in allCodons:
    print("       ReadingFrame {} ({}) : {}".format(i[1], [*i[0]][0], i[0][[*i[0]][0]]))


### C
all_prots=all_orfs_ord3(seq,120)
all_prots=[i[0] for i in all_prots]
# 7
all_prots_seq= str(all_prots)\
    .replace("[","")\
    .replace("]","")\
    .replace("'","")\
    .replace(",","\n")\
    .replace(" ","")
write_seq_to_file(all_prots_seq,"all_potential_proteins.txt")
# 8
all_prots2=all_orfs_ord3(seq,120)
string=""
for i in all_prots2:
    seq_trans = translate_seq(seq,i[1]-1)
    s2 = seq_trans.find(i[0])
    string=string+str(s2*3+i[1])+","+str(len(seq_trans)+s2*3+i[1])+","+str(i[1])+"\n"
write_seq_to_file(string,"orf_coordinates.txt")

def overlap(a,b):
    startA = int(a[0])
    stopA = int(a[1])
    lenA = stopA-startA
    startB = int(b[0])
    stopB = int(b[1])
    if (startA > stopB) or (stopB < startA):
        return 0
    elif startA <= startB and stopB >= stopA:
        return (stopA-startB)/lenA
    elif startA <= startB and stopB < stopA:
        return (stopB - startB) / lenA
    elif startA > startB and stopA > stopB:
        return (stopB - startA)/lenA
    elif startA > startB and stopA < stopB:
        return 1

csv=read_csv(inputFile2)[1:]
csv=[[csv[x][2],csv[x][3],csv[x][6],seq[int(csv[x][2]):int(csv[x][3])]] for x in range(len(csv))]

orf=read_csv("orf_coordinates.txt")
orf=[i+[seq[int(i[0]):int(i[1])]]for i in orf]

print("Max Overlaps for each protein in proteins(...).csv:")
for i in csv:
    l=[]
    for j in orf:
        ol=round(overlap(i,j)*100,2)
        l = l+[ol]
    ma=max(l)
    l=[]
    print("       {:7}  {:3} %".format(i[2],ma))







