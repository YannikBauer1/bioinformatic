import pandas as pd
import read_fasta
seq = "AACATT"

def dna_ohe(seq):
    #    A C G T
    l = [0,0,0,0]
    for i in seq:
        if i == "A":
            l[0]=1
        if i == "C":
            l[1]=1
        if i == "G":
            l[2]=1
        if i == "T":
            l[3]=1
    return l
print(dna_ohe(seq))

def word_to_kmer(word, k):
    dic = {}
    for i in range(0,len(word)-k):
        kmer = word[i:i+k]
        if kmer in dic:
            dic[kmer]=dic[kmer]+1
        else:
            dic[kmer]=1
    return dic
print(word_to_kmer(seq,2))

def frequencyCount(seq,word):
    c=0
    for i in range(0,len(seq)-len(word)):
        if word == seq[i:i+len(word)]:
            c=c+1
    return c

def file_to_kmer_table(file_name):
    dnas = read_fasta.read_Fasta(file_name)
    kmers = []
    for i in range(len(dnas)):
        kmers.append((i,word_to_kmer(dnas[i],3)))
    all_keys = []
    for i in kmers:
        all_keys = all_keys + list(i[1].keys())
    all_keys = list(dict.fromkeys(all_keys))
    for i in range(len(dnas)):

    print(all_keys)
file_to_kmer_table("PS00727.fasta")




