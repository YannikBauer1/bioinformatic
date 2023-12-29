
def ex1(string):
    return "A:"+str(string.count("A"))+" "+\
           "C:"+str(string.count("C"))+" "+\
           "G:"+str(string.count("G"))+" "+\
           "T:"+str(string.count("T"))

print(ex1("ATGCTTCAGAAAGGTCTTACG"))
print(ex1("AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"))

def ex2(list):
    l=[i for i in list if i>0]
    return "min:"+str(min(l))+" "+\
           "avg:"+str(sum(l)/len(l))+" "+\
           "max:"+str(max(l))

print(ex2([4, 5, -6, 6, 10, -15, 0, 5]))