
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

def find_zync_finger(seq):
    from re import search
    regexp = "C.H.[LIVMFY]C.{2}C[LIVMYA]"
    mo = search(regexp, seq)
    if (mo != None):
        return mo.span()
    else:
        return -1
    
def find_prosite(seq, profile):
    from re import search
    regexp = profile.replace("-","")
    regexp = regexp.replace("x",".")
    regexp = regexp.replace("(","{")
    regexp = regexp.replace(")","}")
    mo = search(regexp, seq)
    if (mo != None):
        return mo.span()[0]
    else:
        return -1
    
def test():
    seq = "HKMMLASCKHLLCLKCIVKLG"
    print(find_zync_finger(seq))
    print(find_prosite(seq,"C-x-H-x-[LIVMFY]-C-x(2)-C-[LIVMYA]"))

def task2():
    dic = read_fasta_2dictionary("Q8RXD4.fasta")
    seq = dic[[*dic][0]]
    newSeq=""
    i =0
    while i < len(seq):
        if find_zync_finger(seq[i:])==-1:
            break
        a,b = find_zync_finger(seq[i:])
        newSeq=seq[i:a].lower()+seq[a:b]
        i=b
    newSeq=newSeq+seq[i:].lower()
    print(newSeq)


task2()
#test()