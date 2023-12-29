def scan_word(filename):
    from re import search
    fh = open(filename)
    lines=fh.readlines()
    res=[]
    for l in range(len(lines)):
        regexp = ">.{10}.{0,13}sapiens"
        if search(regexp, lines[l]):
            res.append(l+1)
    fh.close()
    print(res)
    return res
scan_word("PS00727.fasta")