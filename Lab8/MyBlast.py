# -*- coding: utf-8 -*-

class MyBlast:
    '''
    Classe para matrizes de pontos
    '''

    def __init__(self, filename = None, w = 3):
        '''
        Construtor
        '''
        if filename is not None:
            self.db = self.readDatabase(filename)
        else:
            self.db = []
        self.w = w

    def readDatabase(self, filename):
        """From file with sequences line by line read the sequences to a list"""
        fh = open(filename, "r")
        lines = fh.readlines()
        #lines = [i.reaplace("\n","") for i in lines]
        return lines

    def addSequenceDB(self, seq):
        """Add an extra sequence to DB"""
        self.db.append(seq)

    def buildMap (self, query):
        w = self.w
        res = {}
        for i in range(len(query) - w + 1): subseq = query[i:i + w]
        if subseq in res:
            res[subseq].append(i)
        else:
            res[subseq] = [i]
        return res

    def getHits (self, seq, query):
        m = self.buildMap(query)
        w = self.w
        res = []  # list of tuples
        for i in range(len(seq)-w+1):
            subseq = seq[i:i + w]
            if subseq in m:
                l = m[subseq]
                for ind in l:
                    res.append((ind, i))
        print(res)
        return res

    def getHits2 (self, seq, query,mismatches):
        m = self.buildMap(query)
        w = self.w
        res = []  # list of tuples
        if mismatches>w:
            return False
        for i in range(len(seq)-w+1):
            subseq = seq[i:i + w]
            if subseq in m:
                l = m[subseq]
                for ind in l:
                    res.append((ind, i))
            #[for i in list(m.keys())]
        return res

    def extendsHit (self, seq, hit, query):
        stq, sts = hit[0], hit[1]
        matfw = 0
        k=0
        bestk = 0
        while 2*matfw >= k and stq+self.w+k < len(query) and sts+self.w+k < len(seq):
            if query[stq+self.w+k] == seq[sts+self.w+k]:
                matfw+=1
                bestk = k+1
            k += 1
        size = self.w + bestk

        k = 0
        matbw = 0
        bestk = 0
        while 2*matbw >= k and stq > k and sts > k:
            if query[stq-k-1] == seq[sts-k-1]:
                matbw+=1
                bestk = k+1
            k+=1
        size += bestk

        return (stq-bestk, sts-bestk, size, self.w+matfw+matbw)

    def hitBestScore(self, seq, query):
        hits = self.getHits(seq, query)
        bestScore = -1.0
        best = ()
        for h in hits:
            ext = self.extendsHit(seq, h, query)
            score = ext[3]
            if score > bestScore or (score== bestScore and ext[2] < best[2]):
                bestScore = score
                best = ext
        return best


    def bestAlignment (self, query):
        self.buildMap(query)
        bestScore = -1.0
        res = (0,0,0,0,0)
        for k in range(0,len(self.db)):
            bestSeq = self.hitBestScore(self.db[k], query)
            if bestSeq != ():
                score = bestSeq[3]
                if score > bestScore or (score== bestScore and bestSeq[2] < res[2]):
                    bestScore = score
                    res = bestSeq[0], bestSeq[1], bestSeq[2], bestSeq[3], k
        if bestScore < 0: return ()
        else: return res

def test1():
    mb = MyBlast("seqBlast.txt", 11)
    query = "gacgcctcgcgctcgcgcgctgaggcaaaaaaaaaaaaaaaaaaaatcggatagctagctgagcgctcgatagcgcgttcgctgcatcgcgtatagcgctgaagctcccggcgagctgtctgtaaatcggatctcatctcgctctatcct"
    r = mb.bestAlignment(query)
    print(r)
#test1()


def test2():
    mb = MyBlast("seqBlast.txt", 11)
    query2 = "cgacgacgacgacgaatgatg"
    r = mb.bestAlignment(query2)
    print(r)
test2()
