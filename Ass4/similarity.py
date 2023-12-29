import sys

# Trabalho Grupo D:
#   Paulo Sousa 201805067
#   Sofia Malpique Lopes 201704877
#   Yannik Bauer 201805440



#-------------------------------------------------------------#
#----------------Functions from classes-----------------------#
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

def read_submat_file (filename):
    """read substitution matrix from file """
    sm = {}
    f = open(filename, "r")
    line = f.readline()
    tokens = line.split("\t")
    ns = len(tokens)
    alphabet = []
    for i in range(0, ns):
        alphabet.append(tokens[i][0])
    for i in range(0,ns):
        line = f.readline();
        tokens = line.split("\t");
        for j in range(0, len(tokens)):
            k = alphabet[i]+alphabet[j]
            sm[k] = int(tokens[j])
    return sm

def score_pos(c1, c2, sm, g):
    """score of a position (column)"""
    if c1 == "-" or c2 == "-":
        return g
    else:
        return sm[c1 + c2]

def score_align(seq1, seq2, sm, g):
    """ score of the whole alignment; iterate through the two sequences
    sum the score of each position and return its sum; assume sequences are of equal length

    """
    res = 0
    for i in range(len(seq1)):
        res = res + score_pos(seq1[i], seq2[i], sm, g)
    return res

def score_affinegap(seq1, seq2, sm, g, r):
    ''' calculates the score of alignment based on : affine_gap(len) = g + r*len
    if the gap is open (first occurrence) sum value g; if gap continues sum r to each new gap position;
    if there is no gap use the substitution matrix for the score.
    '''
    res = 0
    ingap1 = False  # two f are true when inside gap sequences
    ingap2 = False
    for i in range(len(seq1)):
        if seq1[i] == "-":
            # gap is already open; add r
            if ingap1:
                res += r
            else:
                # gap is open for the first time; add g
                ingap1 = True
                res += g
        elif seq2[i] == "-":
            # gap is already open; add r
            if ingap2:
                res += r
            else:
                # gap is open for the first time; add g
                ingap2 = True
                res += g
        else:
            # no gaps; use substitution matrix
            if ingap1: ingap1 = False
            if ingap2: ingap2 = False
            res += sm[seq1[i] + seq2[i]]
    return res

def identity(seq1, seq2, alphabet = "ACGT"):
    '''calculate the identity score between seq1 and seq2 '''
    score=0
    for i in range(len(seq1)):
        if seq1[i]==seq2[i]:
            score = score +1
    return score/len(seq1)

def identity2(seq1, seq2, alphabet = "ACGT"):
    '''calculate the identity score between seq1 and seq2 '''
    score=0
    for i in range(len(seq1)):
        if seq1[i]!=seq2[i]:
            score = score +1
    return score

## global alignment
def needleman_Wunsch (seq1, seq2, sm, g):
    """Global Alignment"""
    S = [[0]]
    T = [[0]]
    # initialize gaps in rows
    for j in range(1, len(seq2)+1):
        S[0].append(g * j)
        T[0].append(3)  # horizontal move: 3
    # initialize gaps in cols
    for i in range(1, len(seq1)+1):
        S.append([g * i])
        T.append([2])  # vertical move: 2
    # apply the recurrence to fill the matrices
    for i in range(0, len(seq1)):
        for j in range(len(seq2)):
            s1 = S[i][j] + score_pos (seq1[i], seq2[j], sm, g); # diagonal
            s2 = S[i][j+1] + g  # vertical
            s3 = S[i+1][j] + g # horizontal
            S[i+1].append(max(s1, s2, s3)) # na matrix score add max value
            T[i+1].append(max3t(s1, s2, s3))
    return (S, T)
def max3t (v1, v2, v3):
    """Provides the integer to fill in T"""
    if v1 > v2:
        if v1 > v3: return 1
        else: return 3
    else:
        if v2 > v3: return 2
        else: return 3
def recover_align (T, seq1, seq2):
    # alignment are two strings
    res = ["", ""]
    i = len(seq1)
    j = len(seq2)
    while i>0 or j>0:
        if T[i][j]==1: # diagonal move
            res[0] = seq1[i-1] + res[0] # add to align of seq1 a symbol from seq1(i-1)
            res[1] = seq2[j-1] + res[1] # add to align of seq2 a symbol from seq2(i-1)
            i -= 1
            j -= 1
        elif T[i][j] == 3: # horizontal move
            res[0] = "-" + res[0]   # insert gap na seq 1
            res[1] = seq2[j-1] + res[1] # insert symbol from seq2
            j -= 1
        else: # vertical move
            res[0] = seq1[i-1] + res[0] # insert symbol from seq1
            res[1] = "-" + res[1] # insert gap na seq 2
            i -= 1
    return res

## local alignment
def smith_Waterman (seq1, seq2, sm, g):
    """Local alignment"""
    S = [[0]]
    T = [[0]]
    maxscore = 0
    # first row filled with zero
    for j in range(1, len(seq2)+1):
        S[0].append(0)
        T[0].append(0)
    # first column filled with zero
    for i in range(1, len(seq1)+1):
        S.append([0])
        T.append([0])
    for i in range(0, len(seq1)):
        for j in range(len(seq2)):
            s1 = S[i][j] + score_pos (seq1[i], seq2[j], sm, g);
            s2 = S[i][j+1] + g
            s3 = S[i+1][j] + g
            b = max(s1, s2, s3)
            if b <= 0:
                S[i+1].append(0)
                T[i+1].append(0)
            else:
                S[i+1].append(b)
                T[i+1].append(max3t(s1, s2, s3))
                if b > maxscore:
                    maxscore = b
    return (S, T, maxscore)
def recover_align_local (S, T, seq1, seq2):
    """recover one of the optimal alignments"""
    res = ["", ""]
    """determine the cell with max score"""
    i, j = max_mat(S)
    """terminates when finds a cell with zero"""
    while T[i][j]>0:
        if T[i][j]==1:
            res[0] = seq1[i-1] + res[0]
            res[1] = seq2[j-1] + res[1]
            i -= 1
            j -= 1
        elif T[i][j] == 3:
            res[0] = "-" + res[0];
            res[1] = seq2[j-1] + res[1]
            j -= 1
        elif T[i][j] == 2:
            res[0] = seq1[i-1] + res[0]
            res[1] = "-" + res[1]
            i -= 1
    return res
def max_mat(mat):
    """finds the max cell in the matrix"""
    maxval = mat[0][0]
    maxrow = 0
    maxcol = 0
    for row in range(len(mat)):
        for col in range(len(mat[0])):
            if mat[row][col]>maxval:
                maxval=mat[row][col]
                maxrow=row
                maxcol=col
    return (maxrow,maxcol)

def print_mat2 (mat,ids):
    l=max([len(i) for i in ids])
    print("{x:^{size}}".format(size=l,x=""), end="")
    for i in ids:
        print("{x:^{size}}".format(x=i,size=l),end="")
    print("\n",end="")
    for i in range(len(ids)):
        print("{x:^{size}}".format(x=ids[i], size=l), end="")
        for j in range(len(mat[i])):
            print("{x:^{size}}".format(x=mat[i][j], size=l), end="")
        print("\n", end="")



#-------------------------------------------------------------#
#----------------------Exercicios-----------------------------#
#-------------------------------------------------------------#

inputFile = ""
try:
    inputFile = sys.argv[1].strip()
except:
    sys.exit("Failed because of missing input File!!!")
dic = read_fasta_2dictionary(inputFile)
p1 = "YP_009724390.1"
seq1 = dic[p1]
sm = {}
try:
    sm = read_submat_file("blosum62.mat")
except:
    sm = {'AA': 4, 'AR': -1, 'AN': -2, 'AD': -2, 'AC': 0, 'AQ': -1, 'AE': -1, 'AG': 0, 'AH': -2, 'AI': -1, 'AL': -1, 'AK': -1, 'AM': -1, 'AF': -2, 'AP': -1, 'AS': 1, 'AT': 0, 'AW': -3, 'AY': -2, 'AV': 0, 'RA': -1, 'RR': 5, 'RN': 0, 'RD': -2, 'RC': -3, 'RQ': 1, 'RE': 0, 'RG': -2, 'RH': 0, 'RI': -3, 'RL': -2, 'RK': 2, 'RM': -1, 'RF': -3, 'RP': -2, 'RS': -1, 'RT': -1, 'RW': -3, 'RY': -2, 'RV': -3, 'NA': -2, 'NR': 0, 'NN': 6, 'ND': 1, 'NC': -3, 'NQ': 0, 'NE': 0, 'NG': 0, 'NH': 1, 'NI': -3, 'NL': -3, 'NK': 0, 'NM': -2, 'NF': -3, 'NP': -2, 'NS': 1, 'NT': 0, 'NW': -4, 'NY': -2, 'NV': -3, 'DA': -2, 'DR': -2, 'DN': 1, 'DD': 6, 'DC': -3, 'DQ': 0, 'DE': 2, 'DG': -1, 'DH': -1, 'DI': -3, 'DL': -4, 'DK': -1, 'DM': -3, 'DF': -3, 'DP': -1, 'DS': 0, 'DT': -1, 'DW': -4, 'DY': -3, 'DV': -3, 'CA': 0, 'CR': -3, 'CN': -3, 'CD': -3, 'CC': 9, 'CQ': -3, 'CE': -4, 'CG': -3, 'CH': -3, 'CI': -1, 'CL': -1, 'CK': -3, 'CM': -1, 'CF': -2, 'CP': -3, 'CS': -1, 'CT': -1, 'CW': -2, 'CY': -2, 'CV': -1, 'QA': -1, 'QR': 1, 'QN': 0, 'QD': 0, 'QC': -3, 'QQ': 5, 'QE': 2, 'QG': -2, 'QH': 0, 'QI': -3, 'QL': -2, 'QK': 1, 'QM': 0, 'QF': -3, 'QP': -1, 'QS': 0, 'QT': -1, 'QW': -2, 'QY': -1, 'QV': -2, 'EA': -1, 'ER': 0, 'EN': 0, 'ED': 2, 'EC': -4, 'EQ': 2, 'EE': 5, 'EG': -2, 'EH': 0, 'EI': -3, 'EL': -3, 'EK': 1, 'EM': -2, 'EF': -3, 'EP': -1, 'ES': 0, 'ET': -1, 'EW': -3, 'EY': -2, 'EV': -2, 'GA': 0, 'GR': -2, 'GN': 0, 'GD': -1, 'GC': -3, 'GQ': -2, 'GE': -2, 'GG': 6, 'GH': -2, 'GI': -4, 'GL': -4, 'GK': -2, 'GM': -3, 'GF': -3, 'GP': -2, 'GS': 0, 'GT': -2, 'GW': -2, 'GY': -3, 'GV': -3, 'HA': -2, 'HR': 0, 'HN': 1, 'HD': -1, 'HC': -3, 'HQ': 0, 'HE': 0, 'HG': -2, 'HH': 8, 'HI': -3, 'HL': -3, 'HK': -1, 'HM': -2, 'HF': -1, 'HP': -2, 'HS': -1, 'HT': -2, 'HW': -2, 'HY': 2, 'HV': -3, 'IA': -1, 'IR': -3, 'IN': -3, 'ID': -3, 'IC': -1, 'IQ': -3, 'IE': -3, 'IG': -4, 'IH': -3, 'II': 4, 'IL': 2, 'IK': -3, 'IM': 1, 'IF': 0, 'IP': -3, 'IS': -2, 'IT': -1, 'IW': -3, 'IY': -1, 'IV': 3, 'LA': -1, 'LR': -2, 'LN': -3, 'LD': -4, 'LC': -1, 'LQ': -2, 'LE': -3, 'LG': -4, 'LH': -3, 'LI': 2, 'LL': 4, 'LK': -2, 'LM': 2, 'LF': 0, 'LP': -3, 'LS': -2, 'LT': -1, 'LW': -2, 'LY': -1, 'LV': 1, 'KA': -1, 'KR': 2, 'KN': 0, 'KD': -1, 'KC': -3, 'KQ': 1, 'KE': 1, 'KG': -2, 'KH': -1, 'KI': -3, 'KL': -2, 'KK': 5, 'KM': -1, 'KF': -3, 'KP': -1, 'KS': 0, 'KT': -1, 'KW': -3, 'KY': -2, 'KV': -2, 'MA': -1, 'MR': -1, 'MN': -2, 'MD': -3, 'MC': -1, 'MQ': 0, 'ME': -2, 'MG': -3, 'MH': -2, 'MI': 1, 'ML': 2, 'MK': -1, 'MM': 5, 'MF': 0, 'MP': -2, 'MS': -1, 'MT': -1, 'MW': -1, 'MY': -1, 'MV': 1, 'FA': -2, 'FR': -3, 'FN': -3, 'FD': -3, 'FC': -2, 'FQ': -3, 'FE': -3, 'FG': -3, 'FH': -1, 'FI': 0, 'FL': 0, 'FK': -3, 'FM': 0, 'FF': 6, 'FP': -4, 'FS': -2, 'FT': -2, 'FW': 1, 'FY': 3, 'FV': -1, 'PA': -1, 'PR': -2, 'PN': -2, 'PD': -1, 'PC': -3, 'PQ': -1, 'PE': -1, 'PG': -2, 'PH': -2, 'PI': -3, 'PL': -3, 'PK': -1, 'PM': -2, 'PF': -4, 'PP': 7, 'PS': -1, 'PT': -1, 'PW': -4, 'PY': -3, 'PV': -2, 'SA': 1, 'SR': -1, 'SN': 1, 'SD': 0, 'SC': -1, 'SQ': 0, 'SE': 0, 'SG': 0, 'SH': -1, 'SI': -2, 'SL': -2, 'SK': 0, 'SM': -1, 'SF': -2, 'SP': -1, 'SS': 4, 'ST': 1, 'SW': -3, 'SY': -2, 'SV': -2, 'TA': 0, 'TR': -1, 'TN': 0, 'TD': -1, 'TC': -1, 'TQ': -1, 'TE': -1, 'TG': -2, 'TH': -2, 'TI': -1, 'TL': -1, 'TK': -1, 'TM': -1, 'TF': -2, 'TP': -1, 'TS': 1, 'TT': 5, 'TW': -2, 'TY': -2, 'TV': 0, 'WA': -3, 'WR': -3, 'WN': -4, 'WD': -4, 'WC': -2, 'WQ': -2, 'WE': -3, 'WG': -2, 'WH': -2, 'WI': -3, 'WL': -2, 'WK': -3, 'WM': -1, 'WF': 1, 'WP': -4, 'WS': -3, 'WT': -2, 'WW': 11, 'WY': 2, 'WV': -3, 'YA': -2, 'YR': -2, 'YN': -2, 'YD': -3, 'YC': -2, 'YQ': -1, 'YE': -2, 'YG': -3, 'YH': 2, 'YI': -1, 'YL': -1, 'YK': -2, 'YM': -1, 'YF': 3, 'YP': -3, 'YS': -2, 'YT': -2, 'YW': 2, 'YY': 7, 'YV': -1, 'VA': 0, 'VR': -3, 'VN': -3, 'VD': -3, 'VC': -1, 'VQ': -2, 'VE': -2, 'VG': -3, 'VH': -3, 'VI': 3, 'VL': 1, 'VK': -2, 'VM': 1, 'VF': -1, 'VP': -2, 'VS': -2, 'VT': 0, 'VW': -3, 'VY': -1, 'VV': 4}
g = -8

def ex1(ali):
    sim_scores = {}
    for i in dic.keys():
        if i == p1:
            continue
        seq = dic[i]
        res=()
        if ali=="global":
            res = needleman_Wunsch(seq1, seq, sm, g)
        else:
            res = smith_Waterman(seq1, seq, sm, g)
        S = res[0]
        sim_scores[i] = S[len(seq1)][len(seq)]
    sim_socres = dict(sorted(sim_scores.items(), key=lambda item: item[1], reverse=True))
    print("Similarity list for sequences with ", p1)
    for i in sim_socres.keys():
        print("  {:15} {:6}".format(i, sim_socres[i]))

def ex3(ali):
    score_M = []
    mismatch_M = []
    seqs = list(dic.values())
    for i in range(len(dic)):
        score_M.append([])
        mismatch_M.append([])
        for j in range(len(dic)):
            if j < i:
                res=()
                S=[]
                T=[]
                if ali=="global":
                    res = needleman_Wunsch(seqs[i], seqs[j], sm, g)
                    S = res[0]
                    T = res[1]
                    score = S[len(seqs[i])][len(seqs[j])]
                else:
                    res = smith_Waterman(seqs[i], seqs[j], sm, g)
                    S = res[0]
                    T = res[1]
                    score = res[2]
                alig = recover_align(T, seqs[i], seqs[j])
                mismatches = identity2(alig[0], alig[1])
                score_M[i].append(score)
                mismatch_M[i].append(mismatches)
    print("Similarity matrix:")
    print_mat2(score_M, list(dic.keys()))
    print("")
    print("Mismatch matrix:")
    print_mat2(mismatch_M, list(dic.keys()))
    return score_M, mismatch_M

def matrix_to_dic(M,ids):
    dic={}
    for i in ids:
        dic[i]={}
    for i in range(len(M)):
        for j in range(len(M[i])):
            dic[ids[i]][ids[j]] = M[i][j]
            dic[ids[j]][ids[i]] = M[i][j]
    return dic

def changeMatrix(D,O):
    newM=[]
    for i in range(len(O)):
        newM.append([])
        for j in range(i):
            newM[i].append(D[O[i]][O[j]])
    return newM


### 1
print("-------Start global alignment---------")
print("")
ex1("global")
print("")


### 2
print("Como YP_009724390.1 e QHO60594.1 estão iguais, faz sentido ser a sequencia mais similar.\n"
      "QVT76606.1 também é de 2020 como YP_009724390.1 por isso também faz sentido ser muito similar.\n"
      "Parece que YP_009724390.1 mutou-se de QHR63300.2 de 2013 porque a sequencia é muito parecida.\n"
      "AVP78042.1, AAP13441.1, AAP41037.1, AAU04646.1 também sao parcidas mas parece que sao mutações mais\n"
      "antigas.\n"
      "YP_009047204.1 é muito diferente porque é só um respiratory syndrome related to covid.\n"
      "AMK59677.1 também é do humano como YP_009724390.1 mas é muito diferente e por isso se trata se calhar\n"
      "duma recombinação.")

### 3
score_M, mismatch_M = ex3("global")

### 4
print("")
print("Olhando para as matrizes parece que YP_009724390.1, QHO60594.1, QHR63300.2 e QVT76606.1 são\n"
      "muito parecido um ao outro. YP_009724390.1, QHO60594.1 e QVT76606.1 apareceram em 2020,\n"
      "QHR63300.2 já apareceu em 2013 mas por causa da semalhanca parece que os tres são mutações de\n"
      "QHR63300.2.\n"
      "AAP13441.1, AAU04646.1 e AAP41037.1 também parecem ser muito semelhante um ao outro, o que faz\n"
      "sentido porque os tres são de 2002-2004.\n"
      "AVP78042.1 é mais ou menos parecido com todos mencionados em cima portante é se calhar uma mutação\n"
      "mais antiga.\n"
      "AMK59677.1 e YP_009047204.1 são muito diferente de todos, no caso de YP_009047204.1 é porque é só\n"
      "um Middle East respiratory syndrome-related coronavirus portanto faz sentido ser muito diferente,\n"
      "AMK59677.1 vem se calhar de uma recombinação.")
print("")
print("Matriz ajustado com os grupos relacionados:")
print("")

newOrder=["AAP13441.1","AAU04646.1","AAP41037.1","QHR63300.2","YP_009724390.1","QHO60594.1",
          "QVT76606.1","AVP78042.1","AMK59677.1","YP_009047204.1"]
print("Similarity matrix:")
D=matrix_to_dic(score_M,list(dic.keys()))
newM=changeMatrix(D,newOrder)
print_mat2(newM,newOrder)
print("")
print("Mismatch matrix:")
D=matrix_to_dic(mismatch_M,list(dic.keys()))
newM=changeMatrix(D,newOrder)
print_mat2(newM,newOrder)


### 5
print("")
print("")
print("-------Start local alignment---------")
print("")
ex1("local")
print("")
score_M, mismatch_M = ex3("local")
print("")
print("Matriz ajustado com os grupos relacionados:")
print("")
print("Similarity matrix:")
D=matrix_to_dic(score_M,list(dic.keys()))
newM=changeMatrix(D,newOrder)
print_mat2(newM,newOrder)
print("")
print("Mismatch matrix:")
D=matrix_to_dic(mismatch_M,list(dic.keys()))
newM=changeMatrix(D,newOrder)
print_mat2(newM,newOrder)




