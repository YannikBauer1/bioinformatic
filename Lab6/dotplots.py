#############################
### DOTPLOT FUNCTIONS
#############################

# create matrix given dimensions: number of rows and columns filled with zeros
def create_mat(nrows, ncols):
    mat = []
    for i in range(nrows):
        mat.append([])
        for j in range(ncols):
            mat[i].append(0)
    return mat

# basic dotplot algorithm: fills with ones coincident characters
def dotplot(seq1, seq2):
    ''' Create a matrix based on the the input sequences and fill the cells that correspond to a match
    '''
    mat = create_mat(len(seq1),len(seq2))
    for i in range(len(seq1)):
        for j in range(len(seq2)):
            if seq1[i] == seq2[j]:
                mat[i][j] = 1
    return mat

# extended dotplot with window and stringency parameters
def extended_dotplot (seq1, seq2, window, stringency):
    mat = create_mat(len(seq1), len(seq2))
    start = int(window/2)
    for i in range(start,len(seq1)-start):
        for j in range(start, len(seq2)-start):
            matches = 0
            l = j - start
            for k in range(i-start, i+start+1):
                if seq1[k] == seq2[l]: matches += 1
                l += 1
                if matches >= stringency: mat[i][j] = 1
    return mat

# prints dotplot
def print_dotplot(mat, s1, s2):
    ''' Create a function to visualize the matrix
    Print each row in a different line
    if there is a match print the symbol "*" otherwise blankspace
    Print the symbols of seq2 as columns and seq1 as rows
    use the sys.stdout.write(...) to output
    '''
    import sys
    sys.stdout.write(" " + s2+"\n")
    for i in range(len(s1)):
        sys.stdout.write(s1[i])
        for j in mat[i]:
            if j == 0:
                sys.stdout.write(" ")
            elif j == 1:
                sys.stdout.write("*")
        sys.stdout.write("\n")


def test_diagonal_length(mat, istart, jstart):
    # given the starting indices on the row and column
    # check along the diagonal that starts in istart and jstart
    # the longest sub-sequences of matches; return this value
    ini_i = len(mat)-istart
    ini_j = len(mat[0])-jstart
    ini = 0
    if ini_i < ini_j:
        ini = ini_i
    else:
        ini = ini_j
    longest=0
    rang=0
    for i in range(ini):
        if mat[istart+i][jstart+i] == 1:
            rang = rang +1
        else:
            if rang > longest:
                longest = rang
            rang=0
    return longest

def test():
    s1 = "CGATATAGATT"
    s2 = "TATATAGTAT"
    mat1 = dotplot(s1, s2)
    print_dotplot(mat1, s1, s2)
    print(test_diagonal_length(mat1, 2, 3 ))
    print("")
    mat2 = extended_dotplot(s1, s2, 5, 4)
    print_dotplot(mat2, s1, s2)

def test_w_seq():
    s1 = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"
    s2 = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"
    mat1 = dotplot(s1, s2)
    print_dotplot(mat1, s1, s2)
    print(test_diagonal_length(mat1, 2, 3 ))
    print("")
    mat2 = extended_dotplot(s1, s2, 5, 4)
    print_dotplot(mat2, s1, s2)


if __name__ == "__main__":
    #test()
    test_w_seq()