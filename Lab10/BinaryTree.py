class BinaryTree:

    def __init__(self, val, dist = 0, left = None, right = None):
        self.value = val
        self.distance = dist
        self.left = left
        self.right = right

    def get_cluster(self):
       # returns a list of all leaves for that tree
       trees=[self]
       l=[]
       size = self.size()[1]
       while len(l)<size:
           trees2=[]
           for i in trees:
               if i.value != -1:
                   l.append(i.value)
               else:
                   trees2 = trees2 + [i.left,i.right]
           trees=trees2
       return l

    def print_tree(self):
        self.print_tree_rec(0, "Root")
    
    def print_tree_rec (self, level, side):
        tabs = ""
        for i in range(level): tabs += "\t"
        if self.value >= 0:
            print(tabs, side, " - value:", self.value)
        else:
            print(tabs, side, "- Dist.: ", self.distance)
            if (self.left != None): 
                self.left.print_tree_rec(level+1, "Left")
            if (self.right != None): 
                self.right.print_tree_rec(level+1, "Right")

    def size(self):
        ''' size of the tree: returns two values
        - number of internal nodes of the tree
        - number of leaves'''
        numleaves = 0
        numnodes = 0
        if self.value >= 0:
            numleaves = 1
        else: 
            if (self.left != None):
                resl = self.left.size()
            else: resl = (0,0)
            if (self.right != None):  
                resr = self.right.size() 
            else: resr = (0,0)
            numnodes += (resl[0] + resr[0] + 1)
            numleaves += (resl[1] + resr[1])
        return numnodes, numleaves

    def exists_leaf(self, leafnum):
        '''
        returns true or false if leafnum appears in the leaves of the tree
        '''
        l = self.get_cluster()
        #print("adsfasdf")
        #print(l)
        if l.count(leafnum)>0:
            return True
        else:
            return False
    
    def common_ancestor(self, leaf1, leaf2):
       ''' Return simplest tree that contains leaf1, leaf2
       '''
       if self.value >=0:
           return None
       if self.left.exists_leaf(leaf1):
            if self.left.exists_leaf(leaf2):
                return self.left.common_ancestor(leaf1,leaf2)
            if self.right.exists_leaf(leaf2):
                return self
            return None
       if self.right.exists_leaf(leaf1):
            if self.right.exists_leaf(leaf2):
                return self.right.common_ancestor(leaf1,leaf2)
            if self.left.exists_leaf(leaf2):
                return self
       return None

    def distance_leaves(self, leafnum1, leafnum2):
       ''' distance between leafnum1 and leafnum2 using the common ancestor function.
       d(leaf1, leaf2) = 2 * height(common_ancest(leaf1, leaf2))
       '''
       a = self.common_ancestor(leafnum1,leafnum2)
       return 2*a.distance

def test():
    # leaf
    a = BinaryTree(1)
    b = BinaryTree(2)
    c = BinaryTree(3)
    d = BinaryTree(4)
    # internal nodes
    e = BinaryTree(-1, 2.0, b, c)
    f = BinaryTree(-1, 1.5, d, a)
    g = BinaryTree(-1, 4.5, e, f)
    g.print_tree()
    print(g.get_cluster())
    # testing exercise 3
    print(g.size())
    print(g.exists_leaf(1))
    print(g.exists_leaf(5))
    g.common_ancestor(1,4).print_tree()
    print(g.distance_leaves(1,4))
    print(g.distance_leaves(1,2))

if __name__ == '__main__':
    test()

def test_aula():
    # leaves
    a = BinaryTree(2)
    b = BinaryTree(4)
    c = BinaryTree(3)
    d = BinaryTree(1)
    # internal nodes
    e = BinaryTree(-1, 2.0, a, b)
    f = BinaryTree(-1, 3.0, e, c)
    g = BinaryTree(-1, 4.0, f, d)
    g.print_tree()