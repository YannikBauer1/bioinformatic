
def pergunta1(k1,k2):
    return (k1**2+k2**2)**0.5
#print(pergunta1(5,10))

def pergunta2(l,u):
    l=int(l)+1
    u=int(u)
    l=[x for x in range(l,u+1)]
    return sum(l)
#print(pergunta2(0.9,7.9))

def pergunta3(s):
    l = [x for x in s]
    l.sort()
    i = 0
    dic={}
    while i !=len(l):
        c=l.count(l[i])
        dic[l[i]]=c
        i=i+c
    return dic
#print(pergunta3("asdfasdfasdzzzy1"))

def pergunta4(s):
    return s.upper()
def pergunta4_2():
    r=input("Insert a string:")
    return r.upper()
pergunta4_2()