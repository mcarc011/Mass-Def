#%%
import numpy as np
import itertools
from time import time

#%%
def Sduality(X,F,p,W):
    ProjV = np.zeros(len(X))
    ProjM = np.meshgrid(ProjV,ProjV)[0]
    ProjM[p,p] = 1
    Xn = X - np.dot(X,ProjM) + np.transpose(np.dot(X,ProjM)) - np.transpose(np.dot(np.transpose(X),ProjM))
    Xn += np.dot(np.transpose(X),ProjM)
    Mn = np.transpose(np.dot(Xn,np.dot(ProjM,Xn))) 

    labels = ['X','Y','Z','M','A','B','C','D','E','F','G']
    def derivative(x,f):
        df = []
        for term in f:
            if x in term:
                term = term.split(',')
                term.remove(x)
                df += [','.join(term)]
        return df
    
    def redefine(r,w,switch=False):
        wt = w.copy()
        if len(r[0].split(','))==1:
            m = r[0]
            t = r[1]
        if len(r[1].split(','))==1:
            m = r[1]
            t = r[0]
        if len(t.split(','))==2:
            R1,R2 = t.split(',')
        for wi,term in enumerate(wt):
            tsplit = term.split(',')
            cycle = False
            for pi in range(len(tsplit)):
                tempterm = tsplit[-pi:] + tsplit[:-pi]
                if R1+R2 in ''.join(tempterm):
                    cycle = True
                    break
                if switch and R2+R1 in ''.join(tempterm):
                    R2,R1 = R1,R2
                    cycle = True
                    break
            if cycle:
                tempterm.remove(R2)
                nterm = ','.join(tempterm)
                wt[wi] = nterm.replace(R1,m)
        return wt


    for i in range(len(Mn)):
        for j in range(len(Mn)):
            if Mn[i,j] !=0:
                for l in labels:
                    meson = l+str(i+1)+str(j+1)
                    X1 = l+str(j+1)+str(p+1)
                    X2 = l+str(p+1)+str(i+1)
                    if meson and X1 and X2 not in ','.join(W):
                        break 
                W +=[meson+','+X1+','+X2]

                for l in labels:
                    R1 = l+str(i+1)+str(p+1)
                    R2 = l+str(p+1)+str(j+1)
                    W = redefine([R1+','+R2,meson],W)

    for wi,term in enumerate(W):
        if len(term.split(','))==2:
            m1,m2 = term.split(',')
            W = redefine(derivative(m1,W),W,switch=True)
            W = redefine(derivative(m2,W),W,switch=True)
                
    return Xn,F,W

W = wp4b.replace('−','+')
W = W.replace('\n',' ')
W = W.replace(' ','')
W = W.replace('X',',X')
W = [w[1:] for w in W.split('+')]
Xt,Ft,W = Sduality(p4b,p4b-p4b,5,W)
W
#%%

def Swap(M:np.array, t:tuple):
	Mt = M.copy()
	Mt[t[0]],Mt[t[1]] = M[t[1]],M[t[0]]
	Mt = np.transpose(Mt)
	Mc = Mt.copy()
	Mt[t[0]], Mt[t[1]] = Mc[t[1]], Mc[t[0]]
	return np.transpose(Mt)

def TupFind(L):
	maps = []
	n = 0
	while n!=len(L)-1:
		mt = []
		for m in range(n,len(L)):
			mt += [(n,m)]
		maps += [mt]
		n+=1
	templist = list(itertools.product(*maps))
	return [[(t[0],t[1]) for t in tmap if t[0]!=t[1]] for tmap in templist]

def FindPhases(X: np.array,F: np.array) -> np.array:
    DualityWeb = [(X,F)]
    TrialityMaps = []

    def equivalent(X1, F1, X2, F2,counter=False):
        if np.array_equal(X1,X2) and np.array_equal(F1,F2):
            return True

        ANodes = [(sorted(X1[node]),sorted(np.transpose(X1)[node]),sorted(F1[node])) for node in range(len(X1))]
        BNodes = [(sorted(X2[node]),sorted(np.transpose(X2)[node]),sorted(F2[node])) for node in range(len(X1))]

        try:
            Xt, Ft = X1.copy(), F1.copy()
            for i in range(len(X1)):
                if BNodes[i] != ANodes[i]:
                    b = BNodes[i]
                    j = ANodes[i:].index(b)+i
                    Xt,Ft = Swap(Xt,(j,i)),Swap(Ft,(j,i))
                    ANodes = [(sorted(Xt[node]),sorted(np.transpose(Xt)[node]),sorted(Ft[node])) for node in range(len(X1))]
        except:
            return False 

        if np.array_equal(Xt,X2) and np.array_equal(Ft,F2):
            return True
        
        #final test
        Aswaps = {}
        ANodes = [(sorted(Xt[node]),sorted(np.transpose(Xt)[node]),sorted(Ft[node])) for node in range(len(X1))]
        for i,a in enumerate(ANodes):
            for j,b in enumerate(ANodes):
                if a==b and i!=j:
                    if str(a) not in Aswaps:
                        Aswaps[str(a)] = [i,i,j]
                    if i not in Aswaps[str(a)]:
                        Aswaps[str(a)] += [i]
                    if j not in Aswaps[str(a)]:
                        Aswaps[str(a)] += [j]

        temptylist = [val for val in Aswaps.values()]
        temptlist = [] 
        for tem in temptylist:
            mapto = TupFind(np.array(tem))
            tp = [] 
            for mapper in mapto:
                tp += [[[tem[m[0]],tem[m[1]]] for m in mapper]]
            temptlist += [tp]

        for tlist in itertools.product(*temptlist):
            Xp,Fp = Xt.copy(),Ft.copy()
            for tlistswap in itertools.permutations(tlist):
                for titer in tlistswap: 
                    for step in titer:
                        t = (step[0],step[1])
                        Xp,Fp = Swap(Xp,t),Swap(Fp,t)
                        if np.array_equal(Xp,X2) and np.array_equal(Fp,F2):
                            return True
        return False

    def inweb(xt, ft, dweb, findn=False):
        for d,phase in enumerate(dweb):
            e1 = equivalent(xt,ft,phase[0],phase[1])
            e2 = equivalent(np.transpose(xt),ft,phase[0],phase[1])
            if e1 or e2:# or e3 or e4:
                if findn:
                    return d
                return True
        return False
    
    def anomalies(x1,f1):
        for i in range(len(x1)):
            chirals = np.sum(x1[i]) + np.sum(np.transpose(x1)[i])
            fermis = np.sum(f1[i])
            if chirals - fermis !=2:
                return True
        return False

    for phase in DualityWeb:
        Xt,Ft,Wt = phase[0].copy(),phase[1].copy(),phase[2].copy()

        print('Dweb Length: '+str(len(DualityWeb)),end="\r")

        if len(DualityWeb)>100:
            break
        
        # if anomalies(Xt,Ft):
        #     print('Anomaly Found')
        #     break

        for n in range(len(Xt)):
            if np.sum(np.transpose(Xt)[n])+np.sum(Xt[n])==4:
                Xi,Fi,Wi = Sduality(Xt, Ft, n, Wt)
                if np.trace(Xi)==0 and np.trace(Fi)==0:
                    TrialityMaps += [[(Xt,Ft),(Xi,Fi),n+1]]
                    if not inweb(Xi,Fi,DualityWeb):
                        DualityWeb += [(Xi,Fi,Wi)]

        

    TrialityTuples = [(inweb(t[0][0],t[0][1],DualityWeb,findn=True),
        inweb(t[1][0],t[1][1],DualityWeb,findn=True)) for t in TrialityMaps]
    print('\n'+str(len(DualityWeb))+' Phase\n')
    return DualityWeb,TrialityMaps,TrialityTuples

wp4b = '''X1 2X2 5X5 1 + X1 3X3 4X4 1 + X1 4X4 7X7 1 + X2 4X4 5X5 2 + X3 5X5 6X6 3
− X1 2X2 4X4 1 − X1 3X3 7X7 1 − X1 4X4 5X5 1 − X2 3X3 5X5 2 − X2 5X5 6X6 2
+ X2 3X3 7X7 6X6 2 − X3 4X4 7X7 6X6 3'''

test =np.array([
[0,1,1,0,0,0],
[0,0,1,1,0,0],
[0,0,0,1,1,0],
[0,0,0,0,1,1],
[1,0,0,0,0,1],
[1,1,0,0,0,0]
])

testf = test-test

arrinfo  = []
a = True
for t in wp4b:
    try:
        int(t)
        if a:
            a = False
            aval = int(t)
        else:
            a = True
            arrinfo += [[aval,int(t)]]
    except:
        pass

n = max(max(arrinfo))
incmatrix = [[0 for j in range(n)] for i in range(n)]
for arr in arrinfo:
    incmatrix[arr[0]-1][arr[1]-1] =1
p4b = np.array(incmatrix)

W = wp4b.replace('−','+')
W = W.replace('\n',' ')
W = W.replace(' ','')
W = W.replace('X',',X')
W = [w[1:] for w in W.split('+')]

dweb = FindPhases(p4b,p4b-p4b)

#%%