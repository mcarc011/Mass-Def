#%%
import numpy as np
import itertools
from time import time

#%%
def Sduality(X,F,p,Wi):
    W = Wi.copy()
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

    delW = []
    addedmeson = []
    relabels = {}
    for i in range(len(Mn)):
        for j in range(len(Mn)):
            if Mn[i,j] !=0:
                for la in labels:
                    meson = la+str(i+1)+str(j+1)
                    X1 = la+str(j+1)+str(p+1)
                    X2 = la+str(p+1)+str(i+1)
                    if meson not in ','.join(W) and X1 not in ','.join(W) and X2 not in ','.join(W):
                        break 
                delW +=[meson+','+X1+','+X2]
                addedmeson += [meson]

                for l1 in labels:
                    R1 = l1+str(i+1)+str(p+1)
                    for l2 in labels:
                        R2 = l2+str(p+1)+str(j+1)
                        for term in W:
                            save = False
                            term = term.split(',')
                            for wi in range(len(W)):
                                cycleterm = term[wi:]+term[:wi]
                                if R1+','+R2 in ','.join(cycleterm):
                                    save = True
                            if save: 
                                newlabel = ','.join(cycleterm)
                                newlabel = newlabel.replace(R1+','+R2,meson)
                                relabels[','.join(term)] = newlabel
    for wi in range(len(W)):
        if W[wi] in relabels:
            W[wi] = relabels[W[wi]]
    W += delW

    rules = []
    for term in W:
        if len(term.split(','))==2:
            m1,m2 = term.split(',')
            Mn[int(list(m1)[1])-1,int(list(m1)[2])-1] -=1
            Mn[int(list(m2)[1])-1,int(list(m2)[2])-1] -=1
            rules += [derivative(m1,W)]
            rules += [derivative(m2,W)]

    for ri in rules:
        sizes = [len(s.split(',')) for s in ri]
        rit = [(len(s.split(',')),s) for s in ri]
        rit.sort()
        if 1 in sizes:
            wtemp = '+'.join(W)
            wtemp = wtemp.replace(rit[0][1],rit[1][1])
            W = wtemp.split('+')

    # for ai, aterm in enumerate(W):
    #     aterm = aterm.split(',')
    #     for bi,bterm in enumerate(W):
    #         bterm = bterm.split(',')
    #         sameterm = False
    #         if len(aterm)==len(bterm):
    #             sameterm = True
    #             for a in aterm:
    #                 if a not in bterm:
    #                     sameterm = False
    #                     break
    #         if ai!=bi and sameterm:
    #             print(W[ai],W[bi])
    #             W.remove(W[ai])
    #             W.remove(W[bi])
    #             break

    for m in addedmeson:
        if m not in ','.join(W):
            mt=[int(s) for s in m if s.isdigit()]
            Mn[mt[0]-1,mt[1]-1]-=1
            Xn[mt[1]-1,mt[0]-1]-=1

    return Xn,F, (rules,W)

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

def FindPhases(X: np.array,F: np.array,W:list) -> np.array:
    DualityWeb = [(X,F,W)]
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

dweb = Sduality(p4b,p4b-p4b,5,W)
a,b,c = dweb
#dweb2 = Sduality(a,b,5,c)

ansW = 'X13X34X41+X14X47X71+X24X45X52+X12X24X41+X13X37X71+X14X45X51+X23X37Y72+X34X47Y73+Y72X26X67+Y73X36X67+X51X12X26X65+X52X23X36X65'
ansW = ansW.replace('X',',X')
ansW = ansW.replace('Y',',Y')
ansW = [w[1:] for w in ansW.split('+')]
#%%