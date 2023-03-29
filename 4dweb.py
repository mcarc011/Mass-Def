import numpy as np
import itertools
from time import time

def Sduality(X,F,p):
    ProjV = np.zeros(len(X))
    ProjM = np.meshgrid(ProjV,ProjV)[0]
    ProjM[p,p] = 1
    Xn = X - np.dot(X,ProjM) + np.transpose(np.dot(X,ProjM)) - np.transpose(np.dot(np.transpose(X),ProjM))
    Xn += np.dot(np.transpose(X),ProjM)
    Xn += np.transpose(np.dot(Xn,np.dot(ProjM,Xn))) 

    for i in range(len(X)):
        for j in range(len(X)):
            pairs = min([Xn[i,j],Xn[j,i]])
            Xn[i,j] -= pairs
            Xn[j,i] -= pairs
    return Xn,Xn-Xn



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
	DualityWeb = []
	TrialityMaps = []
	PermutationMaps = []
	Tuplemaps = TupFind(X)

	#%%
	def findpermutations(X1, F1, tmaps):
		relabels = [(X1,F1)]
		
		for tlist in tmaps:
			Xt,Ft = X1.copy(),F1.copy()
			for t in tlist:
				Xt, Ft = Swap(Xt, t), Swap(Ft, t)
			relabels += [(Xt,Ft)]
		return relabels
	#%%

	def inweb(xt,ft,pmaps,findn=False):
		for d,testphase in enumerate(pmaps):
			if np.sum(xt)==np.sum(testphase[0][0]) and np.sum(ft)==np.sum(testphase[0][1]):
				test = True
				if len(testphase) > 10000:
					X1,F1,X2,F2 = xt,ft, testphase[0][0],testphase[0][1]
					ANodes = [(sorted(X1[node]),sorted(np.transpose(X1)[node]),sorted((F1[node]))) for node in range(len(X1))]
					BNodes = [(sorted(X2[node]),sorted(np.transpose(X2)[node]),sorted((F2[node]))) for node in range(len(X1))]

					for node in ANodes:
						if node not in BNodes:
							test = False

				if test:
					for phase in testphase:
						if np.array_equal(xt,phase[0]) and np.array_equal(ft,phase[1]):
							if findn:
								return d
							return True
		return False

	qweb = [(X, F)]
	PermutationMaps += [findpermutations(X,F,Tuplemaps)]

	count = 0

	while len(qweb)!=0:
		phase = qweb[0]
		qweb = qweb[1:]

		Xt,Ft = phase[0].copy(),phase[1].copy()
		DualityWeb += [(Xt,Ft)]
		print('Dweb Length: '+str(len(DualityWeb)) + '| Qweb Length: '+str(len(qweb)),end="\r")
		
		for n in range(len(Xt)):
			if np.sum(np.transpose(Xt)[n]) + np.sum(Xt[n])==4:
				Xi,Fi = Sduality(Xt, Ft, n)
				count += 1
				TrialityMaps += [[(Xt,Ft),(Xi,Fi),n+1]]
				if not inweb(Xi,Fi,PermutationMaps) and not inweb(np.transpose(Xi),Fi,PermutationMaps):
					PermutationMaps += [findpermutations(Xi,Fi,Tuplemaps)]
					qweb += [(Xi,Fi)]
					

	TrialityTuples = [(inweb(t[0][0],t[0][1],PermutationMaps,findn=True),
		inweb(t[1][0],t[1][1],PermutationMaps,findn=True)) for t in TrialityMaps]
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

dweb = FindPhases(p4b,p4b-p4b)[0]

#%%