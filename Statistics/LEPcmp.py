import numpy as np
import matplotlib.pyplot as plt
import LEPfunction as lf
reload(lf)
#
#Bad writed scripts for a estimation of coincidences of loops and enhancer-promoter interactions
#
base_dir = 'C:/Desktop/Arbeit/AllData/'
cell_type = [ 'Mcyte', 'Mphage']#, 'Tcell']
enh_path = base_dir+'Enhancer/hg19.'
suf = '.enhancer.tissue.bed'
out = '.dist'
loop_dir = 'C:/Desktop/juicebox/'
tr_loop = '.loops/merged_loops'
resolution = 25000

for ct in cell_type:
	X = []
	Y = np.zeros(10)
	rY = np.zeros((100,10))
	D = np.zeros(10)
	M = np.zeros(10)
	L = lf.readLoops(loop_dir+ct+tr_loop,resolution,'chr')
	E = lf.readLoops(enh_path+ct+suf,resolution,'')
	dL = set([])
	dE = set([])
	print 'cell type:', ct
	print 'num loops: %i, num enhancers %i' % (len(L),len(E))
	
	sE = {}
	dE = {}
	for e in E: 
		if sE.has_key(e[0]):
			if sE[e[0]][0] > e[1]: sE[e[0]][0] = e[1]
			elif sE[e[0]][1] < e[3]: sE[e[0]][1] = e[3]
			else: pass
		else:  sE[e[0]] = [min(e[1],e[3]),max(e[1],e[3])]
		if dE.has_key(e[0]): dE[e[0]].append(abs(e[3]-e[1]))
		else: dE[e[0]] = [ abs(e[3]-e[1]) ]
	for e in E:
		h = []
		for l in L:
			if e[0] == l[0] and abs(e[1] - l[1]) < 10 and e[2] == l[2] and abs(e[3]-l[3]) < 10:
				i = max(abs(e[1] - l[1]),abs(e[3] - l[3]))
				h.append(i)
			elif e[0] == l[0] and abs(e[3] - l[1]) < 10 and e[2] == l[2] and abs(e[1]-l[3]) < 10:
				i = max(abs(e[3] - l[1]),abs(e[1] - l[3]))
				h.append(i)
			else: pass
		if len(h) > 0: 
			Y[min(h)] += 1
			#X.append((e[0],e[1]*resolution,e[1]*resolution+resolution,min(h)))
			#X.append((e[0],e[3]*resolution,e[3]*resolution+resolution,min(h)))
	#X.sort()
	for r in range(100):
		rE = []
		for key in sE:
			sS = np.random.randint(sE[key][0],high=sE[key][1]+1, size=len(dE[key]))
			for i in range(len(sS)): rE.append((key,sS[i],key,sS[i]+dE[key][i]))
		for e in rE:
			h = []
			for l in L:
				if e[0] == l[0] and abs(e[1] - l[1]) < 10 and e[2] == l[2] and abs(e[3]-l[3]) < 10:
					i = max(abs(e[1] - l[1]),abs(e[3] - l[3]))
					h.append(i)
				elif e[0] == l[0] and abs(e[3] - l[1]) < 10 and e[2] == l[2] and abs(e[1]-l[3]) < 10:
					i = max(abs(e[3] - l[1]),abs(e[1] - l[3]))
					h.append(i)
				else: pass
			if len(h) > 0: 
				rY[r][min(h)] += 1
	for j in range(1,10): 
		Y[j] += Y[j-1]
		rY[:,j] += rY[:,j-1]
	for i in range(100): M += rY[i]/100
	for i in range(100): D += (rY[i]-M)**2/99
	plt.cla()
	plt.plot(range(10),Y,'r-',lw=3)
	plt.plot(range(10),M, c = '0.5',lw=3)
	plt.fill_between(range(10),M-D,M+D, color = '0.9', facecolor='0.9')
	plt.savefig(enh_path+ct+out + '.png', dpi=400)
	plt.clf()
	#file = open(enh_path+ct+out+'.bedGraph','w')
	#for i in X: print >> file, '%s\t%i\t%i\t%i' % (i[0],i[1],i[2],i[3])
	#file.close()
	print 'observed:', Y
	print 'mean:', M
	print 'disp:', D