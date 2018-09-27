#module with base functions

def readGene(fname,resolution):
	refGene = {}
	file = open(fname, 'r')
	lines = file.readlines()
	file.close()
	for line in lines:
		parse = line.split('\t')
		if parse[3] == '+': refGene[parse[12]] = parse[2],int(parse[4])/resolution
		else: refGene[parse[12]] = parse[2],int(parse[5])/resolution
	return refGene
	
def readEnhancer(fname,resolution,score):
	refEn = {}
	file = open(fname, 'r')
	lines = file.readlines()
	file.close()
	for line in lines:
		parse = line.split('\t')
		if parse[0] == 'Enhancer': 
			key = parse[1],(int(parse[2])+int(parse[3]))/(2*resolution)
			refEn[key] = []
		else:
			if parse[2] >=score: refEn[key].append(parse[1])
	return refEn

def readTissue(fname,resolution):
	tsEn = {}
	file = open(fname, 'r')
	lines = file.readlines()
	file.close()
	for line in lines:
		parse = line.split('\t')
		key = parse[0],(int(parse[1])+int(parse[2]))/(2*resolution)
		tsEn[key] = 1
	return tsEn

def readLoops(fname,resolution,pref):
	tsLoop = []
	file = open(fname, 'r')
	lines = file.readlines()
	file.close()
	for line in lines[1:]:
		parse = line.split('\t')
		tsLoop.append( (pref+parse[0],(int(parse[1])+int(parse[2]))/(2*resolution),pref+parse[3],(int(parse[4])+int(parse[5]))/(2*resolution) ) )
	tsLoop.sort()
	return tsLoop

def filtringTissue(refEn,tsEn):
	fltEn = {}
	for key in tsEn:
		try: fltEn[key] = refEn[key]
		except KeyError: pass
	return fltEn

def interactions(refGene,refEn,res1,res2,frame):
	PromEnInter = []
	for key in refEn:
		#print refEn[key]
		for i in refEn[key]:
			try: concor = refGene[i][0],refGene[i][1]*res1/res2,key[0],key[1]*res1/res2
			except KeyError: pass
			else:
				if abs(concor[1] - concor[3]) > frame and len(concor[0]) < 6 and len(concor[2]) < 6: PromEnInter.append( concor )
				else: pass
	PromEnInter.sort()
	return PromEnInter