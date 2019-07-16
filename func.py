import numpy as np

def readFeatures(path, resolution, **kwargs):
	f = open(path,'r')
	lines = f.readlines()
	f.close()
	try: filter = kwargs['filter']
	except KeyError: filter = False
	if kwargs['feature'] == 'promoter':
		ftr = {}
		for line in lines:
			parse = line.split()
			if parse[3] == '+': 
				try: ftr[parse[2]].add( int(parse[4])/resolution )
				except KeyError: ftr[parse[2]] = set([ int(parse[4])/resolution ])
			else:
				try: ftr[parse[2]].add( int(parse[5])/resolution )
				except KeyError: ftr[parse[2]] = set([ int(parse[5])/resolution ])
	elif kwargs['feature'] == 'enhancer':
		ftr = {}
		for line in lines:
			parse = line.split()
			try: ftr[parse[0]].add( (int(parse[1])+int(parse[2]))/(2*resolution) )
			except KeyError: ftr[parse[0]] = set([ (int(parse[1])+int(parse[2]))/(2*resolution) ])
	elif kwargs['feature'] == 'real':
		ftr = {}
		chr = kwargs['chr']
		try: rng = kwargs['range']
		except KeyError: rng = 100
		for i in range(len(lines)-1,-1,-1):
			parse = lines[i].split()
			coor = [ int(parse[0])/resolution, int(parse[1])/resolution ]
			coor.sort()
			if abs(coor[1] - coor[0]) < rng: 
					key = (chr, coor[0], chr, coor[1])
					try: 
						filter[key]
						ftr[key] = float(parse[2])
					except KeyError: pass
			del lines[i]
	elif kwargs['feature'] == 'model':
		ftr = {},{}
		chr = kwargs['chr']
		try: rng = kwargs['range']
		except KeyError: rng = 100
		for line in lines:
			parse = line.split()
			coor = [ int(parse[0])/resolution, int(parse[1])/resolution ]
			coor.sort()
			if abs(coor[1] - coor[0]) < rng: ftr[0][ (chr, coor[0], chr, coor[1]) ] = float(parse[2])
			if abs(coor[1] - coor[0]) < rng: ftr[1][ (chr, coor[0], chr, coor[1]) ] = float(parse[3])
	elif kwargs['feature'] == 'pre':
		ftr = {}
		for line in lines:
			parse = line.split()
			ftr[ ( parse[1], int(parse[2])/resolution, parse[5], int(parse[6])/resolution ) ] = float(parse[8]) 
	else: pass
	del lines
	return ftr

def metric(Contacts1,Contacts2,region):
	M = 0
	k = 0
	for key in region:
		try:
			c1 = Contacts1[key]
			c2 = Contacts2[key]
			if c1 > 0 and c2 > 0:
				M += abs( np.log2(c1) - np.log2(c2) )
				k += 1
			else: pass
		except KeyError: pass
	try: M = M/k
	except ZeroDivisionError: M = 0
	return M

def difFinder(Contacts1, Contacts2, resolution, out, **kwargs):
	try: method = kwargs['method']
	except KeyError: return "method!"
	try: frame = kwargs['frame']
	except KeyError: frame = 10
	try: high = kwargs['high']
	except KeyError: frame = 10
	Keys = set([])
	Metric = []
	n=0
	for key in Contacts1:
		Keys.add(key[:2])
		Keys.add(key[2:])
	for key in Contacts2:
		Keys.add(key[:2])
		Keys.add(key[2:])
	Keys=list(Keys)
	if method == 'triangle':
		for key in Keys:
			region = []
			for i in range(key[1]-frame,key[1]+frame+1):
				for j in range(i,key[1]+frame+1): region.append((key[0],i,key[0],j))
			Metric.append( ( key[0],key[1]*resolution,key[1]*resolution+resolution,metric(Contacts1,Contacts2,region) ) )
			if n == 50: print region
			n += 1
	elif method == 'strip':
		for key in Keys:
			region = []
			for i in range(key[1]-frame,key[1]+1):
				for j in range(i-high,i+1 ): region.append((key[0],i,key[0],j))
			for i in range(key[1],key[1]+frame+1):
				for j in range(i,i+high+1 ): region.append((key[0],i,key[0],j))
			Metric.append( ( key[0],key[1]*resolution,key[1]*resolution+resolution,metric(Contacts1,Contacts2,region) ) )
			if n == 50: print region
			n += 1
	else: pass
	Metric.sort()
	f = open(out,'w')
	for m in Metric: print >> f, '%s\t%i\t%i\t%.5f' % m
	f.close()
	return Metric

def EPpairer(prm,enh,rng):
	pairs = {}
	Keys = list(set(prm.keys()) & set(enh.keys()))
	Keys.sort()
	for key in Keys:
		if len(key) > 5: pass
		else:
			print key,
			for i in prm[key]:
				for j in enh[key]:
					d = abs(i-j)
					if d<rng and d>2:
						if i < j: pairs[(key,i,key,j)] = 0
						else: pairs[(key,j,key,i)] = 0
	return pairs

def OEdumping(epp,*args):
	c = []
	for i in epp:
		cs = []
		for j in args:
			try: cs.append( j[i])
			except KeyError: cs.append( 0 )
		cs.append( abs(i[3]-i[1]) )
		c.append( cs )
		#print len(c),len(cs)
	print len(cs)
	c = np.array(c)
	c = c[((c[:,0] != 0) & (c[:,1] != 0) & (c[:,2] != 0) & (c[:,3] != 0) & (c[:,4] != 0) & (c[:,5] != 0)),:]
	c = np.log2(c)
	return c
	