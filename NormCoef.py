import os
import numpy as np
base = 'C:/Desktop/juicebox/' #base folder contains subfolder with contacts in verbose format
path = 'Coef/',#'Dactyly/',#GM12878/replicate1/','GM12878/replicate2/' - list of subfolders
resolution = 5000 #resolution of contacts
#
#
#
H = {}
for p in path:
	files = os.listdir(base+p)
	for file in files:
		parse = file.split('.') # the name of file with contacts must contain unique identificator delimited by points
		smpName = parse[2] # idetificator
		res = int(parse[3][:-2])*1000
		print file
		if H.has_key(smpName) == True: pass
		else: H[smpName] = {}
		f = open(base+p+file,'r')
		lines = f.readlines()
		f.close()
		for i in range(len(lines)-1,-1,-1):
			parse = lines[i].split()
			del lines[i]
			try: pos = int(parse[0]), int(parse[1])
			except IndexError: pass
			else:
				c = float(parse[2])
				key = pos[0]/res
				if c != np.nan:
					if H[smpName].has_key(key) == True: H[smpName][key] += c
					else: H[smpName][key] = c
				else: pass
				# Uncomment next block if you want to calculate the TRUE normed coefficient.
				# key = pos[1]/resolution
				# if c != np.nan:
					# if H[smpName].has_key(key) == True: H[smpName][key] += c
					# else: H[smpName][key] = c
				# else: pass
			del lines

	for smpName in H:
		L = []
		for key in H[smpName].keys():
			L.append(H[smpName][key])
			del H[smpName][key]
		f = open(base+'coefficient.%s.txt' % smpName, 'w')
		print >> f, 'cell_type\tcoeff\tmedian\tstd'
		print >> f, '%s\t%.0f\t%.0f\t%.0f' % ( smpName, np.nanmean(L), np.nanmedian(L),np.nanstd(L) )
		f.close()
