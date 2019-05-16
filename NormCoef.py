import os
import numpy as np
path = 'C:/Desktop/juicebox/Coef/' #path to folder contained contact matrix in juicebox verbose format (use -v to dump)
resolution = 5000 #resolution of contacts
#
#The name of contact matrices must have the next structure: SAMPLE_NAME.other_user_information
#
H = {}
files = os.listdir(path)
for file in files:
	parse = file.split('.') # the name of file with contacts must contain unique identificator delimited by points
	smpName = parse[0] # unique name of sample
	print(file)
	if smpName in H: pass
	else: H[smpName] = {}
	f = open(path+file,'r')
	lines = f.readlines()
	f.close()
	for i in range(len(lines)-1,-1,-1):
		parse = lines[i].split()
		del lines[i]
		try: pos = int(parse[0]), int(parse[1])
		except IndexError: pass
		else:
			c = float(parse[2])
			key = pos[0]/resolution
			if c != np.nan:
				if key in H[smpName]: H[smpName][key] += c
				else: H[smpName][key] = c
			else: pass
	del lines
for smpName in H:
	L = []
	Keys = list(H[smpName].keys())
	for key in Keys:
		L.append(H[smpName][key])
		del H[smpName][key]
	f = open(path+'coefficient.%s.txt' % smpName, 'w')
	print( 'cell_type\tcoeff\tmedian\tstd',file=f)
	print( '%s\t%.0f\t%.0f\t%.0f' % ( smpName, np.nanmean(L), np.nanmedian(L),np.nanstd(L) ), file=f )
	f.close()
