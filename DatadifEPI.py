import numpy as np
import matplotlib.pyplot as plt
import func
reload(func)

base = 'C:/Desktop/Arbeit/AllData/Genome/'
prm = '/RefGene/mm10.refGene.txt' # pass to promoters
#prm = '/RefGene/hg19.refGene.txt'
enh = '/Enhancer/mm/mm10.enhancer.total.bed' # pass to enhancers
# enh = '/Enhancer/hg/hg19.enhancer.total.bed'
bin = 5000
contacts = 'contacts/Hepat/hepat.control','/contacts/NPC/npc.control','contacts/Hepat/hepat.rep1','contacts/Hepat/hepat.rep2', # pass to control oe dumped by juicer
#contacts = '/contacts/Mcyte/Mcyte.chr1.contacts','K562'
models = 'contacts/Hepat/hepat.data','contacts/NPC/npc.data', # pass to model oe dumped by juicer
# contacts = 'contacts/Mcyte/Mcyte.control', 'contacts/K562/K562.control','contacts/K562/K562.rep1.','contacts/K562/K562.rep2.'
# models = 'contacts/Mcyte/Mcyte.data','contacts/K562/K562.data'
suf = 'oe.none','oe.kr' # extention of files
rng = 1500000/bin

p = func.readFeatures(base+prm, bin, feature='promoter') # read promoters
e1 = func.readFeatures(base+enh, bin, feature='enhancer') # read enhancers

EPP = func.EPpairer(p,e1,rng) # pairing promoters and enhancers

c1 = {}
c2 = {}
c1m = {}
c2m = {}
c1r1 = {}
c1r2 = {}
c = range(1,23)
c.append('X') # chromosome list
print
# reading OE for EPP
for i in c:
	chrm = 'chr' + str(i)
	print chrm
	try:
		c1.update( func.readFeatures( ('%s/%s.%s.%i.%s') % (base,contacts[0],chrm,bin,suf[0]), bin, feature='real',chr=chrm,range=rng,filter=EPP) )
		c2.update( func.readFeatures( ('%s/%s.%s.%i.%s') % (base,contacts[1],chrm,bin,suf[0]), bin, feature='real',chr=chrm,range=rng,filter=EPP) )
		c1m.update( func.readFeatures( ('%s/%s.%s.%i.%s') % (base,models[0],chrm,bin,suf[0]), bin, feature='real',chr=chrm,range=rng,filter=EPP) )
		c2m.update( func.readFeatures( ('%s/%s.%s.%i.%s') % (base,models[1],chrm,bin,suf[0]), bin, feature='real',chr=chrm,range=rng,filter=EPP) )
		c1r1.update( func.readFeatures( ('%s/%s.%s.%i.%s') % (base,contacts[2],str(i),bin,suf[1]), bin, feature='real',chr=chrm,range=rng,filter=EPP) )
		c1r2.update( func.readFeatures( ('%s/%s.%s.%i.%s') % (base,contacts[3],str(i),bin,suf[1]), bin, feature='real',chr=chrm,range=rng,filter=EPP) )
	except IOError: pass
