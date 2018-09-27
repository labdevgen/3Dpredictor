import LEPfunction as lf
reload(lf)
#Filtring enhancer by tissue or cell specificity and
#finding enhancer-promoter interaction based on enhancer-gene dependences
base_dir = 'C:/Desktop/Arbeit/AllData/'
gene_path = base_dir+'Genome/hg19.refGene.txt'
enh_path = base_dir+'Enhancer/hg19'
ts_enh = '.Mphage','.Mcyte','.Tcell'
suf = '.enhancer'

resolution = 2500
resolution2 = 25000
refGene = lf.readGene(gene_path,resolution)
refEn = lf.readEnhancer(enh_path+suf,resolution,0.5)
for i in range(3):
	tsEn = lf.readTissue(enh_path+ts_enh[i]+suf+'.bed',resolution)
	print len(tsEn)
	fltEn = lf.filtringTissue(refEn,tsEn)
	print len(fltEn)
	
	PEI = lf.interactions(refGene,fltEn,resolution,resolution2,5)
	file = open(enh_path+ts_enh[i]+suf+'.tissue.bed','w')
	print >> file, 'chr1\tx1\tx2\tchr2\ty1\ty2'
	for p in PEI: print >> file, '%s\t%i\t%i\t%s\t%i\t%i' % ( p[0],p[1]*resolution2,p[1]*resolution2+resolution2,p[2],p[3]*resolution2,p[3]*resolution2+resolution2 )
	file.close()