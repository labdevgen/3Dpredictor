**How to generate data using new nn modules**

I.) where to get data

1.) download genomic sequence (letters) in fasta format for hg38 genome:

http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

unzip is not necessary, but recommended for faster loading

2.) download ChipSeq data in binary format. 
There are many of them. For the beginning you can use following:

https://www.encodeproject.org/files/ENCFF473IZV/@@download/ENCFF473IZV.bigWig --> CTCF

https://www.encodeproject.org/files/ENCFF986PCY/@@download/ENCFF986PCY.bigWig --> H3K27

https://www.encodeproject.org/files/ENCFF795YDZ/@@download/ENCFF795YDZ.bigWig --> polR2a

3.) download binary hic-file for H1-ES cells:

https://data.4dnucleome.org/files-processed/4DNFI2TK7L2F/

II)

0.) requered imports:

'''python

    from bigWigFileReader import bigWigReader
    from fastaFileReader import fastaReader
    from hicFileReader import hicReader
    from shared import Interval, Genome

'''

1.) to generate genome sequence:

'''python

    faReader = fastaReader(path="../input/hg38/hg38.fa",useOnlyChromosomes=["chr1"])
    # path - path to fasta file(s) with sequence
    # useOnlyChromosomes - a list of chromsomes to use. Other chromosomes in file(s) will be ignored
    # excludeChromosomes - a list of chromsomes to exclude. These chromosomes will not be load from file

    res = faReader.get_interval(Interval("chr1",start,end))
    # this will extract fragment of seq from chrosome with the name chr1 from start to end (excluding end)
    # result will be a numpy array with numbers [0,1,2,3,4] corresponding to DNA letters (A,T,G,C,N)

    sizes = faReader.get_chr_sizes()
    # sizes = dict{"chr1":len_of_chr_1, "chr2":len_of_chr_2 ... }
'''

2.) to generate ChipSeq-data, storing whole dataset in-memory:

'''python

    path = "../input/ENCFF966IHQ.bigWig"
    bwReader = bigWigReader(path,genome=faReader)
    bwReader.readData(inMemory=True)

    # path - name of file with row ChipSeq-data
    # faReader - fastaReader object (see above)
    # note that this will load only data for those chromosomes which are loaded in faReader

    result = bwReader.get_interval(Interval("chr1",start,end))
    # result will be a numpy array with floats corresponding to ChipSeq signals for each letter in the interval

3.) to generate ChipSeq-data, reading data each time from file:

Same as above, but use:

'''python

    bwReader.readData(inMemory=False)

'''

4.) to load Hi-C contacts:

'''python

    hic = hicReader(path, genome=faReader, resolution = 5000)
    # path and genome - same as above
    # resolution - size of the bin of Hi-C matrix.
    # I recommend to use 100000 for tests (low memory, but not much biological scense), 1000 or 5000 for real trainings
    # Again, only data for chromsomes in genome object will be loaded
    hic.read_data()
    result = hic.get_contact(Interval("chr1",start,end)) # single float value, NaN or None
                                                         # NaN means contact was not defined (NaN in original data)
                                                         # None probably means contact was equal to 0 (absent in original sparce matrix)
    result = hic.get_chr_contact("chr1") # returns sparse matrix of the whole chrm as pandas dataframe
    # fields: ["st", "en", "count"]
'''