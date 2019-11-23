

def test_bigWig():
    from bigWigFileReader import bigWigReader
    from shared import Interval
    bwReader = bigWigReader("https://www.encodeproject.org/files/ENCFF966IHQ/@@download/ENCFF966IHQ.bigWig","Test")
    bwReader.readData(inMemory=False)
    res = bwReader.get_interval(Interval("chr1",1000,2000))
    print (res)

test_bigWig()