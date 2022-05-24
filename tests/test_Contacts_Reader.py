import sys
import os
head_folder_path = os.path.dirname(os.path.abspath(sys.argv[0]))+"/3Dpredictor"
source_path = os.path.dirname(os.path.abspath(sys.argv[0]))+"/3Dpredictor/source"
source_path2 = os.path.dirname(os.path.abspath(sys.argv[0]))+"/3Dpredictor/nn/source"
sys.path.append(source_path)
sys.path.append(source_path2)
from Contacts_reader import ContactsReader
from shared import Interval
import logging

interval = Interval("chr5", 74000000, 76400000)
binsize=5000
mindist = binsize*2+1
maxdist = 1500000
deletion = Interval("chr5", 75852814, 75881252)

contacts_reader = ContactsReader()
contacts_reader.generate_contacts_for_region(interval=interval, binsize=binsize, maxdist=maxdist, mindist=mindist)
print(contacts_reader.data["chr5"].keys())
print(contacts_reader.data["chr5"])
contacts_reader.delete_region(deletion)
print(contacts_reader.data["chr5"])
contacts = contacts_reader.get_contacts(Interval("chr5", 75200000, 76400000),mindist=mindist,maxdist=maxdist)
print(contacts)
