from nucleotides import *
import time
from primer import *
import primer3
import multiprocessing
from multiprocessing.managers import BaseManager
import json

tseq1 = 'CTCTTCTCTTCTCTTCTCTCTTCTCTTCTCTCTCTTCTCTCTTCTCTTCTTCTTCTTCTCTC'

# class Mymanager(BaseManager):
#     pass

# Mymanager.register('primersets',primersets)
# Mymanager.register('dict',dict)
# manager = Mymanager()

# if __name__ == "__main__":
#     manager.start()
#     dict2 = {'a':1,'b':2}
#     dict1 = manager.dict(dict2)
#     dict2['a'] = 3
#     print(dict1)

# solution = [1,2,3,4,5]
# ind = 4
# z = zip(solution,ind)

# existsprimer = {}
# existsquery = {'ACGT':12}

# fileprimer = open('fileprimer.txt','w')
# jsprimer = json.dumps(existsprimer)
# fileprimer.write(jsprimer)
# fileprimer.close()

# filequery = open('filequery.txt','w')
# jsquery = json.dumps(existsquery)
# filequery.write(jsquery)
# fileprimer.close()

class Mymanager(BaseManager):
    pass

Mymanager.register('primersets',primersets)
# Mymanager.register('dict',dict)
manager = Mymanager()



if __name__ == '__main__':

    # a = multiprocessing.Manager()
    # # a.start()
    # remotedict = a.dict()
    # remotedict['a'] = 1
    # # for key in remotedict:
    # #     print(remotedict[key])
    # for key in remotedict.keys():
    #     print(remotedict[key])
    manager.start()
    fasta = 'conseq.fasta'
    seqdict = readfasta(fasta)
    seq = seqdict['Consensus']
    a = manager.primersets(1,len(seq))
    b = manager.primersets(1,len(seq))
    a.getPrimerSeq(seq)
    b.getPrimerSeq(seq)
    print(a.getcode())
    print(b.getcode())
    # a,b = hybrid_primer(a,b)
    print(a.getcode())
    print(b.getcode())
