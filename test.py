from nucleotides import *
import time
from primer import *
import primer3

tseq1 = 'CTCTTCTCTTCTCTTCTCTCTTCTCTTCTCTCTCTTCTCTCTTCTCTTCTTCTTCTTCTCTC'

# seqinfo(tseq1)
# print(type(100*conCG(tseq1)))

if __name__ == "__main__":

    time1 = time.time()

    fasta = 'conseq.fasta'
    seqdict = readfasta(fasta)
    seq = seqdict['Consensus']

    a = primersets(1,len(seq))
    a.getPrimerSeq(seq)
    a.fitness()
    print(a.blastscore)
    print(a)
    time2 = time.time()
    print(time2-time1)