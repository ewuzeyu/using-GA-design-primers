from nucleotides import *
import time
from primer import *
import primer3

tseq1 = 'CTCTTCTCTTCTCTTCTCTCTTCTCTTCTCTCTCTTCTCTCTTCTCTTCTTCTTCTTCTCTC'

# seqinfo(tseq1)
# print(type(100*conCG(tseq1)))

time1 = time.time()
# a = readfasta(r'C:\Users\admin\Downloads\phageMS2.fasta')
# seqt2 = a['NC_001417.2 phage MS2 genome']

seqinfo(tseq1)
print(primer3.calcTm(tseq1))
# a = []
# for i in range(100000):
#     a.append(primersets())
# for i in range(100):
#     a = primersets()
# print(a.primers)
#     print(degbase_generator(a.primers,6))

# l1 = [1,2,3,4]
# l2 = [2,2,2,2]
# print(multiply(l1,l2))
# with open('log.txt',mode='a') as log:
# fasta = 'conseq.fasta'
# seqdict = readfasta(fasta)
# seq = seqdict['Consensus']
#     for i in range(10):
#         a = primersets(1,len(seq))
#         a.getPrimerSeq(seq)
#         if a.Filter() == 1:
#             log.write(a.__str__()+'\n'+str(a.fitness())+'\n')
#             time3 = time.time()
#             print(time3-time1)
        # for key in a.primer_seq:
        #     log.write(str(a.primer_seq[key].count('C')+a.primer_seq[key].count('G')) + '\n')
# pr = a.primer_seq
# for key in pr:
#     seqinfo(pr[key],key)
#     print('========--------****--------========')
# print(a)

# a = primersets()
# b = primersets()

# print(a.init_primers)
# print(b.init_primers)
# hybrid_primer(a,b)

# print(a.init_primers)
# print(b.init_primers)

# a = primersets(1,7000)
# a.printself()
# b = duplicate_primer(a)

# print(a.mutate())

# a.printself()
# b.printself()


# fasta = 'conseq.fasta'
# seqdict = readfasta(fasta)
# seq = seqdict['Consensus']

# myprimer = {'start':107,'p1':21,'p2':1,'p3':21,'p4':24,'p5':20,'gap':30,'p10':20,'p9':20,'p8':20,'p7':8,'p6':18}

# my = primersets(existprimer=myprimer)
# my.getPrimerSeq(seq)
# my.fitness()
# print(my)
# solution = []
# filtered = []
# ranked = []

# for i in range(5):
#     t = 0
#     while t != 1:
#         a = primersets(1,len(seq))
#         a.getPrimerSeq(seq)
#         t = a.Filter()
#     solution.append(a)
#     # print(a.primer_seq)


# for i in range(5):
#     a = random.randint(0,4)
#     b = a
#     while b == a:
#         b = random.randint(0,4)
#     p1 = duplicate_primer(solution[a])
#     p2 = duplicate_primer(solution[b])
#     hybrid_primer(p1,p2)
#     p1.getPrimerSeq(seq)
#     filtered.append(p1)

# for r in solution:
#     if random.random()<0.5:
#         m = duplicate_primer(r)
#         m.mutate()
#         m.getPrimerSeq(seq)
#         ranked.append(m)
# print('yuanshi')
# for t1 in solution:
#     print(t1)
# print('zajiao')
# for t2 in filtered:
#     print(t2)
# print('bianyi')
# for t3 in ranked:
#     print(t3)




    # print(ps.fitness())

# print('----')
# ranked = quickSort(solution)
# print(filtered)
# for k in range(len(filtered)):
#     print(filtered[k])
# for c in ranked:
#     print(c)


time2 = time.time()
print(time2-time1)