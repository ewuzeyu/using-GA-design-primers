#-*- coding = utf-8
import os
import math

# 计算Tm值
def ptm(query):
    SanHan = {'AA': -8.4, 'AT': -6.5, 'AC': -8.6, 'AG': -6.1, 'TA': -6.3, 'TC': -7.7, 'TG': -7.4, 'TT': -8.4, 'CG': -10.1, 'CA': -7.4, 'CT': -6.1, 'CC': -6.7, 'GA': -7.7, 'GT': -8.6, 'GC': -11.1, 'GG': -6.7}
    SanSang = {'AA': -23.6, 'AT': -18.8, 'AC': -23, 'AG': -16.1, 'TA': -18.5, 'TC': -20.3, 'TG': -19.3, 'TT': -23.6, 'CG': -25.5, 'CA': -19.3, 'CT': -16.1, 'CC': -15.6, 'GA': -20.3, 'GT': -23, 'GC': -28.4, 'GG': -15.6}
    kcl = 50
    mg = 2
    dntp = 0.2
    primer = 0.4
    adjusted = 16.6*math.log10((kcl+120*math.sqrt(mg-4*dntp))/1000)
    San_Han = 0
    San_Sang = 0
    for i in range(len(query)-1):
        NN = query[i:i+2]
        San_Han += SanHan[NN]
        San_Sang += SanSang[NN]

    if query.count('C') + query.count('G') > 0:
        San_Sang += -5.9
    else:
        San_Sang += 0.6
    
    San = San_Han *1000/(San_Sang +1.987*math.log(primer*10**-6))-273.15+adjusted-3
    return San

# 计算blast得分
def blastn(query):
    # query = r'GGACAAACGTCATAACTAGC'
    with open(r'C:\Users\admin\Desktop\zhuanhuan\Blast\mytemp.fasta','w') as myfa:
        myfa.write('>myquery\n')
        myfa.write(query)

    # blastcmd = 'blastn -query C:/Users/admin/Desktop/zhuanhuan/Blast/mytemp.fasta -db C:/Users/admin/Desktop/zhuanhuan/Blast/FCVfull -evalue 1 -task blastn -outfmt "6 delim=  qlen qstart qend mismatch"'
    # blastcmd = 'blastn -query C:/Users/admin/Desktop/zhuanhuan/Blast/mytemp.fasta -db C:/Users/admin/Desktop/zhuanhuan/Blast/FCVfull -evalue 1 -task blastn -outfmt "6 delim=  sstart send qstart qend qlen saccver"'
    sum = 0

    with os.popen(r'blastn -query C:/Users/admin/Desktop/zhuanhuan/Blast/mytemp.fasta -db C:/Users/admin/Desktop/zhuanhuan/Blast/FCVfull -evalue 1 -task blastn -outfmt "6 delim=  qlen qstart qend mismatch"') as blast:
        r = blast.readlines()
        # with open('get.txt','w') as fasta:
        #     for line in r:
        #         l = line.split()
        #         sstart = int(l[0])
        #         send = int(l[1])
        #         qstart = int(l[2])
        #         qend = int(l[3])
        #         qlen = int(l[4])
        #         sid = l[5]
        #         if sstart < send:
        #             small = sstart-(qstart-1)
        #             large = send+(qlen-qend)
        #         else:
        #             small = send-(qlen-qend)
        #             large = sstart+(qstart-1)
        #         fasta.write(sid+'\t'+str(small)+'\t'+str(large)+'\n')
        for line in r:
            l = line.split()
            qlen = int(l[0])
            qstart = int(l[1])
            qend = int(l[2])
            mism = int(l[3])
            score = 1 - 0.5*(qstart - 1)/qlen - 1*(qlen - qend)/qlen - 0.3*mism/qlen
            # print(score)
            sum = sum + score


    print(query+':'+str(sum)+'  Tm:'+str(ptm(query)))

t1 = ['A','C','T']
t2 = ['T','C']
t3 = ['A','G']

for i in range(3):
    for j in range(2):
        for k in range(2):
            blastn('CG'+t1[i]+'GCAAG'+t2[j]+'GTTA'+t3[k]+'CTTGAC')
            # print(ptm('CG'+t1[i]+'GCAAG'+t2[j]+'GTTA'+t3[k]+'CTTGAC'))

