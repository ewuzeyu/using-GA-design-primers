#-*- coding = utf-8
# 引物处理

import random
from nucleotides import Tm,ifdimer,conCG
import os
import multiprocessing
import threading

PRIMER_LENGTH = (18,22)
PRIMER_23GAP = (0,60)
PRIMER_CGAP = (120,160)
PRIMER_12CLEN = (40,60)

def dictdeepcopy(d:dict) -> dict:
    a = {}
    for key in d:
        a[key] = d[key]
    return a

def listdeepcopy(l:list) -> list:
    a = []
    for i in range(len(l)):
        a.append(l[i])
    return a

def randomlist(start:str, stop:str, length:str) -> list:
    '''生成随机从start到stop共length长的数组'''
    start, stop = (int(start), int(stop)) if start <= stop else (int(stop), int(start))
    length = int(abs(length)) if length else 0
    random_list = []
    for i in range(length):
        random_list.append(random.randint(start, stop))
    return random_list

def dec2vigse(n, x = 26):
    #n为待转换的十进制数，x为机制，取值为2-16
    a=['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
    b=[]
    while True:
        s=n // x  # 商
        y=n % x  # 余数
        b=b+[y]
        if s==0:
            break
        n=s
    b.reverse() # 辗转相除法
    ans = ''
    for i in b:
        ans += a[i]
    return ans

def vigse2dec(n, x = 26):
    a = {'a':0,'b':1,'c':2,'d':3,'e':4,'f':5,'g':6,'h':7,'i':8,'j':9,'k':10,'l':11,'m':12,'n':13,'o':14,'p':15,'q':16,'r':17,'s':18,'t':19,'u':20,'v':21,'w':22,'x':23,'y':24,'z':25}
    result = 0
    for i in range(len(n)):
        result *= x
        result += a[n[i]]
    return result

def primer2zip(p:dict) -> str:
    res = ''
    start = dec2vigse(p['start'])
    start = 'a'*(3-len(start)) + start
    res = res + start
    gap = dec2vigse(p['gap'])
    gap = 'a'*(2-len(gap)) + gap
    res = res + gap
    a1 = ['p1','p3','p5','p6','p8','p10']
    a2 = ['p2','p4','p7','p9']
    for t in a1:
        res = res + dec2vigse(p[t])
    for t in a2:
        pd = dec2vigse(p[t])
        res = res + 'a'*(2-len(pd)) + pd
    return res

def zip2primer(s:str) -> dict:
    a = {}
    b1 = ['p1','p3','p5','p6','p8','p10']
    b2 = ['p2','p4','p7','p9']
    a['start'] = vigse2dec(s[0:3])
    a['gap'] = vigse2dec(s[3:5])
    for i in range(len(b1)):
        a[b1[i]] = vigse2dec(s[i+5:i+6])
    for i in range(len(b2)):
        a[b2[i]] = vigse2dec(s[i*2+11:i*2+13])
    return a

def multiply(l1:list,l2:list) -> list:
    if len(l1) != len(l2):
        return 0
    else:
        res = []
        for (a,b) in zip(l1,l2):
            res.append(a*b)
        return res
    
def primer_generator(start:int,end:int) -> dict:
    '''
    生成一随机引物组'''
    primer = {}
    primer['p1'] = random.randint(PRIMER_LENGTH[0],PRIMER_LENGTH[1])
    primer['p3'] = random.randint(PRIMER_LENGTH[0],PRIMER_LENGTH[1])
    primer['p5'] = random.randint(PRIMER_LENGTH[0],PRIMER_LENGTH[1])
    primer['p6'] = random.randint(PRIMER_LENGTH[0],PRIMER_LENGTH[1])
    primer['p8'] = random.randint(PRIMER_LENGTH[0],PRIMER_LENGTH[1])
    primer['p10'] = random.randint(PRIMER_LENGTH[0],PRIMER_LENGTH[1])
    primer['p2'] = random.randint(PRIMER_23GAP[0],PRIMER_23GAP[1])
    primer['p7'] = random.randint(PRIMER_23GAP[0],PRIMER_23GAP[1])
    primer['p4'] = random.randint(PRIMER_12CLEN[0]-primer['p3'],PRIMER_12CLEN[1]-primer['p3'])
    primer['p9'] = random.randint(PRIMER_12CLEN[0]-primer['p8'],PRIMER_12CLEN[1]-primer['p8'])
    primer['gap'] = random.randint(PRIMER_CGAP[0]-primer['p4']-primer['p5']-primer['p9']-primer['p10'],
                                PRIMER_CGAP[1]-primer['p4']-primer['p5']-primer['p9']-primer['p10'])
    if primer['gap'] < 0:
        primer['gap'] = 0
    primer['start'] = random.randint(start, end-primer['p3']-primer['p4']-primer['p5']-primer['p8']-primer['p9']-primer['p10']-primer['p1']-primer['p2']-primer['p6']-primer['p7']-primer['gap'])
    return primer

def degbase_generator(primer:dict,number:int = 6) -> list:
    totallen = primer['p1'] + primer['p3'] + primer['p5'] + primer['p6'] + primer['p8'] + primer['p10']
    a = randomlist(1,totallen,number)
    mask = randomlist(0,1,number)
    prere = multiply(a,mask)
    return prere

def readfasta(filepath:str) -> dict:
    with open(file=filepath,mode='r') as f:
        a = {}
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                name = line[1:]
                a[name] = ''
            else:
                a[name] += line
    return a

# 计算blast得分
def blastn(query:str, q) -> float:
    # return random.randint(100,600)
    # query = r'GGACAAACGTCATAACTAGC'
    # print('-',end='')
    with open('C:\\Users\\admin\\Desktop\\zhuanhuan\\Blast\\qprimer\\'+query+'.fasta','w') as myfa:
        myfa.write('>myquery\n')
        myfa.write(query)

    # blastcmd = 'blastn -query C:/Users/admin/Desktop/zhuanhuan/Blast/mytemp.fasta -db C:/Users/admin/Desktop/zhuanhuan/Blast/FCVfull -evalue 1 -task blastn -outfmt "6 delim=  qlen qstart qend mismatch"'
    # blastcmd = 'blastn -query C:/Users/admin/Desktop/zhuanhuan/Blast/mytemp.fasta -db C:/Users/admin/Desktop/zhuanhuan/Blast/FCVfull -evalue 1 -task blastn -outfmt "6 delim=  sstart send qstart qend qlen saccver"'
    sum = 0

    with os.popen('blastn -query C:/Users/admin/Desktop/zhuanhuan/Blast/qprimer/'+query+'.fasta -db C:/Users/admin/Desktop/zhuanhuan/Blast/FCVfull -evalue 1 -task blastn -outfmt "6 delim=  qlen qstart qend mismatch"') as blast:
        r = blast.readlines()
        # print(r)
        for line in r:
            l = line.split()
            qlen = int(l[0])
            qstart = int(l[1])
            qend = int(l[2])
            mism = int(l[3])
            score = 1 - 0.5*(qstart - 1)/qlen - 1*(qlen - qend)/qlen - 0.3*mism/qlen
            # print(score)
            sum = sum + score
    # return sum
    # q.send(sum)
    q.put(sum)

class primersets:
    '''
    引物组'''
    
    

    def __init__(self,start:int =1 , end:int =7000, degbase:int = 6, existprimer:dict = None, existdegbase:list = None) -> None:
        '''初始化引物，随机生成'''
        self.init_primers = {}
        self.sstart = 1
        self.send = 7000
        if existprimer:
            self.init_primers = existprimer
        else:
            self.sstart = start
            self.send = end
            self.init_primers = primer_generator(start,end)
        self.code = primer2zip(self.init_primers)
        self.degbase = []
        if existdegbase:
            self.degbase = existdegbase
        else:
            self.degbase = degbase_generator(self.init_primers,degbase)
        return

    def getPrimerSeq(self,toSeq:list) -> None:
        self.order_primers = {}
        self.primer_seq = {}
        order = ['start','p1','p2','p3','p4','p5','gap','p10','p9','p8','p7','p6']
        lasti = ''
        for i in order:
            if i == 'start':
                self.order_primers[i] = self.init_primers[i]
            else:
                self.order_primers[i] = self.order_primers[lasti] + self.init_primers[i]
            lasti = i

        self.primer_seq['F3'] = toSeq[self.order_primers['start']-1:self.order_primers['p1']-1]
        self.primer_seq['F2'] = toSeq[self.order_primers['p2']-1:self.order_primers['p3']-1]
        self.primer_seq['F1c'] = toSeq[self.order_primers['p4']-1:self.order_primers['p5']-1]
        self.primer_seq['B1c'] = toSeq[self.order_primers['gap']-1:self.order_primers['p10']-1]
        self.primer_seq['B2'] = toSeq[self.order_primers['p9']-1:self.order_primers['p8']-1]
        self.primer_seq['B3'] = toSeq[self.order_primers['p7']-1:self.order_primers['p6']-1]
        return
    
    def Filter(self) -> int:
        '''过滤自身，若CG含量、Tm值、二聚物（尚未完成）均符合标准则返回1。其他则返回其他数字'''
        if self.primer_seq:
            for key in self.primer_seq:
                if conCG(self.primer_seq[key]) <0.4 or conCG(self.primer_seq[key]) >0.6:
                    return 2
                elif Tm(self.primer_seq[key]) > 66 or Tm(self.primer_seq[key]) <50:
                    return 3
                elif ifdimer(self.primer_seq[key]):
                    return 4
                else:
                    return 1
        else:
            return -1
        
    def init_verify(self) -> int:
        '''验证自身引物，若合规则返回1，不合规则返回0'''
        if self.init_primers['p3']+self.init_primers['p4'] > 60 or self.init_primers['p3']+self.init_primers['p4'] < 40:
            return 0
        elif self.init_primers['p8']+self.init_primers['p9'] > 60 or self.init_primers['p8']+self.init_primers['p9'] < 40:
            return 0
        elif self.init_primers['p4']+self.init_primers['p5']+self.init_primers['p9']+self.init_primers['p10']+self.init_primers['gap'] >160 or  \
            self.init_primers['p4']+self.init_primers['p5']+self.init_primers['p9']+self.init_primers['p10']+self.init_primers['gap']<120:
            return 0
        else:
            return 1

    def fitness(self):
        '''计算自己的适合度，之后存入self.blastscore'''
        try:
            self.blastscore
        except AttributeError:
            self.blastscore = 0
            # print('-',end='')
            job = []
            # mpp = []
            mpq = multiprocessing.Queue()
            for key in self.primer_seq:
                # if key in exist:
                #     mpq.put(exist[key])
                #     continue
            #     self.blastscore += blastn(self.primer_seq[key])
            # return self.blastscore
                # parent_conn, child_conn = multiprocessing.Pipe()
                # thisp = multiprocessing.Process(target=blastn,args=(self.primer_seq[key],key,mpq))
                thisp = threading.Thread(target=blastn,args=(self.primer_seq[key],mpq))
                # print(self.primer_seq[key])
                job.append(thisp)
                thisp.start()
                # mpp.append(parent_conn)

            for p in job:
                p.join()
            
            blascore = 0
            for p in job:
                blascore += mpq.get()
                # print(blascore)

            self.blastscore = blascore
            # print('-',end='')
            return self.blastscore
        else:
            return self.blastscore

        
    def mutate(self) -> int:
        '''自身变异,成功则返回1，失败返回0'''
        order = ['p1','p3','p5','p6','p8','p10','p2','p7','p4','p9','gap','start']
        ifmut = randomlist(0,1,12)
        for i in range(12):
            if ifmut[i]:
                if i < 6:
                    self.init_primers[order[i]] = random.randint(PRIMER_LENGTH[0],PRIMER_LENGTH[1])
                elif i < 8:
                    self.init_primers[order[i]] = random.randint(PRIMER_23GAP[0],PRIMER_23GAP[1])
                elif i == 8:
                    self.init_primers[order[i]] = random.randint(PRIMER_12CLEN[0]-self.init_primers['p3'],PRIMER_12CLEN[1]-self.init_primers['p3'])
                elif i == 9:
                    self.init_primers[order[i]] = random.randint(PRIMER_12CLEN[0]-self.init_primers['p8'],PRIMER_12CLEN[1]-self.init_primers['p8'])
                elif i == 10:
                    a = random.randint(PRIMER_CGAP[0]-self.init_primers['p4']-self.init_primers['p5']-self.init_primers['p9']-self.init_primers['p10'],
                                    PRIMER_CGAP[1]-self.init_primers['p4']-self.init_primers['p5']-self.init_primers['p9']-self.init_primers['p10'])
                    if a>0:
                        self.init_primers[order[i]] = a
                    else:
                        self.init_primers[order[i]] = 0
                else:
                    self.init_primers['start'] = random.randint(self.sstart,
                                        self.send-self.init_primers['p3']-self.init_primers['p4']-self.init_primers['p5']-self.init_primers['p8']-self.init_primers['p9']-self.init_primers['p10']-self.init_primers['p1']-self.init_primers['p2']-self.init_primers['p6']-self.init_primers['p7']-self.init_primers['gap'])
        self.code = primer2zip(self.init_primers)
        return 1

    def printself(self) -> None:
        print(self.init_primers)

    def getblastscore(self) -> float:
        return self.blastscore
    
    def getinit_primers(self):
        return self.init_primers

    def getdegbase(self):
        return self.degbase

    def codefresh(self):
        self.code = primer2zip(self.init_primers)
        return
    
    def getcode(self):
        return self.code

    def getprimersets(self):
        return 'fitness:%.2f | %d:%s,%d:%s,%d:%s,%d:%s,%d:%s,%d:%s' % (self.blastscore,self.order_primers['start'],self.primer_seq['F3'], \
                                                    self.order_primers['p2'],self.primer_seq['F2'], \
                                                    self.order_primers['p4'],self.primer_seq['F1c'], \
                                                    self.order_primers['gap'],self.primer_seq['B3'], \
                                                    self.order_primers['p9'],self.primer_seq['B2'], \
                                                    self.order_primers['p7'],self.primer_seq['B1c'])

    def __str__(self) -> str:
        return 'fitness:%.2f | %d:%s,%d:%s,%d:%s,%d:%s,%d:%s,%d:%s' % (self.blastscore,self.order_primers['start'],self.primer_seq['F3'], \
                                                    self.order_primers['p2'],self.primer_seq['F2'], \
                                                    self.order_primers['p4'],self.primer_seq['F1c'], \
                                                    self.order_primers['gap'],self.primer_seq['B3'], \
                                                    self.order_primers['p9'],self.primer_seq['B2'], \
                                                    self.order_primers['p7'],self.primer_seq['B1c'])
    # def degbase_generator(self,)





def hybrid_degbase(pa:primersets,pb:primersets) -> list:
    pass