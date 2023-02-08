#-*- coding = utf-8
# 主程序
from nucleotides import *
import time
from primer import *

time1 = time.time()
f = open(file='log2.txt',mode='a')
# 初始化
p1 = primersets()
batch = 100 # 批次数
elite = 10 # 精英数
TOTALCYCLE = 10 # 总循环数

# 导入序列
fasta = 'conseq.fasta'
seqdict = readfasta(fasta)
seq = seqdict['Consensus']

# 前处理
solution = [p1 for i in range(batch)]
hybrid = [p1 for i in range(batch)]
mutant = [p1 for i in range(batch)]

# 生成最初解

for i in range(batch):
    t = 0
    a = primersets(1,len(seq))
    while t != 1:
        a = primersets(1,len(seq))
        a.getPrimerSeq(seq)
        t = a.Filter()
    # a.fitness()
    solution[i] = a

for tejifje in range(TOTALCYCLE):
    # 计算适应度
    print('计算适应度')
    for ps in solution:
        ps.fitness()

    # 排序
    ranked = quickSort(solution)

    f.write('-------'+str(tejifje)+'-------\n')
    # 输出部分信息
    for c in range(5):
        f.write(ranked[c].__str__()+'\n')
    print('杂交')
    # 杂交
    i = 0
    while True:
        if i < batch:
            ran1 = random.randint(0,elite-1)
            ran2 = ran1
            while ran2 == ran1:
                ran2 = random.randint(0,elite-1)
            p1 = duplicate_primer(solution[ran1])
            p2 = duplicate_primer(solution[ran2])
            if hybrid_primer(p1,p2):
                p1.getPrimerSeq(seq)
                hybrid[i] = p1
                i += 1
        else:
            break
    print('变异')
    # 变异
    i = 0
    while True:
        if i < batch:
            m = duplicate_primer(solution[i])
            if m.mutate():
                m.getPrimerSeq(seq)
                if m.init_verify():
                    m.getPrimerSeq(seq)
                    mutant[i] = m
                    i += 1
        else:
            break
    print('合并')
    # 合并
    for i in range(batch):
        if i < elite:
            solution[i] = ranked[i]
        elif i < (batch-elite)/2:
            solution[i] = hybrid[random.randint(0,batch-1)]
        elif i < batch - elite:
            solution[i] = mutant[random.randint(0,batch-1)]
        else:
            t = 0
            a = primersets(1,len(seq))
            while t != 1:
                a = primersets(1,len(seq))
                a.getPrimerSeq(seq)
                t = a.Filter()
            solution[i] = a

f.close()
time2 = time.time()
print(time2 - time1)