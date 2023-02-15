#-*- coding = utf-8
# 主程序
from nucleotides import *
import time
from primer import *
import multiprocessing
from multiprocessing.managers import BaseManager
import json

# 读取保存的引物组历史数据和序列查询数据
fileprimer = open('fileprimer.txt','r')
jsprimer = fileprimer.read()
existsprimer = json.loads(jsprimer)
fileprimer.close()

class Mymanager(BaseManager):
    pass

Mymanager.register('primersets',primersets)
# Mymanager.register('dict',dict)
manager = Mymanager()


def duplicate_primer(dup:primersets) -> primersets:
    a = manager.primersets(existprimer=dictdeepcopy(dup.getinit_primers()), existdegbase=listdeepcopy(dup.getdegbase()))
    return a

def prtmcal(pr):
    # pr.printself()
    return pr.getcode(),pr.fitness()

def saveinexist(callback,primer=existsprimer):
    existsprimer[callback[0]] = callback[1]
    

def hybrid_primer(pa:primersets,pb:primersets):
    '''两个引物杂交，成功返回1，重复5次不成功返回0'''
    order = ['start','p1','p2','p3','p4','p5','gap','p10','p9','p8','p7','p6']
    t = 5
    pra = pa.getinit_primers()
    prb = pb.getinit_primers()
    # print(ifswap)
    while t > 0:
        ifswap = randomlist(0,1,12)
        for i in range(12):
            if ifswap[i]:
                temp = pra[order[i]]
                pra[order[i]] = prb[order[i]]
                prb[order[i]] = temp
        psa = manager.primersets(existprimer=pra)
        psb = manager.primersets(existprimer=prb)
        if psa.init_verify() and psb.init_verify():
            psa.codefresh()
            psb.codefresh()
            # print(psa.getcode())
            # print(psb.getcode())
            return psa,psb
        else:
            t -= 1
    return 0,0

if __name__ == "__main__":

    manager.start()
    time1 = time.time()
    f = open(file='log5.txt',mode='a')
    # 初始化
    p1 = manager.primersets()
    batch = 10 # 批次数
    elite = 5 # 精英数
    TOTALCYCLE = 10 # 总循环数
    
    

    # filequery = open('filequery.txt','r')
    # jsquery = filequery.read()
    # existsquery = manager.dict(json.loads(jsquery))
    # filequery.close()
    # print(existsquery)
    # mngr = multiprocessing.Manager() # 进程管理器，传递py对象
    # manager = BaseManager()
    # manager.register('primersets',primersets)
    # manager.start()

    # 导入序列
    fasta = 'conseq.fasta'
    seqdict = readfasta(fasta)
    seq = seqdict['Consensus']

    # 前处理
    # solution = mngr.list()
    solution = [p1 for i in range(batch)]
    hybrid = [p1 for i in range(batch)]
    mutant = [p1 for i in range(batch)]

    # 生成最初解

    for i in range(batch):
        t = 0
        # a = primersets(1,len(seq))
        while t != 1:
            a = manager.primersets(1,len(seq))
            a.getPrimerSeq(seq)
            t = a.Filter()
        # a.fitness()
        solution[i] = a
        # existsprimer[a.getcode()]
        # a.printself()

    for tejifje in range(TOTALCYCLE):
        # 计算适应度
        print('-----第'+str(tejifje)+'次-----')
        print('计算适应度')
        # for ps in solution:
        #     mtpool.apply_async(prtmcal,args=(ps,))
            # ps.fitness()
        mtpool = multiprocessing.Pool(4) # 进程管理池
        for p in solution:
            mtpool.apply_async(prtmcal,(p,),callback=saveinexist)
        # mtpool.map_async(prtmcal,solution)
        mtpool.close()
        mtpool.join()
        print(existsprimer)
        # print('zuse')
        # 排序
        ranked = quickSort(solution)
        print(time.time()-time1)
        f.write('-------'+str(tejifje)+'-------\n')
        # 输出部分信息
        for c in range(elite):
            f.write(ranked[c].getprimersets()+'\n')
        print('杂交')
        # 杂交
        i = 0
        while True:
            if i < batch:
                ran1 = random.randint(0,elite-1)
                # print('try 1')
                ran2 = random.randint(0,elite-1)
                while ran2 == ran1:
                    ran2 = random.randint(0,elite-1)
                p1 = duplicate_primer(solution[ran1])
                p2 = duplicate_primer(solution[ran2])
                p1,p2 = hybrid_primer(p1,p2)
                if p1:
                    # print(p1.getcode())
                    if p1.getcode() not in existsprimer:
                        p1.getPrimerSeq(seq)
                        hybrid[i] = p1
                        i += 1
                        # print('succsses 1')
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
                    if m.init_verify() and m.getcode() not in existsprimer:
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
                a = manager.primersets(1,len(seq))
                while t != 1:
                    a = manager.primersets(1,len(seq))
                    a.getPrimerSeq(seq)
                    t = a.Filter() and a.getcode() not in existsprimer
                solution[i] = a

    # print(existsprimer)
    # print(existsquery)

    fileprimer = open('fileprimer.txt','w')
    jsprimer = json.dumps(existsprimer)
    fileprimer.write(jsprimer)
    fileprimer.close()

    # anodict = {}
    # for key in existsquery:
    #     anodict[key] = existsquery[key]
    # filequery = open('filequery.txt','w')
    # jsquery = json.dumps(anodict)
    # filequery.write(jsquery)
    # filequery.close()
    f.close()
    time2 = time.time()
    print(time2 - time1)