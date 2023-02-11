#-*- coding = utf-8
# 主程序
from nucleotides import *
import time
from primer import *
import multiprocessing
from multiprocessing.managers import BaseManager

class Mymanager(BaseManager):
    pass

Mymanager.register('primersets',primersets)
manager = Mymanager()


def duplicate_primer(dup:primersets) -> primersets:
    a = manager.primersets(existprimer=dictdeepcopy(dup.getinit_primers()), existdegbase=listdeepcopy(dup.getdegbase()))
    return a

def prtmcal(pr):
    # pr.printself()
    pr.fitness()

if __name__ == "__main__":

    manager.start()
    time1 = time.time()
    f = open(file='log5.txt',mode='a')
    # 初始化
    p1 = manager.primersets()
    batch = 10 # 批次数
    elite = 5 # 精英数
    TOTALCYCLE = 5 # 总循环数
    
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
        # a.printself()

    for tejifje in range(TOTALCYCLE):
        # 计算适应度
        print('计算适应度')
        # for ps in solution:
        #     mtpool.apply_async(prtmcal,args=(ps,))
            # ps.fitness()
        mtpool = multiprocessing.Pool(4) # 进程管理池
        mtpool.map_async(prtmcal,solution)
        mtpool.close()
        mtpool.join()
        # print('zuse')
        # 排序
        ranked = quickSort(solution)

        f.write('-------'+str(tejifje)+'-------\n')
        # 输出部分信息
        for c in range(elite):
            f.write(ranked[c].getprimersets()+'\n')

    f.close()
    time2 = time.time()
    print(time2 - time1)