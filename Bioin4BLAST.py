#-*- coding = utf-8
# 这是处理LAMP引物设计的程序，同时引入BLAST结果作为考量工具。

# 初始化区
import os
import math
import re

# 全局变量


# 函数定义区
def init_seq(in_seq:str) -> str:
    '''
    处理输入序列，只留下大写ATCG
    '''
    up_seq = in_seq.upper()
    get = re.findall(r'[ATCG]*',up_seq)
    res = ''.join(get)
    return res

def reverse(seq:str) -> str:
    '''
    获得反向序列
    '''
    return seq[::-1]

def complementory(seq:str) -> str:
    '''
    获得互补序列
    '''
    # def com(char):
    #     N = {'A':'T','T':'A','C':'G','G':'C'}
    #     return N[char]
    temp_num = seq.replace('A','1').replace('T','2').replace('C','3').replace('G','4')
    comp = temp_num.replace('1','T').replace('2','A').replace('3','G').replace('4','C')
    # comp = map(com,seq)
    return comp

def rev_com(seq:str) -> str:
    '''
    获得反向互补序列
    '''
    return reverse(complementory(seq))

def dG(seq:str) -> float:
    '''
    计算序列自由能,kJ'''
    H, S = dHS(seq)
    T = Tm2(seq) + 273.15
    G = (H*1000 - T*S) /239
    return G

def dH(seq:str) -> float:
    '''
    计算序列焓变,kcal'''
    San_Han, San_Sang = dHS(seq)
    return San_Han

def dS(seq:str) -> float:
    '''
    计算序列熵变,cal'''
    San_Han, San_Sang = dHS(seq)
    return San_Sang

def dHS(seq:str) -> float:
    '''
    计算序列焓变和熵变'''
    SanSang = {'AA': -23.6, 'AT': -18.8, 'AC': -23, 'AG': -16.1, 'TA': -18.5, 'TC': -20.3, 'TG': -19.3, 'TT': -23.6, 'CG': -25.5, 'CA': -19.3, 'CT': -16.1, 'CC': -15.6, 'GA': -20.3, 'GT': -23, 'GC': -28.4, 'GG': -15.6}
    SanHan = {'AA': -8.4, 'AT': -6.5, 'AC': -8.6, 'AG': -6.1, 'TA': -6.3, 'TC': -7.7, 'TG': -7.4, 'TT': -8.4, 'CG': -10.1, 'CA': -7.4, 'CT': -6.1, 'CC': -6.7, 'GA': -7.7, 'GT': -8.6, 'GC': -11.1, 'GG': -6.7}
    San_Han = 0
    San_Sang = 0
    for i in range(len(seq)-1):
        NN = seq[i:i+2]
        San_Sang += SanSang[NN]
        San_Han += SanHan[NN]
    if seq.count('C') + seq.count('G') > 0:
        San_Sang += -5.9
    else:
        San_Sang += 0.6
    return San_Han, San_Sang

def Tm2(seq:str, kcl = 50, mg = 2, dntp = 0.2, primer = 0.4) -> float:
    adjusted = 16.6*math.log10((kcl+120*math.sqrt(mg-4*dntp))/1000)
    San_Han, San_Sang = dHS(seq)
    San = San_Han *1000/(San_Sang +1.987*math.log(primer*10**-6))-273.15+adjusted-3
    return San

def Tm(seq:str) -> float:
    '''
    计算序列Tm值
    '''
    SanHan = {'AA': -8.4, 'AT': -6.5, 'AC': -8.6, 'AG': -6.1, 'TA': -6.3, 'TC': -7.7, 'TG': -7.4, 'TT': -8.4, 'CG': -10.1, 'CA': -7.4, 'CT': -6.1, 'CC': -6.7, 'GA': -7.7, 'GT': -8.6, 'GC': -11.1, 'GG': -6.7}
    SanSang = {'AA': -23.6, 'AT': -18.8, 'AC': -23, 'AG': -16.1, 'TA': -18.5, 'TC': -20.3, 'TG': -19.3, 'TT': -23.6, 'CG': -25.5, 'CA': -19.3, 'CT': -16.1, 'CC': -15.6, 'GA': -20.3, 'GT': -23, 'GC': -28.4, 'GG': -15.6}
    kcl = 50
    mg = 2
    dntp = 0.2
    primer = 0.4
    adjusted = 16.6*math.log10((kcl+120*math.sqrt(mg-4*dntp))/1000)
    San_Han = 0
    San_Sang = 0
    for i in range(len(seq)-1):
        NN = seq[i:i+2]
        San_Han += SanHan[NN]
        San_Sang += SanSang[NN]

    if seq.count('C') + seq.count('G') > 0:
        San_Sang += -5.9
    else:
        San_Sang += 0.6
    
    San = San_Han *1000/(San_Sang +1.987*math.log(primer*10**-6))-273.15+adjusted-3
    return San

def printseqinfo(seq:str):
    '''
    打印序列信息'''
    s = init_seq(seq)
    print()

def nG(seq:str) -> float:
    '''新的计算自由能的方法,kcal'''
    Gibs = 1.96
    DG = {
        'AA' : -1 ,
        'AT' : -0.88 ,
        'TA' : -0.58 ,
        'CA' : -1.45 ,
        'GT' : -1.44 ,
        'CT' : -1.28 ,
        'GA' : -1.3 ,
        'CG' : -2.17 ,
        'GC' : -2.24 ,
        'GG' : -1.84 ,
        'TT' : -1 ,
        'TG' : -1.45 ,
        'AC' : -1.44 ,
        'AG' : -1.28 ,
        'TC' : -1.3 ,
        'CC' : -1.84 ,
    }
    for i in range(len(seq)-1):
        NN = seq[i:i+2]
        Gibs += DG[NN]
    if seq[-1] == 'T' or seq[-1] == 'A':
        Gibs += 0.05
    return Gibs


# 主程序


# 测试区
s = 'CATGCTAGCATTGCTAGCTAGCGCTAGCA'
s2 = 'cagtcgGCATCGjkaCgTiojo./'
seq = init_seq(s2)
# print(seq)
# print(Tm(init_seq(s2)))
# print(complementory(s2))
# print(dH(seq))
# print(dS(seq))
# print(Tm2(seq))
# print(Tm(seq))
# print(dG(seq))
s3 = ['GCCCAAACAACGACGATCGG','TAAAACCAGCATCCGTAGCCT',
'GCACGTTCTCCAACGGTGCT','GGTTGCTTGTTCAGCGAACT',
'CCAGAGAGGAGGTTGCCAA',
'TGCAGGATGCAGCGCCTTA',
'TGTCATGGGATCCGGATGTT',
'CAATAGAGCCGCTCTCAGAG']
for s4 in s3:
    print(init_seq(s4) + '  |  dG = '+ str(nG(init_seq(s4))) + '  |  Tm = ' + str(Tm2(init_seq(s4))))

