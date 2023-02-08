#-*- coding = utf-8
import re

def init_seq(in_seq:str) -> str:
    '''
    处理输入序列, 只留下大写ATCG
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
    temp_num = seq.replace('A','1').replace('T','2').replace('C','3').replace('G','4')
    comp = temp_num.replace('1','T').replace('2','A').replace('3','G').replace('4','C')
    return comp

def rev_com(seq:str) -> str:
    '''
    获得反向互补序列
    '''
    return reverse(complementory(seq))

se = 'AATTCCAGCCCGGAGCCTCAAT'

rec = complementory(se)

length = len(se)

# def cov(a:int) -> str:
#     if a == 1:
#         return '\\'
#     else:
#         return ' '

a = [[i for i in range(length)] for j in range(length)]

for i in range(length):
    for j in range(length):
        if (i-j)**2 < 9:
            a[i][j] = ' '
        else:
            a[i][j] = ('/' if se[i]==rec[j] else 'O')


for i in range(length+1):
    st = ''
    for j in range(length+1):
        if i == 0:
            if j == 0:
                st = st + ' '
            else:
                st = st + se[j-1]
        else:
            if j == 0:
                st = st + se[i-1]
            else:
                st = st + a[i-1][j-1]
    print(st)
