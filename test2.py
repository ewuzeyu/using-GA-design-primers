import random
import primer3
from primer3 import thermoanalysis

seq = 'ACATGCATCGGGCTATTGCGGCTAGCTAGCTAGCTGCTAC'

en = primer3.bindings.calcEndStability()
Ha = primer3.calcHairpin(seq)

print(Ha.structure_found)

# class tes:
#     def __init__(self) -> None:
#         self.num = random.randint(1,100)

# l = []
# for i in range(50):
#     l.append(tes())

# for i in l:
#     print(i.num)
# print('排序中')
# def quickSort(listx):
#     if len(listx)<=1:
#         return listx
#     pivot = listx[len(listx)//2]
#     listl = [x for x in listx if x.num < pivot.num]
#     listm = [x for x in listx if x.num ==pivot.num]
#     listr = [x for x in listx if x.num > pivot.num]
#     left = quickSort(listl)
#     right = quickSort(listr)
#     return left + listm + right

# a = quickSort(l)

# for i in l:
#     print(i.num)
# print('---')
# for i in a:
#     print(i.num)