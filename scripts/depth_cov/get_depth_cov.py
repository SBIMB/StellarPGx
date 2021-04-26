#!/usr/bin/env python3

import os
import sys

res_file = sys.argv[1]
depth_file = sys.argv[2]

f = open(res_file, "r")
g = open(depth_file, "r")

orig = []
pos = []

for line in f:

    line = line.strip().split()
    orig.append(line)
   
    a = line[0] + ":" + line[1]
    b = line[0] + ":" + line[2]
    c = a + "-" + b
    
    list1 = c.split("-")

    pos.append(list1)


data_list = []

for line in g:
    line = line.strip().split()

    if ":" not in line[0]:
        line[0] = line[0] + ":" + line[1]
        line.pop(1)

    else:
        pass

    data_list.append(line)
    

flag = False

reg_list = []

for i in range(len(pos)):

    reg_sublist = []
    pos_cov = []

    for elem in data_list:

        if elem[0] == (pos[i][0]):
            ind1 = data_list.index(elem)
            
        if elem[0] == (pos[i][1]):
            ind2 = data_list.index(elem)

    reg_sublist = data_list[ind1:ind2 + 1]


    for elem in reg_sublist:
        pos_cov.append(int(elem[1]))

    reg_list.append(sum(pos_cov))

out_cov = []

for elem in orig:
    ind = orig.index(elem)
    elem.append(str(reg_list[ind]))
    print("\t".join(elem))
