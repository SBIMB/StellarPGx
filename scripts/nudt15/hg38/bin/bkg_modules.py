#!/usr/bin/env python3


import os
import sys


def get_backgroud_alleles(database, core_vars):

    dbs = []
    dbs_temp = []

    core_vars_list = core_vars.split(";")
    core_temp1 = core_vars_list[-1][:-4]
    core_temp2 = core_vars_list[0][:-4]

    for line in open(database, "r"):
        line = line.strip().split("\t")
        dbs.append(line)

    for record in dbs:
        temp_rec = record[1]
        
        if core_temp1 and core_temp2 in temp_rec:
            dbs_temp.append(record)

            
    scores = []
    candidates = []
    cand_vars = []

    for elem in dbs_temp:
        candidates.append(elem[0])
        record_core_var = elem[1].split(";")
        cand_vars.append(record_core_var)

        counter = 0

        for i in record_core_var:
            if i in core_vars_list:
                counter += 3
            elif i[:-4] in core_vars:
                counter += 1
            else:
                counter += -2

        scores.append(counter)

    cand_diplos = []
    diplo_vars2 = []

    if len(scores) == 0:
        diplo1 = '1.v1_1.v1'
        allele_res = '*1/*1'

    else:
        max_score = max(scores)

        indices = [i for i, x in enumerate(scores) if x == max_score or x == max_score - 1]

        for i in indices:
            diplo = candidates[i]
            diplo_vars1 = len(cand_vars[i])
            cand_diplos.append(diplo)
            diplo_vars2.append(diplo_vars1)

        min_index = diplo_vars2.index(min(diplo_vars2))

        diplo1 =  cand_diplos[min_index]


        res1 = [i for i in range(len(diplo1)) if diplo1.startswith("_", i)]
        res2 = [i for i in range(len(diplo1)) if diplo1.startswith(".", i)]
        hap1 = "*" + str (diplo1[:res2[0]])
        hap2 = "*" + str (diplo1[res1[0]+1:res2[1]])
        allele_res =  hap1 + "/" + hap2 

    return [allele_res, diplo1];
