#!/usr/bin/env python3

import os
import sys
import math
        

def get_total_CN(cov_file):

    all_reg =[]
    for line in open(cov_file, "r"):
        line = line.strip().split()
        all_reg.append(line)

    av_nudt15_cov = float(all_reg[0][3])/(float(all_reg[0][2]) - float(all_reg[0][1]))
    av_vdr_cov = float(all_reg[1][3])/(float(all_reg[1][2]) - float(all_reg[1][1]))
    av_egfr_cov = float(all_reg[2][3])/(float(all_reg[2][2]) - float(all_reg[2][1]))
    av_e1_int2 = float(all_reg[3][3])/(float(all_reg[3][2]) - float(all_reg[3][1]))
    av_int2_e3 = float(all_reg[4][3])/(float(all_reg[4][2]) - float(all_reg[4][1]))

    av_ctrl_cov = (av_vdr_cov + av_egfr_cov)/2

    comp_av = av_nudt15_cov/av_ctrl_cov
    temp_cn = 2 * comp_av
    total_cn = round(temp_cn)


    return [str(int(total_cn)), round(av_nudt15_cov), round(av_ctrl_cov), str(av_e1_int2), str(av_int2_e3)];


def del_test(sv_del):

    if os.stat(sv_del).st_size == 0:
        return "None"

    else:
        for line in open(sv_del, "r"):
            if "COVERAGE" in line:
                line = line.strip().split()

                ABHom = line[-1]
                ABHet = line[-2]
                GT = line[2]
                DP = int(line[3])
        
                if float(ABHom) == 1.0:
                    return "*(full_gene_del)/*(full_gene_del)"
                elif float(ABHom) == -1.0:
                    return "*(full_gene_del)"
            else:
                pass


hap_adv_list = []
hap_t1 = []


def del_adv_test(hap_dbs, cand_allele1, cand_allele2, test_allele1, test_allele2, core_vars):
    g = open(hap_dbs, "r")
    for line in g:
        line = line.strip().split()
        hap_adv_list.append(line)

    a1 = core_vars.split(";")

    for i in a1:
        if i[-3:] == "0/1":
            hap_t1.append(i[:-4])                                                                                                                                

    for elem in hap_adv_list:
        if elem[1] == cand_allele1:
            list_t1 = (elem[2]).split(';')

        if elem[1] == cand_allele2:
            list_t2 = (elem[2]).split(';')

    if hap_t1[0] in list_t1:
        return test_allele1

    elif hap_t1[0] in list_t2:
        return test_allele2


het_hom_list = []
het_hom_list_new = []

def dup_test_init(sv_dup, av_cov):
    for line in open(sv_dup, "r"):
        if "COVERAGE" in line:
            continue
        elif "AGGREGATED" in line:
            continue

        else:
            fields = line.strip().split()
            het_hom_list.append(fields)

    test_list1 = []

    for i in het_hom_list:
        test_list1.append(int(i[2]))

    av_read_cov = sum(test_list1)/len(test_list1)
    norm_cov = (av_cov + av_read_cov)/2

    for i in het_hom_list:
        supp_reads = round(float(i[-2])*int(i[2]))
        i.append(round(supp_reads/norm_cov, 3))
        i.append(supp_reads)
        het_hom_list_new.append(i)


    return (het_hom_list_new)


hap_def_list = []
allele_cn_list = []

def dup_test_cn_3_4(sv_dup, hap_dbs, cand_allele1, cand_allele2, test_allele1, test_allele2, c_num, av_cov, in_list):

    g = open(hap_dbs, "r")
    for line in g:
        line = line.strip().split()
        hap_def_list.append(line)
        
            
    test_list1 = []
    test_list2 = []
    het_list = []


    for i in in_list:
        if i[1] == "0/1":
            het_list.append(i)

    for i in het_list:
        test_list1.append(i[0])
        test_list2.append(i[-2])


    max_het = max(test_list2)
    max_het_pos = test_list2.index(max_het)
    var = test_list1[max_het_pos]
    
    for elem in hap_def_list:
        if elem[1] == cand_allele1:
            list_3t = elem
            list_3t_2 = list_3t[2].split(';')
            l3 = len(list_3t_2)
            
        if elem[1] == cand_allele2:
            list_4t = elem
            list_4t_2 = list_4t[2].split(';')
            l4 = len(list_4t_2)

    hdb_list = list_3t_2 + list_4t_2


    index_var = hdb_list.index(var)

    if index_var < l3:
        allele_cn_list.append(test_allele1)
        allele_cn_list.append(int(round(max_het*int(c_num))))

    elif index_var >= l3:
        allele_cn_list.append(test_allele2)
        allele_cn_list.append(int(round(max_het*int(c_num))))

    
    if allele_cn_list[0] == test_allele1:
        rt_2 = int(c_num) - allele_cn_list[1]
        allele_cn_list.append(test_allele2)
        allele_cn_list.append(rt_2)

    elif allele_cn_list[0] == test_allele2:
        rt_2 = int(c_num) - allele_cn_list[1]
        allele_cn_list.append(test_allele1)
        allele_cn_list.append(rt_2)

    if allele_cn_list[1] == 0:
        res_dip = allele_cn_list[0] + "/" + allele_cn_list[2] + "x" + str(allele_cn_list[3] - 1)

    elif allele_cn_list[3] == 0:
        res_dip = allele_cn_list[2] + "/" + allele_cn_list[0] + "x" + str(allele_cn_list[1] - 1)

    elif allele_cn_list[1] == 1:
        res_dip = allele_cn_list[0] + "/" + allele_cn_list[2] + "x" + str(allele_cn_list[3])

    elif allele_cn_list[3] == 1:
        res_dip = allele_cn_list[2] + "/" + allele_cn_list[0] + "x" + str(allele_cn_list[1])  

    elif allele_cn_list[1] == 2:
        res_dip = allele_cn_list[0] + "x2" + "/" + allele_cn_list[2] + "x" + str(allele_cn_list[3])

    elif allele_cn_list[3] == 2:
        res_dip = allele_cn_list[2] + "x2" + "/" + allele_cn_list[0] + "x" + str(allele_cn_list[1])


    else:
        res_dip = 'check'

    return res_dip



def dup_test_cn_n(sv_dup, hap_dbs, cand_allele1, cand_allele2, test_allele1, test_allele2, c_num, av_cov, in_list):

    g = open(hap_dbs, "r")
    for line in g:
        line = line.strip().split()
        hap_def_list.append(line)


    test_list1 = []
    test_list2 = []
    het_list = []


    for i in in_list:
        if i[1] == "0/1":
            het_list.append(i)

    for i in het_list:
        test_list1.append(i[0])
        test_list2.append(i[-2])

    max_het = max(test_list2)
    max_het_pos = test_list2.index(max_het)
    var = test_list1[max_het_pos]


    for elem in hap_def_list:
        if elem[1] == cand_allele1:
            list_3t = elem
            list_3t_2 = list_3t[2].split(';')
            l3 = len(list_3t_2)

        if elem[1] == cand_allele2:
            list_4t = elem
            list_4t_2 = list_4t[2].split(';')
            l4 = len(list_4t_2)

    hdb_list = list_3t_2 + list_4t_2

    index_var = hdb_list.index(var)

    if index_var < l3:
        allele_cn_list.append(test_allele1)
        allele_cn_list.append(int(round(max_het*int(c_num)-0.15)))

    elif index_var >= l3:
        allele_cn_list.append(test_allele2)
        allele_cn_list.append(int(round(max_het*int(c_num)-0.15)))


    if allele_cn_list[0] == test_allele1:
        rt_2 = int(c_num) - allele_cn_list[1]
        allele_cn_list.append(test_allele2)
        allele_cn_list.append(rt_2)

    elif allele_cn_list[0] == test_allele2:
        rt_2 = int(c_num) - allele_cn_list[1]
        allele_cn_list.append(test_allele1)
        allele_cn_list.append(rt_2)

    if allele_cn_list[1] == 0:
        res_dip = allele_cn_list[0] + "/" + allele_cn_list[2] + "x" + str(allele_cn_list[3] - 1)

    elif allele_cn_list[3] == 0:
        res_dip = allele_cn_list[2] + "/" + allele_cn_list[0] + "x" + str(allele_cn_list[1] - 1)

    elif allele_cn_list[1] == 1:
        res_dip = allele_cn_list[0] + "/" + allele_cn_list[2] + "x" + str(allele_cn_list[3])

    elif allele_cn_list[3] == 1:
        res_dip = allele_cn_list[2] + "/" + allele_cn_list[0] + "x" + str(allele_cn_list[1])

    elif allele_cn_list[1] == 2:
        res_dip = allele_cn_list[0] + "x2" + "/" + allele_cn_list[2] + "x" + str(allele_cn_list[3])

    elif allele_cn_list[3] == 2:
        res_dip = allele_cn_list[2] + "x2" + "/" + allele_cn_list[0] + "x" + str(allele_cn_list[1])

    elif allele_cn_list[1] == 3:
        res_dip = allele_cn_list[0] + "x3" + "/" + allele_cn_list[2] + "x" + str(allele_cn_list[3])

    elif allele_cn_list[3] == 3:
        res_dip = allele_cn_list[2] + "x3" + "/" + allele_cn_list[0] + "x" + str(allele_cn_list[1])

    elif allele_cn_list[1] == 4:
        res_dip = allele_cn_list[0] + "x4" + "/" + allele_cn_list[2] + "x" + str(allele_cn_list[3])

    elif allele_cn_list[3] == 4:
        res_dip = allele_cn_list[2] + "x4" + "/" + allele_cn_list[0] + "x" + str(allele_cn_list[1])


    else:
        res_dip = 'check'

    return res_dip


def hybrid_sv1_test1(cov_e1_int2, cov_int2_e3):

    if 0.85 < float(cov_e1_int2)/float(cov_int2_e3) < 1.2:
        return 'norm_var'

    elif 0.45 < float(cov_e1_int2)/float(cov_int2_e3) < 0.75:
        return 'hyb_sv1'

    elif float(cov_e1_int2)/float(cov_int2_e3) < 0.15:
        return 'hyb_sv1_2'

    else:
        return 'norm_var'


def hybrid_new_test1(cov_e1_int2, cov_int2_e3):

    if 0.85 < float(cov_e1_int2)/float(cov_int2_e3) < 1.2:
        return 'norm_var'

    elif 0.45 < float(cov_int2_e3)/float(cov_e1_int2) < 0.75:
        return 'hyb_new'

    elif float(cov_int2_e3)/float(cov_e1_int2) < 0.15:
        return 'hyb_new_2'

    else:
        return 'norm_var'

