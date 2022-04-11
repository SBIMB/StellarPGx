#!/usr/bin/env python3

import os
import sys
import math
        

def get_total_CN(cov_file):

    all_reg =[]
    for line in open(cov_file, "r"):
        line = line.strip().split()
        all_reg.append(line)

    av_2d6_cov = float(all_reg[2][3])/(float(all_reg[2][2]) - float(all_reg[2][1]))
    av_vdr_cov = float(all_reg[3][3])/(float(all_reg[3][2]) - float(all_reg[3][1]))
    av_in1_3pr = float(all_reg[1][3])/(float(all_reg[1][2]) - float(all_reg[1][1]))
    av_ex9_3pr = float(all_reg[0][3])/(float(all_reg[0][2]) - float(all_reg[0][1]))
    av_5pr_in1 = float(all_reg[4][3])/(float(all_reg[4][2]) - float(all_reg[4][1]))
    av_in4_3pr = float(all_reg[5][3])/(float(all_reg[5][2]) - float(all_reg[5][1]))
    av_5pr_in4 = float(all_reg[6][3])/(float(all_reg[6][2]) - float(all_reg[6][1]))    
    av_2d7_ex9 = float(all_reg[7][3])/(float(all_reg[7][2]) - float(all_reg[7][1]))
    av_2d7_in4_in8 = float(all_reg[8][3])/(float(all_reg[8][2]) - float(all_reg[8][1]))
    av_egfr_cov = float(all_reg[9][3])/(float(all_reg[9][2]) - float(all_reg[9][1]))
    av_2d7_ex2_in8 = float(all_reg[10][3])/(float(all_reg[10][2]) - float(all_reg[10][1]))
    av_2d7_5pr_in1 = float(all_reg[11][3])/(float(all_reg[11][2]) - float(all_reg[11][1]))

    av_ctrl_cov = (av_vdr_cov + av_egfr_cov)/2
    # av_ctrl_cov = av_vdr_cov

    comp_av = av_2d6_cov/av_ctrl_cov
    temp_cn = 2 * comp_av + 0.15
    total_cn = round(temp_cn)

    in1_3pr = round(2 * av_in1_3pr/av_ctrl_cov) 
    ex9_3pr = (2 * av_ex9_3pr/av_ctrl_cov)


    return [str(total_cn), round(av_2d6_cov), str(int(in1_3pr)), round(av_ctrl_cov), str(ex9_3pr), round(av_in1_3pr), str(av_in4_3pr), str(av_5pr_in4), str(av_2d7_ex9), str(av_2d7_in4_in8), str(av_2d7_ex2_in8), str(av_2d7_5pr_in1)];


samp_gt = ""
samp_gt_hap1 = ""

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
                    return "*5/*5"
                elif float(ABHom) == -1.0:
                    return "*5"
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
    test_list3 = []
    het_list = []


    for i in in_list:
        if i[1] == "0/1":
            het_list.append(i)

    # return het_list

    # if len(het_list) > 1 and het_list[0][0] == "42522613~G>C":
    #     het_list.pop(0)
    # else:
    #     pass

    for i in het_list:
        test_list1.append(i[0])
        test_list2.append(i[-2])
        test_list3.append(float(i[-4]))

    max_het = max(test_list2)
    
    # if max_het > 1:
    #     max_het = min(test_list2)
    # else:
    #     pass

    max_het_pos = test_list2.index(max_het)
    var = test_list1[max_het_pos]

    if max_het > 1:
        max_het = test_list3[max_het_pos]
    elif max_het > test_list3[max_het_pos]:
        max_het = test_list3[max_het_pos]
    else:
        pass


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

    # return hdb_list

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


    # return allele_cn_list

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

    elif allele_cn_list[3] < -1:
        res_dip =  allele_cn_list[2] + "/" + allele_cn_list[0] + "x" + str(3)

    elif allele_cn_list[1] < -1:
        res_dip =  allele_cn_list[0] + "/" + allele_cn_list[2] + "x" + str(3)


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

    # if len(het_list) > 1 and het_list[0][0] == "42522613~G>C":
    #     het_list.pop(0)
    # else:
    #     pass


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

    elif allele_cn_list[3] < 0:
        res_dip =  allele_cn_list[2] + "/" + allele_cn_list[0] + "x" + str(allele_cn_list[1] + allele_cn_list[3] - 1)

    elif allele_cn_list[1] < 0:
        res_dip = allele_cn_list[0] + "/" + allele_cn_list[2] + "x" + str(allele_cn_list[3] + allele_cn_list[1] - 1)


    else:
        res_dip = 'check'

    return res_dip


def hybrid_test_68(sv_dup, c_num, av_cov, cn_in1_3pr1, in_list):
    
    test_list1 = []
    test_list2 = []
    test_list3 = []


    for i in in_list:
        test_list1.append(i[0])
        test_list2.append(abs(float(i[-2])))
        test_list3.append(i[-1])

    index1 = test_list1.index('42526694~G>A')
    index2 = test_list1.index('42524947~C>T')
    
    val_68 = test_list3[index1]
    val_4 = test_list3[index2]
    
    rt = val_68/val_4

    
    if rt <= 1.4:
        return 'norm_dup'

    elif rt > 1.4:
        return 'hyb_68'

       
    else:
        return 'norm_dup'


def hyb_test_5_68_4(sv_del, in1_3pr1_float, av_cov):
    test_del = []
    for line in open(sv_del, "r"):
        if "COVERAGE" in line:
            test_del.append(line.strip())

    t1 = 2 * in1_3pr1_float/av_cov

    if len(test_del) == 0 and (1.6 < t1 < 2.8):
        return 'norm_art'

    elif len(test_del) > 0 and t1 < 1.6:
        return 'del_hyb'


def hybrid_test_36(sv_dup, cn, av_cov, cn_ex9_3pr):

    if int(round(float(cn_ex9_3pr))) == int(cn):
        return 'norm_dup'

    elif ((int(cn) - 1) - 0.35) < float(cn_ex9_3pr) < ((int(cn) - 1) + 0.5):
        return 'hyb_36_10'
    elif (int(cn) - 2) <= float(cn_ex9_3pr) < (int(cn) - 2 + 0.65):
        return 'hyb_36_36'


def hybrid_test_36_single(sv_dup, cn, av_cov, cn_ex9_3pr):

    if int(round(float(cn_ex9_3pr))) == int(cn):
        return 'norm_star10'

    elif ((int(cn) - 1) - 0.35) < float(cn_ex9_3pr) < ((int(cn) - 1) + 0.5):
        return 'hyb_36_single'

    else:
        return 'norm_star10'




def hybrid_test_36_mod(sv_dup, cn, av_cov, cn_ex9_3pr):

    if int(round(float(cn_ex9_3pr))) == int(cn):
        return 'norm_mt'
                                                                                                  
    elif ((int(cn) - 1) - 0.3) < float(cn_ex9_3pr) < ((int(cn) - 1) + 0.5):
        return 'hyb_36_10'
    elif ((int(cn) - 2) - 0.3) <= float(cn_ex9_3pr) < (int(cn) - 2 + 0.7):
        return 'hyb_36_36'


def hybrid_test_36_multi(sv_dup, cn, av_cov, cn_ex9_3pr):

    if int(round(float(cn_ex9_3pr))) == int(cn):
        return 'norm_mt'

    elif ((int(cn) - 1) - 0.05) < float(cn_ex9_3pr) < ((int(cn) - 1) + 0.5):
        return 'hyb_36_10'
    elif (int(cn) - 2) <= float(cn_ex9_3pr) < (int(cn) - 2 + 0.95):
        return 'hyb_36_36'
    elif (int(cn) - 3) <= float(cn_ex9_3pr) < (int(cn) - 3 + 0.95):
        return 'hyb_36_36_36'
    else:
        return 'check'


def hybrid_test_36_multi_10(sv_dup, cn, av_cov, cn_ex9_3pr, cn_star10):

    if int(round(float(cn_ex9_3pr))) == int(cn):
        return 'norm_mt'

    elif float(cn_ex9_3pr) < ((int(cn) - 1) + 0.5):
        cn_star36 = int(cn) - int(round(float(cn_ex9_3pr)))
        adj_cn_star10 = int(cn_star10) - cn_star36

        if cn_star36 == 1:
            return '*36+*10x' + str(adj_cn_star10)
        else:
            return '*36x' + str(cn_star36) + '+*10x'  + str(adj_cn_star10)

    else:
        return 'check'


def hybrid_13_2_v1(cov_in4_3pr, cov_5pr_in4):

    if 0.85 < float(cov_in4_3pr)/float(cov_5pr_in4) < 1.2:
        return 'norm_var'
    elif 0.45 < float(cov_in4_3pr)/float(cov_5pr_in4) < 0.75:
        return 'hyb_13_2'
    else:
        return 'norm_var'


def hybrid_13_2_v2(cov_2d7_ex2_in8, cov_2d7_5pr_in1):

    if 0.85 < float(cov_2d7_ex2_in8)/float(cov_2d7_5pr_in1) < 1.2:
        return 'norm_var'
    elif 0.45 < float(cov_2d7_ex2_in8)/float(cov_2d7_5pr_in1) < 0.75:
        return 'hyb_13_2_v2'
    else:
        return 'norm_var'



def tandem_90_1(in_list, alt_allele, cn):

    test_list1 = []
    test_list2 = []
    test_list3 = []

    for i in in_list:
        test_list1.append(i[0])
        test_list2.append(abs(float(i[-2])))
        test_list3.append(i[-1])


    if len(test_list1) > 1:
        index1 = test_list1.index('42525100~T>C')
        a = test_list3[index1]
        test_list3.pop(index1)
        b = max(test_list3)

        c = round(b/a)

        if int(cn) == 3 and c == 1:
            res = alt_allele + "/" + "*90+*1"

        elif int(cn) == 3 and c > 1:
            res = alt_allele + "x2" + "/" + "*90"

        elif int(cn) == 4 and c == 2:
            res = alt_allele + "x2" + "/" + "*90+*1"

        elif int(cn) == 4 and c >= 3:
            res = alt_allele + "x3" + "/" + "*90"

    else:
        val1 = test_list2[0]
        val2 = round(val1 * int(cn))

        if int(cn) == 3 and val2 == 1:
            res = '*1/*90+*1'

    return res


def tandem_57_10(in_list, alt_allele, cn):

    test_list1 = []
    test_list2 = []
    test_list3 = []

    for i in in_list:
        test_list1.append(i[0])
        test_list2.append(abs(float(i[-2])))
        test_list3.append(i[-1])

    if len(test_list1) > 1:
        index1 = test_list1.index('42525908~G>A')
        a = test_list3[index1]
        test_list3.pop(index1)
        index2 = test_list1.index('42526694~G>A')
        m = test_list3[index2]
        test_list3.pop(index2)
        b = max(test_list3)

        c = round(b/a)
        p = round(m/a)

        if int(cn) == 3 and c == 1 and p > 1:
            res = alt_allele + "/" + "*57+*10"

        elif int(cn) == 3 and c > 1 and p == 1:
            res = alt_allele + "x2" + "/" + "*57"

        elif int(cn) == 4 and c == 2 and p > 1:
            res = alt_allele + "x2" + "/" + "*57+*10"

        elif int(cn) == 4 and c >= 3 and p == 1:
            res = alt_allele + "x3" + "/" + "*57"

        elif int(cn) == 4 and p == 1 and alt_allele == '*10':
            res = "*57+*10" + "/" + "*57+*10"

        elif int(cn) == 4 and p > 1 and alt_allele == '*10':
            res = "*10x2" + "/" + "*57+*10"


        else:
            res = alt_allele + "/" + "*57+*10"

    return res


def hybrid_test_83_single(sv_dup, cn, av_cov, cn_ex9_3pr):

    if int(round(float(cn_ex9_3pr))) == int(cn):
        return 'norm_star39'

    elif ((int(cn) - 1) - 0.35) < float(cn_ex9_3pr) < ((int(cn) - 1) + 0.5):
        return 'hyb_83_single'

    else:
        return 'norm_star39'


def hybrid_test_83(sv_dup, cn, av_cov, cn_ex9_3pr):

    if int(round(float(cn_ex9_3pr))) == int(cn):
        return 'norm_star39'

    elif ((int(cn) - 1) - 0.35) < float(cn_ex9_3pr) < ((int(cn) - 1) + 0.5):
        return 'hyb_83'

    else:
        return 'norm_star39'

