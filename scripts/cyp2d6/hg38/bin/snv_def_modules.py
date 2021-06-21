#!/usr/bin/env python3

import os
import sys


def get_core_variants(infile, cn):
    core_vars = []
    for line in open(infile, "r"):
        line = line.strip()
        core_vars.append(line)
    core_vars = ";".join(sorted(core_vars))

    if int(cn) == 1:
        core_vars = core_vars.replace("~0/1", "~1/1")

    return core_vars

def get_all_vars_gt(infile_full_gt):
    all_vars_gt = []
    for line in open(infile_full_gt, "r"):
        line = line.strip()
        all_vars_gt.append(line)
    all_vars_gt = ";".join(sorted(all_vars_gt))
    return all_vars_gt

def cand_snv_allele_calling(database, infile, infile_full, infile_full_gt, infile_spec, cn):
    

    f = open(infile_spec, "r")

    all_variants = []

    for line in open(infile_full, "r"):
        line.strip()
        all_variants.append(line)

    if os.stat(infile).st_size == 0:
        cand_res = ['1.v1_1.v1']
        allele_res = "*1/*1"
        return ["".join(cand_res), "".join(cand_res), allele_res];
        sys.exit()

    core_variants = get_core_variants(infile, cn)

    all_var_gt = []
    for line in open(infile_full_gt, "r"):
        line = line.strip()
        all_var_gt.append(line)

    
    dbs = []

    for line in open(database, "r"):
        line = line.strip().split("\t")
        dbs.append(line)

    soln_list1 = []
    soln_list2 = []

    for record in dbs:
        record_core_var = record[1].split(";")
        record_core_var = ";".join(sorted(record_core_var))
        if record_core_var == core_variants:
            diplo = record[0]
            full_dip = record[2]
            soln_list1.append(record[0])
            soln_list2.append(record[2])
        else:
            pass


    diff_alleles_check = False

    def chkList(lst):
        if len(lst) < 0 :
            diff_alleles_check = True
        diff_alleles_check = all(ele == lst[0] for ele in lst)

        if(diff_alleles_check):
            return("Equal")
        else:
            return("Not equal")


    if len(soln_list1) == 1:
        diplo = "".join(soln_list1)
        res1 = [i for i in range(len(diplo)) if diplo.startswith("_", i)]
        res2 = [i for i in range(len(diplo)) if diplo.startswith(".", i)]
        hap1 = "*" + str (diplo[:res2[0]])
        hap2 = "*" + str (diplo[res1[0]+1:res2[1]])
        allele_res = hap1 + "/" + hap2
        return [soln_list1, diplo, allele_res];
        #print ("\nSupporting variants:")
        #print ("\n" + core_variants + "\n")


    elif len(soln_list1) == 2:
        diplo1 = soln_list1[0]
        diplo2 = soln_list1[1]
        diplo1_supp_var = soln_list2[0].split(";")
        diplo2_supp_var = soln_list2[1].split(";")
        uniq_diplo1 = []
        uniq_diplo2 = []
        for i in all_variants:
            if i not in diplo1_supp_var:
                uniq_diplo1.append(i)
        
            if i not in diplo2_supp_var:
                uniq_diplo2.append(i)

            
        if len(uniq_diplo1) < len(uniq_diplo2):
            res1 = [i for i in range(len(diplo1)) if diplo1.startswith("_", i)]
            res2 = [i for i in range(len(diplo1)) if diplo1.startswith(".", i)]
            hap1 = "*" + str (diplo1[:res2[0]])
            hap2 = "*" + str (diplo1[res1[0]+1:res2[1]])
            allele_res =  hap1 + "/" + hap2 
            return [soln_list1, diplo1, allele_res];
            #print ("Supporting variants:")
            #print ("\n" + core_variants + "\n")

        elif len(uniq_diplo1) > len(uniq_diplo2):
            res1 = [i for i in range(len(diplo2)) if diplo2.startswith("_", i)]
            res2 = [i for i in range(len(diplo2)) if diplo2.startswith(".", i)]
            hap1 = "*" + str (diplo2[:res2[0]])
            hap2 = "*" + str (diplo2[res1[0]+1:res2[1]])
            allele_res =  hap1 + "/" + hap2 
            return [soln_list1, diplo2, allele_res];
            #print ("Supporting variants:")
            #print ("\n" + core_variants + "\n")

        elif len(uniq_diplo1) == len(uniq_diplo2) and (diplo1 == "4.v11_74.v1" and diplo2 == "4.v12_1.v1"):
            res1 = [i for i in range(len(diplo2)) if diplo2.startswith("_", i)]
            res2 = [i for i in range(len(diplo2)) if diplo2.startswith(".", i)]
            hap1 = "*" + str (diplo2[:res2[0]])
            hap2 = "*" + str (diplo2[res1[0]+1:res2[1]])
            allele_res =  hap1 + "/" + hap2
            return [soln_list1, diplo2, allele_res];
    
        # elif len(uniq_diplo1) == len(uniq_diplo2) and diplo1 == "1.v1_2.v1" and diplo2 == "34.v1_39.v1":
        #     res1 = [i for i in range(len(diplo2)) if diplo2.startswith("_", i)]
        #     res2 = [i for i in range(len(diplo2)) if diplo2.startswith(".", i)]
        #     hap1 = "*" + str (diplo2[:res2[0]])
        #     hap2 = "*" + str (diplo2[res1[0]+1:res2[1]])
        #     allele_res =  hap1 + "/" + hap2 
        #     return [soln_list1, allele_res];
            #print ("Supporting variants:")
            #print ("\n" + core_variants + "\n")

        elif len(uniq_diplo1) == len(uniq_diplo2) and (diplo1 == "4.v1_6.v1" and diplo2 == "4.v4_6.v2") :
            res1 = [i for i in range(len(diplo1)) if diplo1.startswith("_", i)]
            res2 = [i for i in range(len(diplo1)) if diplo1.startswith(".", i)]
            hap1 = "*" + str (diplo1[:res2[0]])
            hap2 = "*" + str (diplo1[res1[0]+1:res2[1]])
            allele_res =  hap1 + "/" + hap2 
            return [soln_list1, diplo1, allele_res];
            #print ("Supporting variants:")
            #print ("\n" + core_variants + "\n")

    
        else:
            tiebreak1 = []
            tiebreak2 = []
            tiebreak3 = []
            score = []
            for line in f:
                line = line.strip().split()
                #print(line)
                if line[2] == core_variants:
                    tiebreak1.append(line[1])
                    tiebreak2.append(line[3])
                    tiebreak3.append(line[0])
            for full_dip in tiebreak2:
                diplo_supp_gt = full_dip.split(";")
                uniq_gt = []
                for i in all_var_gt:
                    if i not in diplo_supp_gt:
                        uniq_gt.append(i)
                score_dip = len(uniq_gt)
                score.append(score_dip)

            min_score = min(score)
            max_score = max(score)

            if chkList(score) == "Equal": # and soln_list1[0] != "1.v1_4.v1" and soln_list1[1] != "39.v1_4.v5":
                if soln_list1[0] == "10.v1_106.v1" and soln_list1[1] == "1.v1_52.v1":
                    if score == [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]:
                        elem = "10.v1_106.v1"
                        res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                        res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                        hap1 = "*" + str (elem[:res2[0]])
                        hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                        result_dip = hap1 + "/" + hap2
                        return [soln_list1, elem, result_dip];

                    else:
                        elem = "1.v1_52.v1"
                        res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                        res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                        hap1 = "*" + str (elem[:res2[0]])
                        hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                        result_dip = hap1 + "/" + hap2
                        return [soln_list1, elem, result_dip];


                elif soln_list1[0] == "119.v1_2.v1" and soln_list1[1] == "1.v1_41.v1":
                    elem = "1.v1_41.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

                elif soln_list1[0] == "106.v1_4.v1" and soln_list1[1] == "4.v8_52.v1":
                    elem = "106.v1_4.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

                elif soln_list1[0] == "10.v1_34.v1" and soln_list1[1] == "1.v1_65.v1":
                    elem = "10.v1_34.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

                elif soln_list1[0] == "114.v1_2.v1" and soln_list1[1] == "14.v1_65.v1":
                    elem = "114.v1_2.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

                elif soln_list1[0] == "1.v1_102.v1" and soln_list1[1] == "2.v1_48.v1":
                    elem = "2.v1_48.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

                elif soln_list1[0] == "1.v1_4.v5" and soln_list1[1] == "34.v1_4.v2":
                    elem = "1.v1_4.v5"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

                elif soln_list1[0] == "10.v1_4.v5" and soln_list1[1] == "4.v2_65.v1":
                    elem = "10.v1_4.v5"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

                elif soln_list1[0] == "1.v1_4.v6" and soln_list1[1] == "39.v1_4.v9":
                    elem = "1.v1_4.v6"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

                elif soln_list1[0] == "10.v1_14.v1" and soln_list1[1] == "114.v1_39.v1":
                    elem = "10.v1_14.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

                elif soln_list1[0] == "10.v1_115.v1" and soln_list1[1] == "132.v1_9.v1":
                    elem = "10.v1_115.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

                elif soln_list1[0] == "139.v1_4.v2" and soln_list1[1] == "1.v1_4.v11":
                    elem = "139.v1_4.v2"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];


                elif soln_list1[0] == "1.v1_2.v1" and soln_list1[1] == "34.v1_39.v1":
                    if score == [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]:
                        amb_soln_set = []
                        temp_set = []
                        temp_set.append("1.v1_2.v1")
                        temp_set.append("34.v1_39.v1")

                        for elem in temp_set:
                            res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                            res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                            hap1 = "*" + str (elem[:res2[0]])
                            hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                            result_dip = hap1 + "/" + hap2
                            amb_soln_set.append(result_dip)

                        if amb_soln_set[0] != amb_soln_set[1]:
                            allele_res =  " or ".join(amb_soln_set)
                        else:
                            allele_res = amb_soln_set[0]

                        return [soln_list1, allele_res];

                    else:
                        elem = "1.v1_2.v1"
                        res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                        res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                        hap1 = "*" + str (elem[:res2[0]])
                        hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                        allele_res = hap1 + "/" + hap2
                        return [soln_list1, elem, allele_res];

                elif soln_list1[0] == "1.v1_4.v5" and soln_list1[1] == "34.v1_4.v2":
                    if score == [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]:
                        amb_soln_set = []
                        temp_set = []
                        temp_set.append("1.v1_4.v5")
                        temp_set.append("34.v1_4.v2")

                        for elem in temp_set:
                            res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                            res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                            hap1 = "*" + str (elem[:res2[0]])
                            hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                            result_dip = hap1 + "/" + hap2
                            amb_soln_set.append(result_dip)

                        if amb_soln_set[0] != amb_soln_set[1]:
                            allele_res =  " or ".join(amb_soln_set)
                        else:
                            allele_res = amb_soln_set[0]

                        return [soln_list1, allele_res];

                    else:
                        elem = "1.v1_4.v5"
                        res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                        res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                        hap1 = "*" + str (elem[:res2[0]])
                        hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                        allele_res = hap1 + "/" + hap2
                        return [soln_list1, elem, allele_res];


                else:

                    amb_soln_set = []
                    for elem in soln_list1:
                        res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                        res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                        hap1 = "*" + str (elem[:res2[0]])
                        hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                        result_dip = hap1 + "/" + hap2
                        amb_soln_set.append(result_dip)
                        # elem_pos = tiebreak1.index(elem)
                        # print ("Solution " + str(elem_pos) + ": " + result_dip)
                        # allele_res =  " or ".join(amb_soln_set) 
                    if amb_soln_set[0] != amb_soln_set[1]:
                        allele_res =  " or ".join(amb_soln_set)
                    else:
                        allele_res = amb_soln_set[0]

                    return [soln_list1, soln_list1, allele_res];



            elif score.count(min_score) > 1 and soln_list1[0] == "1.v1_2.v1" and soln_list1[1] == "34.v1_39.v1":
                if score == [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]:
                    amb_soln_set = []
                    temp_set = []
                    temp_set.append("1.v1_2.v1")
                    temp_set.append("34.v1_39.v1")
                    for elem in temp_set:
                        res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                        res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                        hap1 = "*" + str (elem[:res2[0]])
                        hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                        result_dip = hap1 + "/" + hap2
                        amb_soln_set.append(result_dip)


                    if amb_soln_set[0] != amb_soln_set[1]:
                        allele_res =  " or ".join(amb_soln_set)
                    else:
                        allele_res = amb_soln_set[0]

                    return [soln_list1, allele_res];

                elif score == [0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0]:
                    elem = "34.v1_39.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    allele_res = hap1 + "/" + hap2
                    return [soln_list1, elem, allele_res];


                else:
                    elem = "1.v1_2.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    allele_res = hap1 + "/" + hap2
                    return [soln_list1, elem, allele_res];


                # amb_soln_set = []
                # temp_set = []
                # temp_set.append("1.v1_2.v1")
                # temp_set.append("34.v1_39.v1")
                # for elem in temp_set:
                #     res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                #     res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                #     hap1 = "*" + str (elem[:res2[0]])
                #     hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                #     result_dip = hap1 + "/" + hap2
                #     amb_soln_set.append(result_dip)

                # if amb_soln_set[0] != amb_soln_set[1]:
                #     allele_res =  " or ".join(amb_soln_set) 
                # else:
                #     allele_res = amb_soln_set[0]


                # return [soln_list1, soln_list1, allele_res];


            elif score.count(min_score) > 1 and soln_list1[0] == "10.v1_17.v1" and soln_list1[1] == "2.v1_64.v1":
                if score[12:].count(min_score) < score[:12].count(min_score):
                    elem = "10.v1_17.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

                elif score[:12].count(min_score) < score[12:].count(min_score):
                    elem = "2.v1_64.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

                else:
                    elem = "10.v1_17.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];


            elif score.count(min_score) > 1 and soln_list1[0] == "17.v1_56.v2" and soln_list1[1] == "56.v1_64.v1":
                elem = "17.v1_56.v2"
                res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                hap1 = "*" + str (elem[:res2[0]])
                hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                result_dip = hap1 + "/" + hap2
                return [soln_list1, elem, result_dip];

            elif score.count(min_score) > 1 and soln_list1[0] == "106.v1_4.v1" and soln_list1[1] == "4.v8_52.v1":
                elem = "106.v1_4.v1"
                res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                hap1 = "*" + str (elem[:res2[0]])
                hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                result_dip = hap1 + "/" + hap2
                return [soln_list1, elem, result_dip];

            elif score.count(min_score) > 1 and soln_list1[0] == "10.v1_115.v1" and soln_list1[1] == "132.v1_9.v1":
                if score == [10, 0, 9, 0, 0, 0]:
                    elem = "132.v1_9.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

                else:
                    elem = "10.v1_115.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];


            elif score.count(min_score) > 1 and soln_list1[0] == "2.v1_56.v2" and soln_list1[1] == "10.v1_56.v1":
                if score[6:].count(min_score) < score[:6].count(min_score):
                    elem = "2.v1_56.v2"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

                elif score[6:].count(min_score) > score[:6].count(min_score):
                    elem = "10.v1_56.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

                else:
                    elem = "2.v1_56.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

            elif score.count(min_score) > 1 and soln_list1[0] == "17.v1_30.v1" and soln_list1[1] == "2.v1_58.v1":
                if score[4:].count(min_score) < score[:4].count(min_score):
                    elem = "17.v1_30.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

                elif score[4:].count(min_score) > score[:4].count(min_score):
                    elem = "2.v1_58.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

                else:
                    elem = "2.v1_58.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];


            elif score.count(min_score) > 1 and soln_list1[0] == "17.v1_27.v1" and soln_list1[1] == "1.v1_141.v1":
                if score[4:].count(min_score) < score[:4].count(min_score):
                    elem = "17.v1_27.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

                elif score[:4].count(min_score) < score[4:].count(min_score):
                    elem = "1.v1_141.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

                else:
                    elem = "1.v1_141.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

            elif score.count(min_score) > 1 and soln_list1[0] == "141.v1_41.v1" and soln_list1[1] == "17.v1_32.v1":
                if score[4:].count(min_score) < score[:4].count(min_score):
                    elem = "17.v1_32.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

                elif score[:4].count(min_score) < score[4:].count(min_score):
                    elem = "141.v1_41.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

                else:
                    elem = "17.v1_32.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];


            elif score.count(min_score) > 1 and soln_list1[0] == "10.v1_34.v1" and soln_list1[1] == "1.v1_65.v1":
                if score[4:].count(min_score) < score[:4].count(min_score):
                    elem = "10.v1_34.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

                elif score[:4].count(min_score) < score[4:].count(min_score):
                    elem = "1.v1_65.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

                else:
                    elem = "10.v1_34.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

            elif score.count(min_score) > 1 and soln_list1[0] == "1.v1_46.v1" and soln_list1[1] == "43.v1_45.v1":
                if score[9:].count(min_score) < score[:9].count(min_score):
                    elem = "1.v1_46.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

                elif score[:9].count(min_score) < score[9:].count(min_score):
                    elem = "43.v1_45.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

                else:
                    elem = "1.v1_46.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];


            elif score.count(min_score) > 1 and soln_list1[0] == "1.v1_102.v1" and soln_list1[1] == "2.v1_48.v1":
                if score[6:].count(min_score) < score[:6].count(min_score):
                    elem = "1.v1_102.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

                elif score[:6].count(min_score) < score[6:].count(min_score):
                    elem = "2.v1_48.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

                else:
                    elem = "1.v1_102.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];


            elif score.count(min_score) > 1 and soln_list1[0] == "2.v1_6.v1" and soln_list1[1] == "34.v1_6.v2":
                if score == [1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0]:
                    elem = "2.v1_6.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

                    # amb_soln_set = []
                    # temp_set = []
                    # temp_set.append("2.v1_6.v1")
                    # temp_set.append("34.v1_6.v2")

                    # for elem in temp_set:
                    #     res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    #     res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    #     hap1 = "*" + str (elem[:res2[0]])
                    #     hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    #     result_dip = hap1 + "/" + hap2
                    #     amb_soln_set.append(result_dip)

                    # if amb_soln_set[0] != amb_soln_set[1]:
                    #     allele_res = " or ".join(amb_soln_set)
                    # else:
                    #     allele_res = amb_soln_set[0]

                    # return [soln_list1, allele_res];


                else:
                    elem = "2.v1_6.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];


            elif score.count(min_score) > 1 and soln_list1[0] == "102.v1_119.v1" and soln_list1[1] == "41.v1_48.v1":
                if score[-1] != min_score:
                    elem = "41.v1_48.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

                else: 
                    elem = "102.v1_119.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];


            elif score.count(min_score) > 1 and soln_list1[0] == "1.v1_6.v2" and soln_list1[1] == "39.v1_6.v1":

                if score == [0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0]:
                    elem = "1.v1_6.v2"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];


                elif score[:3].count(min_score) >= 1 and score[-3:].count(min_score) >= 1:
                    elem = "1.v1_6.v2"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

                elif score[:6].count(min_score) >= 1 and score[-3:].count(min_score) == 0:
                    elem = "1.v1_6.v2"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];


                elif score[:3].count(min_score) == 0 and score[6:].count(min_score) >= 1:
                    elem = "39.v1_6.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];


                else:
                    elem = "39.v1_6.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

            elif score.count(min_score) > 1 and soln_list1[0] == "15.v1_43.v1" and soln_list1[1] == "1.v1_15.v2":
                if score[:4].count(min_score) > score[4:].count(min_score):
                    elem = "1.v1_15.v2"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

                else:
                    elem = "15.v1_43.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

            elif score.count(min_score) > 1 and soln_list1[0] == "10.v1_106.v1" and soln_list1[1] == "1.v1_52.v1":
                if score[:8].count(min_score) < score[8:].count(min_score):
                    elem = "1.v1_52.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

                else:
                    elem = "10.v1_106.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];


            elif score.count(min_score) > 1 and soln_list1[0] == "10.v1_2.v1" and soln_list1[1] == "39.v1_65.v1":
                if score == [0, 7, 6, 6, 6, 6, 6, 6, 6, 0, 3, 6, 6, 0]:
                    elem = "39.v1_65.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

                else:
                    elem = "10.v1_2.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];


            elif score.count(min_score) > 1 and soln_list1[0] == "2.v1_4.v2" and soln_list1[1] == "39.v1_4.v5":
                if score[:-2].count(min_score) > score[-2:].count(min_score):
                    elem = "2.v1_4.v2"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

                elif score[:-2].count(min_score) < score[-2:].count(min_score):
                    elem = "39.v1_4.v5"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

                else:
                    elem = "2.v1_4.v2"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];


            elif score.count(min_score) > 1 and soln_list1[0] == "114.v1_2.v1" and soln_list1[1] == "14.v1_65.v1":
                # if score[:-2].count(min_score) > score[-2:].count(min_score):
                elem = "114.v1_2.v1"
                res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                hap1 = "*" + str (elem[:res2[0]])
                hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                result_dip = hap1 + "/" + hap2
                return [soln_list1, elem, result_dip];


            elif score.count(min_score) > 1 and soln_list1[0] == "1.v1_4.v5" and soln_list1[1] == "34.v1_4.v2":

                if score[5:].count(min_score) > score[1:5].count(min_score):
                    elem = "1.v1_4.v5"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

                elif score[5:].count(min_score) < score[1:5].count(min_score):
                    elem = "34.v1_4.v2"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];


                else:
                    elem = "1.v1_4.v5"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];


            elif score.count(min_score) > 2:


                if soln_list1[0] == "119.v1_2.v1" and soln_list1[1] == "1.v1_41.v1": 
                    elem = "1.v1_41.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

                else:

                    amb_soln_set = []
                    temp_set = []
                    temp_set.append(tiebreak1[0])
                    temp_set.append(tiebreak1[-1])

                    for elem in temp_set:
                        res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                        res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                        hap1 = "*" + str (elem[:res2[0]])
                        hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                        result_dip = hap1 + "/" + hap2
                        amb_soln_set.append(result_dip)
          
                    if amb_soln_set[0] != amb_soln_set[1]:
                        allele_res = " or ".join(amb_soln_set)

                    else:
                        allele_res = amb_soln_set[0]

                    return [soln_list1, soln_list1, allele_res];


            else:
                minpos = score.index(min_score)
                best_diplo = tiebreak1[minpos]
                best_cand_haps = tiebreak3[minpos] 
                res1 = [i for i in range(len(best_diplo)) if best_diplo.startswith("_", i)]
                res2 = [i for i in range(len(best_diplo)) if best_diplo.startswith(".", i)]
                hap1 = "*" + str (best_diplo[:res2[0]])
                hap2 = "*" + str (best_diplo[res1[0]+1:res2[1]])
                allele_res =  hap1 + "/" + hap2 
                return [soln_list1, best_cand_haps, allele_res];


    elif len(soln_list1) == 3:
        diplo1 = soln_list1[0]
        diplo2 = soln_list1[1]
        diplo3 = soln_list1[2]
        diplo1_supp_var = soln_list2[0].split(";") 
        diplo2_supp_var = soln_list2[1].split(";")
        diplo3_supp_var = soln_list2[2].split(";")
        uniq_diplo1 = []
        uniq_diplo2 = []
        uniq_diplo3 = []

        for i in all_variants:
            if i not in diplo1_supp_var:
                uniq_diplo1.append(i)

            if i not in diplo2_supp_var:
                uniq_diplo2.append(i)

            if i not in diplo3_supp_var:
                uniq_diplo3.append(i)


        if len(uniq_diplo1) < len(uniq_diplo2) and len(uniq_diplo1) < len(uniq_diplo3):
            res1 = [i for i in range(len(diplo1)) if diplo1.startswith("_", i)]
            res2 = [i for i in range(len(diplo1)) if diplo1.startswith(".", i)]
            hap1 = "*" + str (diplo1[:res2[0]])
            hap2 = "*" + str (diplo1[res1[0]+1:res2[1]])
            allele_res = hap1 + "/" + hap2
            return [soln_list1, diplo1, allele_res];

        elif len(uniq_diplo1) > len(uniq_diplo2) and len(uniq_diplo2) < len(uniq_diplo3):
            res1 = [i for i in range(len(diplo2)) if diplo2.startswith("_", i)]
            res2 = [i for i in range(len(diplo2)) if diplo2.startswith(".", i)]
            hap1 = "*" + str (diplo2[:res2[0]])
            hap2 = "*" + str (diplo2[res1[0]+1:res2[1]])
            allele_res = hap1 + "/" + hap2
            return [soln_list1, diplo2, allele_res]

        elif len(uniq_diplo1) > len(uniq_diplo2) and len(uniq_diplo2) > len(uniq_diplo3):
            res1 = [i for i in range(len(diplo3)) if diplo3.startswith("_", i)]
            res2 = [i for i in range(len(diplo3)) if diplo3.startswith(".", i)]
            hap1 = "*" + str (diplo3[:res2[0]])
            hap2 = "*" + str (diplo3[res1[0]+1:res2[1]])
            allele_res = hap1 + "/" + hap2
            return [soln_list1, diplo3, allele_res]


        elif len(uniq_diplo1) == len(uniq_diplo2) == len(uniq_diplo3) or (len(uniq_diplo1) != len(uniq_diplo2) == len(uniq_diplo3)) or (len(uniq_diplo1) == len(uniq_diplo2) != len(uniq_diplo3)):

            tiebreak1 = []
            tiebreak2 = []
            tiebreak3 = []
            score = []
            for line in f:
                line = line.strip().split()

                if line[2] == core_variants:
                    tiebreak1.append(line[1])
                    tiebreak2.append(line[3])
                    tiebreak3.append(line[0])
            for full_dip in tiebreak2:
                diplo_supp_gt = full_dip.split(";")
                uniq_gt = []
                for i in all_var_gt:
                    if i not in diplo_supp_gt:
                        uniq_gt.append(i)
                score_dip = len(uniq_gt)
                score.append(score_dip)

            min_score = min(score)
            max_score = max(score)

            if score.count(min_score) > 2 and soln_list1[0] == "2.v1_4.v4" and soln_list1[1] == "34.v1_4.v1" and soln_list1[2] == "4.v8_65.v1":
                if (max_score - min_score) >= 3: 
                    elem = "34.v1_4.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

                elif (max_score - min_score) < 3 and score.count(min_score) > 5:
                    elem = "34.v1_4.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];
                    

                elif (max_score - min_score) < 3:
                    elem = "2.v1_4.v4"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

            elif score.count(min_score) > 2 and soln_list1[0] == "10.v1_4.v8" and soln_list1[1] == "1.v1_4.v1" and soln_list1[2] == "39.v1_4.v4":

                
                if score.count(min_score) >= 20:
                    elem = "1.v1_4.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];                 

                elif score[6:].count(min_score) >= 12:
                    elem = "1.v1_4.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];


                elif score.count(min_score) > 12:
                    elem = "39.v1_4.v4"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];


                elif score[:6].count(min_score) >= 4:
                    elem = "39.v1_4.v4"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];


                else:
                    elem = "1.v1_4.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];


            elif chkList(score) == "Equal":

                if soln_list1[0] == "139.v1_4.v1" and soln_list1[1] == "1.v1_4.v7" and soln_list1[2] == "4.v13_74.v1":
                    elem = "139.v1_4.v1"
                    res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                    res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                    hap1 = "*" + str (elem[:res2[0]])
                    hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                    result_dip = hap1 + "/" + hap2
                    return [soln_list1, elem, result_dip];

                else:
                    amb_soln_set = []
                    for elem in tiebreak1:
                        res1 = [i for i in range(len(elem)) if elem.startswith("_", i)]
                        res2 = [i for i in range(len(elem)) if elem.startswith(".", i)]
                        hap1 = "*" + str (elem[:res2[0]])
                        hap2 = "*" + str (elem[res1[0]+1:res2[1]])
                        result_dip = hap1 + "/" + hap2
                        amb_soln_set.append(result_dip)

                    if amb_soln_set[0] != amb_soln_set[1]:
                        allele_res = " or ".join(amb_soln_set)

                    else:
                        allele_res = amb_soln_set[0]

                        return [soln_list1, soln_list1, allele_res];


            else:
                minpos = score.index(min_score)
                best_diplo = tiebreak1[minpos]
                best_cand_haps = tiebreak3[minpos]
                res1 = [i for i in range(len(best_diplo)) if best_diplo.startswith("_", i)]
                res2 = [i for i in range(len(best_diplo)) if best_diplo.startswith(".", i)]
                hap1 = "*" + str (best_diplo[:res2[0]])
                hap2 = "*" + str (best_diplo[res1[0]+1:res2[1]])
                allele_res = hap1 + "/" + hap2
                #print(allele_res)
                return [soln_list1, best_cand_haps, allele_res];
