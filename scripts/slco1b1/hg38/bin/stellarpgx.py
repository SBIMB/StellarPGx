#!/usr/bin/env python3

import os
import sys
import subprocess
from snv_def_modules import *
from sv_modules import *
from bkg_modules import *


print("--------------------------------------------\n")

print("SLCO1B1 Star Allele Calling with StellarPGx\n")

print("--------------------------------------------\n")



database = sys.argv[1]
infile = sys.argv[2]
infile_full = sys.argv[3]
infile_full_gt = sys.argv[4]
infile_spec = sys.argv[5]
sv_del = sys.argv[6]
sv_dup = sys.argv[7]
cov_file = sys.argv[8]
hap_dbs = sys.argv[9]
act_score = sys.argv[10]


cn = get_total_CN(cov_file)[0]

print("Initially computed CN = {}".format(cn))

supp_core_vars = get_core_variants(infile, cn)

print("\nSample core variants:")
print(supp_core_vars)


snv_def_calls = cand_snv_allele_calling(database, infile, infile_full, infile_full_gt, infile_spec, cn)

if snv_def_calls == None:

    bac_alleles = get_backgroud_alleles(database, supp_core_vars)

    if bac_alleles == None:
        print("\nResult:")
        print("Possible novel allele or suballele present: interpret with caution")

    else:
        print("\nCandidate alleles:")
        print("[" + bac_alleles[-1] + "]")

        print("\nResult:")
        print("Possible novel allele or suballele present: interpret with caution; experimental validation and expert review through PharmVar is recommended")
        print("\nLikely background alleles:")
        print("[" + bac_alleles[0] + "]")

    sys.exit()



best_diplos = snv_def_calls[0]

print("\nCandidate alleles:")
print(best_diplos)


snv_def_alleles = snv_def_calls[-1]

if "or" in snv_def_alleles:
    pass
else:
    snv_cand_alleles = snv_def_calls[1]


dip_variants = get_all_vars_gt(infile_full_gt)


print("\nResult:")


av_cov = get_total_CN(cov_file)[1]
cov_e1_int4 = get_total_CN(cov_file)[3]
cov_int4_e9 = get_total_CN(cov_file)[4]


gene_alleles = ""


if snv_def_alleles != '*1/*1' and cn != '0':
    in_list = dup_test_init(sv_dup, av_cov)


# if cn == '2' and snv_def_alleles == '*1/*1':


if cn == '2':

    if 'or' in snv_def_alleles:        
        print (snv_def_alleles)

    else:
        snv_def_alleles = snv_def_alleles.split("/")

        if snv_def_alleles[0] == '*1' or snv_def_alleles[1] == '*1':
            ind_star2 = snv_def_alleles.index('*1')
            ind_other = 1 - ind_star2

            test_49 = hybrid_49_test1(cov_e1_int2, cov_int7_e16)
            test_new_hyb = hybrid_new_test1(cov_e1_int2, cov_int7_e16)


            if test_49 == 'norm_var' and test_new_hyb == 'norm_var':
                gene_alleles = "/".join(snv_def_alleles)
                print(gene_alleles)

            elif test_49 == 'hyb_49':
                gene_alleles = snv_def_alleles[ind_other] + "/" + "*49"
                print(gene_alleles)

            elif test_49 == 'hyb_49_2' and snv_def_alleles == "*1/*1":
                gene_alleles = "*49/*49"
                print(gene_alleles)


        else:
            gene_alleles = "/".join(snv_def_alleles)
            print(gene_alleles)



elif cn == '0':
    del_confirm = del_test(sv_del)
    if del_confirm == '*(full_gene_del)/*(full_gene_del)':
        gene_alleles = '*48/*48'
        print (gene_alleles)
        
    elif del_confirm == '*(full_gene_del)':
        gene_alleles = '*48' + "/" + "*other"
        print(gene_alleles)

    else:
        gene_alleles = "*48/*48"
        print(gene_alleles)


elif cn == '1':
    del_confirm = del_test(sv_del)
 
    if "or" in snv_def_alleles and del_confirm == 'None':
        print (snv_def_alleles + "\t" + "Possible SLCO1B1 gene deletion present")

    elif "or" not in snv_def_alleles and del_confirm == 'None':
        snv_def_alleles = snv_def_alleles.split("/")
        snv_cand_alleles = "".join(snv_cand_alleles)
        snv_cand_alleles = snv_cand_alleles.split("_")

        if snv_def_alleles[0] == snv_def_alleles[1]:
            gene_alleles = snv_def_alleles[0] + "/" + "*48"
            print(gene_alleles)

        elif snv_def_alleles[0] != snv_def_alleles[1]:
            samp_allele1 = del_adv_test(hap_dbs, snv_cand_alleles[0], snv_cand_alleles[1], snv_def_alleles[0], snv_def_alleles[1], supp_core_vars)
            
            gene_alleles = samp_allele1 + "/" + "*48"
            print(gene_alleles)

    else:
        snv_def_alleles = snv_def_alleles.split("/")
        snv_cand_alleles = "".join(snv_cand_alleles)
        snv_cand_alleles = snv_cand_alleles.split("_")

        if snv_def_alleles[0] == snv_def_alleles[1]:
        
            if del_confirm == "*(full_gene_del)/*(full_gene_del)":
                del_confirm = "*48"
            gene_alleles = del_confirm + "/" + snv_def_alleles[0]
            print(gene_alleles)

        elif snv_def_alleles[0] != snv_def_alleles[1]:
            samp_allele1 = del_adv_test(hap_dbs, snv_cand_alleles[0], snv_cand_alleles[1], snv_def_alleles[0], snv_def_alleles[1], supp_core_vars)
    
            if del_confirm == "*(full_gene_del)/*(full_gene_del)":
                del_confirm = "*48"
            gene_alleles = del_confirm + "/" + samp_allele1
            print(gene_alleles)



elif (int(cn) == 3 or int(cn) == 4) and snv_def_alleles != None:

    # in_list = dup_test_init(sv_dup, av_cov)
    # print (snv_def_alleles)
    # print (snv_cand_alleles)
    orig = snv_def_alleles

    if "or" in snv_def_alleles:
        print (snv_def_alleles + "\t" + "Duplication present")

    else:
        snv_def_alleles = snv_def_alleles.split("/")
        snv_cand_alleles = "".join(snv_cand_alleles)
        snv_cand_alleles = snv_cand_alleles.split("_")

        if snv_def_alleles[0] != snv_def_alleles[1]:
            # print("\n" + dup_test(sv_dup, hap_dbs, snv_def_alleles[0], snv_def_alleles[1], cn))
            # print (snv_cand_alleles)
            # print ("\n")
            phased_dup = dup_test_cn_3_4(sv_dup, hap_dbs, snv_cand_alleles[0], snv_cand_alleles[1], snv_def_alleles[0], snv_def_alleles[1], cn, av_cov, in_list)

            phased_dup1 = phased_dup.split("/")


        elif snv_def_alleles[0] == snv_def_alleles[1]:
            
            rt_2 = int(cn) - 1
            
            phased_dup = (snv_def_alleles[0] + "/" + snv_def_alleles[1] + "x" + str(rt_2))


        gene_alleles = phased_dup
        
        print(gene_alleles)


elif int(cn) > 4 and snv_def_alleles != None:

    if "or" in snv_def_alleles:
        print (snv_def_alleles + "\t" + "Duplication present")

    else:
        snv_def_alleles = snv_def_alleles.split("/")
        snv_cand_alleles = "".join(snv_cand_alleles)
        snv_cand_alleles = snv_cand_alleles.split("_")

        if snv_def_alleles[0] != snv_def_alleles[1]:

            phased_dup = dup_test_cn_n(sv_dup, hap_dbs, snv_cand_alleles[0], snv_cand_alleles[1], snv_def_alleles[0], snv_def_alleles[1], cn, av_cov, in_list)
        elif snv_def_alleles[0] == snv_def_alleles[1]:
            rt_2 = int(cn) - 1
            phased_dup = (snv_def_alleles[0] + "/" + snv_def_alleles[1] + "x" + str(rt_2))

        gene_alleles = phased_dup
        print(phased_dup)



elif int(cn) > 2 and snv_def_alleles == None:
    
    print("Possible rare SLCO1B1 hybrid present")
