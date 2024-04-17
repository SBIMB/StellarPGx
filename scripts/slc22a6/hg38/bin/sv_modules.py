#!/usr/bin/env python3

import os
import sys
import math
        

def get_total_CN(cov_file):

    all_reg =[]
    for line in open(cov_file, "r"):
        line = line.strip().split()
        all_reg.append(line)

    av_slc22a6_cov = float(all_reg[0][3])/(float(all_reg[0][2]) - float(all_reg[0][1]))
    av_vdr_cov = float(all_reg[1][3])/(float(all_reg[1][2]) - float(all_reg[1][1]))
    av_egfr_cov = float(all_reg[2][3])/(float(all_reg[2][2]) - float(all_reg[2][1]))

    exon_cov_list = []
    all_reg = all_reg[3:]

    a = list(range(1, len(all_reg)))

    for i in a:
        exon_cov = 'av_e' + str(i)  
        exon_cov = float(all_reg[i-1][3])/(float(all_reg[i-1][2]) - float(all_reg[i-1][1]))
        exon_cov_list.append(exon_cov)
        
    av_ctrl_cov = (av_vdr_cov + av_egfr_cov)/2

    comp_av = av_slc22a6_cov/av_ctrl_cov
    temp_cn = 2 * comp_av
    total_cn = round(temp_cn)


    return [str(int(total_cn)), round(av_slc22a6_cov), round(av_ctrl_cov), exon_cov_list]; 

