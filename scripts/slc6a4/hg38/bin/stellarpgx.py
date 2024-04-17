#!/usr/bin/env python3

import os
import sys
import subprocess
import math
from sv_modules import *


print("--------------------------------------------\n")

print("SLC6A4 Variant Analysis with StellarPGx\n")

print("--------------------------------------------\n")


infile = sys.argv[1]
cov_file = sys.argv[2]


def get_core_variants(infile, cn):
    core_vars = []
    for line in open(infile, "r"):
        line = line.strip()
        core_vars.append(line)
    core_vars = ";".join(sorted(core_vars))

    if int(cn) == 1:
        core_vars = core_vars.replace("~0/1", "~1/1")

    if os.stat(infile).st_size == 0:
        core_vars = "No core SNVs detected; haplotypes equivalent to GRCh38 content"

    return core_vars



cn = get_total_CN(cov_file)[0]
av_ctrl_cov = get_total_CN(cov_file)[2]
exon_cov = get_total_CN(cov_file)[3]

print("Initially computed Copy Number = {}".format(cn))

for i in range(1, len(exon_cov)):
    
    if exon_cov[i-1]/av_ctrl_cov < 0.45 :
        print ('Check exon {} for potential deletion if using high coverage WGS'.format(str(i)))
    else:
        pass

supp_core_vars = get_core_variants(infile, cn)

print("\nSample core variants:")
print(supp_core_vars)


