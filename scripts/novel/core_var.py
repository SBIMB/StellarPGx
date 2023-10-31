#!/usr/bin/env/ python3

import os
import sys

annot_file = sys.argv[1]
gene_symbol = sys.argv[2]
transcript = sys.argv[3]

var_conseq = ["missense", "frameshift", "stop_gained", "splice_donor", "splice_acceptor", "inframe_insertion", "inframe_deletion", "start_lost", "stop_lost"] 

for line in open(annot_file, "r"):
    line = line.strip()

    if line.startswith("#"):
        pass
    
    else:
        for i in var_conseq:
            pattern = i + "|" + gene_symbol + "|" + transcript
            
            if pattern in line:
                print(line)

            else:
                pass
