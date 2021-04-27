#!/usr/bin/env python3

import os
import sys

infile = sys.argv[1]

f = open(infile, "r")

for line in f:
    line = line.strip().split()
    ad = line[-1].split(",")
    abhet = round((int(ad[1])/(int(ad[0]) + int(ad[1]))), 4)
    line[-1] = str(abhet)
    line.append('-1')
    print("\t".join(line))
