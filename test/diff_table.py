#!/usr/bin/python

from __future__ import division
from __future__ import print_function

import sys
import os
import re
import datetime
import zipfile
import tempfile
import argparse
import math
import numpy as np


if sys.argv[1] in ["-h", "--help"] or len(sys.argv) != 5:
    print("Usage:\ndiff_table.py [File A] [File B] [Name of comparison] [Max number of returned lines]\n")
else:
    countOutLines = 1
    maxLines = int(sys.argv[4])
    with open(sys.argv[1], 'r') as txtfileA:
        aFile = txtfileA.read()
    with open(sys.argv[2], 'r') as txtfileB:
        bFile = txtfileB.read()

    if aFile == bFile:
        print(sys.argv[3] + ": [OK]")
    else:
        print(sys.argv[3] + ": [Failed]")
        arrAA = aFile.split('\n')
        arrBB = bFile.split('\n')
        minLen = len(arrAA)
        if len(arrAA) != len(arrBB):
            countOutLines += 1
            print("  Different number of rows: " + str(len(arrAA)) + " != "  + str(len(arrBB)))
        if len(arrAA) > len(arrBB):
            minLen = len(arrBB)
        for i in range(0, minLen):
            if arrAA[i] != arrBB[i]:
                countOutLines += 1
                if countOutLines < maxLines:
                    print("  Difference in row: " + str(i))
                    lineAA = arrAA[i].split('\t')
                    lineBB = arrBB[i].split('\t')
                    minLineLen = len(lineAA)
                    if len(lineAA) > len(lineBB):
                        minLineLen = len(lineBB)
                    if len(lineAA) != len(lineBB):
                        countOutLines += 1
                        if countOutLines < maxLines:
                            print("    Different number of columns: " + str(len(lineAA)) + " != " + str(len(lineBB)))
                    for k in range(0, minLineLen):
                        if lineAA[k] != lineBB[k]:
                            countOutLines += 1
                            if countOutLines < maxLines:
                                print("    Difference in " + str(i) + "-" + str(k) + ": " + lineAA[k] + " != " + lineBB[k])




