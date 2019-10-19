#!/usr/bin/python

from __future__ import division
from __future__ import print_function

import sys
import re
import os


if sys.argv[1] in ["-h", "--help"] or len(sys.argv) != 7:
    print("Usage:\ndiff_table.py [File A] [File B] [Name of comparison] [Max number of returned lines] [Remove tailing tabs] [Show column differences]\n")
else:
    countOutLines = 1
    maxLines = int(sys.argv[4])

    # Open the files
    with open(sys.argv[1], 'r') as txtfileA:
        aFile = txtfileA.read()
    with open(sys.argv[2], 'r') as txtfileB:
        bFile = txtfileB.read()

    # Fix windows newlines
    aFile = re.sub(r"\r\n", "\n", aFile)
    bFile = re.sub(r"\r\n", "\n", bFile)

    # Fix different NaN styles
    aFile = re.sub(r"NaN", "nan", aFile)
    bFile = re.sub(r"NaN", "nan", bFile)

    # Remove tailing tabs
    if sys.argv[5] in ["y", "Y", "yes", "Yes"]:
        aFile = re.sub(r"\t+\n", "\n", aFile)
        bFile = re.sub(r"\t+\n", "\n", bFile)

    # Show column differences
    showColDiff = False
    if sys.argv[6] in ["y", "Y", "yes", "Yes"]:
        showColDiff = True

    if aFile == bFile:
        print(sys.argv[3] + ": [OK]")
    else:
        print(sys.argv[3] + ": [Failed]")
        arrAA = aFile.split('\n')
        arrBB = bFile.split('\n')
        minLen = len(arrAA)
        if len(arrAA) != len(arrBB):
            countOutLines += 1
            print("  Different number of rows: " + str(len(arrAA)) + " != " + str(len(arrBB)))
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
                    if showColDiff:
                        for k in range(0, minLineLen):
                            if lineAA[k] != lineBB[k]:
                                countOutLines += 1
                                if countOutLines < maxLines:
                                    print("  Difference in " + str(i) + "-" + str(k) + ": " + lineAA[k] + " != " + lineBB[k])

        if countOutLines >= maxLines:
            print("    More differences were not shown due to line limit!")




