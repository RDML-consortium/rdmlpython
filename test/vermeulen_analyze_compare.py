#!/usr/bin/python

from __future__ import division
from __future__ import print_function

import os
import json
import sys
import time
import re
import math
import csv
import numpy as np
import scipy.stats as scp

parent_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))

head = ["Standard-Cq", "LinRegPCR (old)", "5PSM", "FPK-PCR", "LRE-Emax", "LRE-E100", "Cy0", "MAK2", "DART", "4PLM", "PCR-Miner", "RDML-Tools"]

startPos = 0
data = []
with open(os.path.join(parent_dir, "experiments/vermeulen/data_vermeulen_methods_results.tsv"), newline='') as tfile:
    table = list(csv.reader(tfile, delimiter="\t"))
    for row in range(startPos * 13 + 1, startPos * 13 + 12):
        data.append([])
        for col in range(1, len(table[row])):
            data[row - startPos * 13 - 1].append(99999999999)
            try:
                vFloat = float(table[row][col])
            except ValueError:
                print("Error " + str(row) + " " + str(col) + " : " + table[row][0] + " " + table[0][col]+ " : " + table[row][col])
            else:
                if math.isfinite(vFloat):
                    data[row - startPos * 13 - 1][col - 1] = vFloat
                else:
                    print("Error finite " + str(row) + " " + str(col) + " : " + table[row][0] + " " + table[0][col]+ " : " +  table[row][col])

for row in range(0, len(data)):
    for col in range(0, len(data[row])):
        if startPos == 0:
            if data[row][col] < 10000.0:
                data[row][col] = 10000.0 / data[row][col]
            else:
                data[row][col] = data[row][col] / 10000.0

with open("temp_vermeulen_report_bias.csv", newline='') as tfile:
    table = list(csv.reader(tfile, delimiter="\t"))
    data.append([])
    for col in range(2, len(table[40])):
        data[11].append(99999999999)
        try:
            vFloat = float(table[40][col])
        except ValueError:
            print("Error: " + table[40][col])
        else:
            if math.isfinite(vFloat):
                if vFloat < 10000.0:
                    vFloat = 10000.0 / vFloat
                else:
                    vFloat = vFloat / 10000.0
                data[11][col - 2] = vFloat
            else:
                print("Error finite: " + table[40][col])

res = scp.rankdata(data, axis=0, method='average')
bias = np.mean(res, axis=1)
bias_mean = np.mean(data, axis=1)
bias_std = np.std(data, axis=1, ddof=1)

startPos = 3
data = []
with open(os.path.join(parent_dir, "experiments/vermeulen/data_vermeulen_methods_results.tsv"), newline='') as tfile:
    table = list(csv.reader(tfile, delimiter="\t"))
    for row in range(startPos * 13 + 1, startPos * 13 + 12):
        data.append([])
        for col in range(1, len(table[row])):
            data[row - startPos * 13 - 1].append(99999999999)
            try:
                vFloat = float(table[row][col])
            except ValueError:
                print("Error " + str(row) + " " + str(col) + " : " + table[row][0] + " " + table[0][col]+ " : " + table[row][col])
            else:
                if math.isfinite(vFloat):
                    data[row - startPos * 13 - 1][col - 1] = vFloat
                else:
                    print("Error finite " + str(row) + " " + str(col) + " : " + table[row][0] + " " + table[0][col]+ " : " +  table[row][col])

with open("temp_vermeulen_report_bias.csv", newline='') as tfile:
    table = list(csv.reader(tfile, delimiter="\t"))
    data.append([])
    for col in range(2, len(table[111])):
        data[11].append(99999999999)
        try:
            vFloat = float(table[111][col])
        except ValueError:
            print("Error: " + table[111][col])
        else:
            if math.isfinite(vFloat):
                data[11][col - 2] = vFloat
            else:
                print("Error finite: " + table[111][col])

res = scp.rankdata(data, axis=0, method='average')
linearity = np.mean(res, axis=1)
linearity_mean = np.mean(data, axis=1)
linearity_std = np.std(data, axis=1, ddof=1)

startPos = 4
data = []
with open(os.path.join(parent_dir, "experiments/vermeulen/data_vermeulen_methods_results.tsv"), newline='') as tfile:
    table = list(csv.reader(tfile, delimiter="\t"))
    for row in range(startPos * 13 + 1, startPos * 13 + 12):
        data.append([])
        for col in range(1, len(table[row])):
            data[row - startPos * 13 - 1].append(99999999999)
            try:
                vFloat = float(table[row][col])
            except ValueError:
                print("Error " + str(row) + " " + str(col) + " : " + table[row][0] + " " + table[0][col]+ " : " + table[row][col])
            else:
                if math.isfinite(vFloat):
                    data[row - startPos * 13 - 1][col - 1] = vFloat
                else:
                    print("Error finite " + str(row) + " " + str(col) + " : " + table[row][0] + " " + table[0][col]+ " : " +  table[row][col])

with open("temp_vermeulen_report_bias.csv", newline='') as tfile:
    table = list(csv.reader(tfile, delimiter="\t"))
    data.append([])
    for col in range(2, len(table[112])):
        data[11].append(99999999999)
        try:
            vFloat = float(table[112][col])
        except ValueError:
            print("Error: " + table[112][col])
        else:
            if math.isfinite(vFloat):
                data[11][col - 2] = vFloat
            else:
                print("Error finite: " + table[112][col])

res = scp.rankdata(data, axis=0, method='average')
precision = np.mean(res, axis=1)
precision_mean = np.mean(data, axis=1)
precision_std = np.std(data, axis=1, ddof=1)

startPos = 5
data = []
with open(os.path.join(parent_dir, "experiments/vermeulen/data_vermeulen_methods_results.tsv"), newline='') as tfile:
    table = list(csv.reader(tfile, delimiter="\t"))
    for row in range(startPos * 13 + 1, startPos * 13 + 12):
        data.append([])
        for col in range(1, len(table[row])):
            data[row - startPos * 13 - 1].append(99999999999)
            try:
                vFloat = float(table[row][col])
            except ValueError:
                print("Error " + str(row) + " " + str(col) + " : " + table[row][0] + " " + table[0][col]+ " : " + table[row][col])
            else:
                if math.isfinite(vFloat):
                    data[row - startPos * 13 - 1][col - 1] = vFloat
                else:
                    print("Error finite " + str(row) + " " + str(col) + " : " + table[row][0] + " " + table[0][col]+ " : " +  table[row][col])

with open("temp_vermeulen_report_difference.csv", newline='') as tfile:
    table = list(csv.reader(tfile, delimiter="\t"))
    data.append([])
    for col in range(3, len(table[97])):
        data[11].append(99999999999)
        try:
            vFloat = float(table[97][col])
        except ValueError:
            print("Error: " + table[97][col])
        else:
            if math.isfinite(vFloat):
                data[11][col - 3] = vFloat
            else:
                print("Error finite: " + table[97][col])

res = scp.rankdata(data, axis=0, method='average')
resolution = np.mean(res, axis=1)
resolution_mean = np.mean(data, axis=1)
resolution_std = np.std(data, axis=1, ddof=1)

final = []
for row in range(0, len(resolution)):
    final.append([bias[row], linearity[row], precision[row], resolution[row]])

res = scp.rankdata(final, axis=0, method='average')
finMean = np.mean(res, axis=1)
resi = scp.rankdata(finMean, axis=0, method='ordinal')
lookUp = {}
for pos in range(0, len(resi)):
    lookUp[resi[pos] - 1] = pos

print("                     Bias         Linearity       Precision       Resolution       Mean rank")
for pos in range(0, len(resi)):
    aa = ''
    bb = ''
    if head[lookUp[pos]] == "RDML-Tools":
        aa = '\033[44m'
        bb = '\033[0m'
    print(aa + head[lookUp[pos]].ljust(17) + "{:5.2f}".format(bias[lookUp[pos]]) + "  (" + "{:2.0f}".format(res[lookUp[pos]][0]) + ")     " 
                                           + "{:5.2f}".format(linearity[lookUp[pos]]) + "  (" + "{:2.0f}".format(res[lookUp[pos]][1]) + ")     " 
                                           + "{:5.2f}".format(precision[lookUp[pos]]) + "  (" + "{:2.0f}".format(res[lookUp[pos]][2]) + ")     " 
                                           + "{:5.2f}".format(resolution[lookUp[pos]]) + "  (" + "{:2.0f}".format(res[lookUp[pos]][3]) + ")        " 
                                           + "{:5.2f}".format(finMean[lookUp[pos]]) + bb)
print("")

final2 = []
for row in range(0, len(resolution)):
    final2.append([bias_mean[row], linearity_mean[row], precision_mean[row], resolution_mean[row]])

res2 = scp.rankdata(final2, axis=0, method='average')
fin2Mean = np.mean(res2, axis=1)
resi2 = scp.rankdata(fin2Mean, axis=0, method='ordinal')
lookUp = {}
for pos in range(0, len(resi2)):
    lookUp[resi2[pos] - 1] = pos

print("                          Bias                     Linearity                   Precision                  Resolution           Mean rank")
for pos in range(0, len(resi)):
    aa = ''
    bb = ''
    if head[lookUp[pos]] == "RDML-Tools":
        aa = '\033[44m'
        bb = '\033[0m'
    print(aa + head[lookUp[pos]].ljust(17) + "{:6.3f}".format(bias_mean[lookUp[pos]]) + " (" + "{:6.3f}".format(bias_std[lookUp[pos]]) + ") (" + "{:2.0f}".format(res2[lookUp[pos]][0]) + ")    " 
                                           + "{:8.5f}".format(linearity_mean[lookUp[pos]]) + " (" + "{:8.5f}".format(linearity_std[lookUp[pos]]) + ") (" + "{:2.0f}".format(res2[lookUp[pos]][1]) + ")    " 
                                           + "{:8.5f}".format(precision_mean[lookUp[pos]]) + " (" + "{:8.5f}".format(precision_std[lookUp[pos]]) + ") (" + "{:2.0f}".format(res2[lookUp[pos]][2]) + ")    " 
                                           + "{:8.5f}".format(resolution_mean[lookUp[pos]]) + " (" + "{:8.5f}".format(resolution_std[lookUp[pos]]) + ") (" + "{:2.0f}".format(res2[lookUp[pos]][3]) + ")       " 
                                           + "{:5.2f}".format(fin2Mean[lookUp[pos]]) + bb)
