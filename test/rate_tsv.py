#!/usr/bin/python

import os
import sys
import math
import csv
import argparse
import numpy as np

posExp = -1
posRun = -1
posSam = -1
posTar = -1
posTD0 = -1
posPCReff = -1
posPCReffsd = -1
posNcopy = -1
expData = {}
runs = []
exps = []
inData = []

parser = argparse.ArgumentParser(description='Merge TSV Files.')
parser.add_argument('-f', '--file', metavar="data.tsv", dest='upfile', help='tsv file to process')
parser.add_argument('-o', '--out', metavar="out.tsv", dest='outfile', help='tsv file for the results')
parser.add_argument('-n', '--n_out', action='store_true', help='print the number of elements')
parser.add_argument('-pm', '--plusminus', action='store_true', help='print the sd as plus/minus symbol')
parser.add_argument('-ph', '--print_header', action='store_true', help='print the header of tsv file')

args = parser.parse_args()
if not args.upfile:
    print("Error no file provided using -f")
    sys.exit(0)

with open(args.upfile, newline='') as tfile:  # add encoding='utf-8' ?
    inData = list(csv.reader(tfile, delimiter="\t"))

if args.print_header:
    for col in range(0, len(inData[0])):
        print(str(col  + 1) + ": " + inData[0][col])

for col in range(0, len(inData[0])):
    if inData[0][col] == "Experiment":
        posExp = col
    if inData[0][col] == "Run":
        posRun = col
    if inData[0][col] == "sample":
        posSam = col
    if inData[0][col] == "target":
        posTar = col
    if inData[0][col] == "TD0":
        posTD0 = col
    if inData[0][col] == "PCR eff":
        posPCReff = col
    if inData[0][col] == "STD PCR eff":
        posPCReffsd = col
    if inData[0][col] == "Ncopy":
        posNcopy = col

if posExp > -1 and posRun > -1:
    for row in range(1, len(inData)):
        if inData[row][posExp] not in expData:
            expData[inData[row][posExp]] = {}
        expData[inData[row][posExp]][inData[row][posRun]] = {}

ww = open(args.outfile, "w")
for exp in expData:
    for run in expData[exp]:
        exps.append(exp)
        runs.append(run)

count = 0
keep_run = True
while keep_run:
    data = {}
    for row in range(1, len(inData)):
        if len(exps) > 0 and inData[row][posExp] != exps[count]:
            continue
        if len(runs) > 0 and inData[row][posRun] != runs[count]:
            continue
        if inData[row][posTD0] == "TD0":
            continue

        if inData[row][posTar] not in data:
            data[inData[row][posTar]] = {}
        if inData[row][posSam] not in data[inData[row][posTar]]:
            data[inData[row][posTar]][inData[row][posSam]] = {}
            data[inData[row][posTar]][inData[row][posSam]]["TD0"] = {"raw": []}
            data[inData[row][posTar]][inData[row][posSam]]["PCR eff"] = {"raw": []}
            data[inData[row][posTar]][inData[row][posSam]]["STD PCR eff"] = {"raw": []}
            data[inData[row][posTar]][inData[row][posSam]]["Ncopy"] = {"raw": []}

        if float(inData[row][posTD0]) > 1.0:
            data[inData[row][posTar]][inData[row][posSam]]["TD0"]["raw"].append(float(inData[row][posTD0]))
        if 1.1 < float(inData[row][posPCReff]) < 10.0:
            data[inData[row][posTar]][inData[row][posSam]]["PCR eff"]["raw"].append(float(inData[row][posPCReff]))
            data[inData[row][posTar]][inData[row][posSam]]["STD PCR eff"]["raw"].append(float(inData[row][posPCReffsd]))
        if float(inData[row][posNcopy]) > 1.0:
            data[inData[row][posTar]][inData[row][posSam]]["Ncopy"]["raw"].append(float(inData[row][posNcopy]))


    for cc1 in data:
        for cc2 in data[cc1]:
            for cc3 in data[cc1][cc2]:
                data[cc1][cc2][cc3]["n"] = len(data[cc1][cc2][cc3]["raw"])
                if data[cc1][cc2][cc3]["n"] > 0:
                    data[cc1][cc2][cc3]["mean"] = np.mean(data[cc1][cc2][cc3]["raw"])
                    data[cc1][cc2][cc3]["std"] = np.std(data[cc1][cc2][cc3]["raw"])
                else:
                    data[cc1][cc2][cc3]["mean"] = -1.0
                    data[cc1][cc2][cc3]["std"] = -1.0

    std_sep = '\t'
    if args.plusminus:
        std_sep = ' \xB1'

    if count > 0:
        if args.n_out:
            ww.write("\t\t")
        if std_sep == '\t':
            ww.write("\t\t\t")
        ww.write("\t\t\t\t\n")

    if args.n_out:
        ww.write("Target\tSample\tn\tTD0" + std_sep + "SD\tPCR Eff" + std_sep + "SD\tn\tNcopy" + std_sep + "SD\n")
    else:
        ww.write("Target\tSample\tTD0" + std_sep + "SD\tPCR Eff" + std_sep + "SD\tNcopy" + std_sep + "SD\n")
    for cc1 in data:
        for cc2 in data[cc1]:
            ww.write(cc1 + "\t" + cc2)

            if args.n_out:
                ww.write("\t" + str(data[cc1][cc2]["TD0"]["n"]))
            ww.write("\t" + "{:.2f}".format(data[cc1][cc2]["TD0"]["mean"]))
            ww.write(std_sep + "{:.2f}".format(data[cc1][cc2]["TD0"]["std"]))

            ww.write("\t" + "{:.4f}".format(data[cc1][cc2]["PCR eff"]["mean"]))
            ww.write(std_sep + "{:.4f}".format(data[cc1][cc2]["STD PCR eff"]["mean"]))

            if args.n_out:
                ww.write("\t" + str(data[cc1][cc2]["Ncopy"]["n"]))
            ww.write("\t" + "{:.1f}".format(data[cc1][cc2]["Ncopy"]["mean"]))
            ww.write(std_sep + "{:.1f}".format(data[cc1][cc2]["Ncopy"]["std"]))



            ww.write("\n")
    count += 1
    if count >= len(exps):
        keep_run = False


ww.close()
