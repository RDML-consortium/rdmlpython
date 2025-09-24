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

parent_dir = os.path.abspath(os.path.join(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir), os.pardir))
sys.path.append(parent_dir)

import rdmlpython as rdml


###########################################################################################
###  This script runs the equivalent calculations as the excel of the paper supplement  ###
###########################################################################################

use_results = False  # True - results from paper, False - use output form run_analyze_vermeulen_rdml.py
set_method = 0  # 0 - Standard-Cq, 1 - LinRegPCR, 2 - 5PSM, 3 - FPK-PCR, 4 - LRE-Emax, 5 - LRE-E100, 6 - Cy0, 7 - MAK2, 8 - DART, 9 - 4PLM, 10 - PCR-Miner

out_file_bias = "temp_vermeulen_report_bias.csv"
out_file_var = "temp_vermeulen_report_variation.csv"
out_file_diff = "temp_vermeulen_report_difference.csv"

ww = open(out_file_bias, "w")
wv = open(out_file_var, "w")
wd = open(out_file_diff, "w")

data = {}
noCq = False

rar_id = 0
rar_well = 1
rar_sample = 2
rar_sample_type = 3
rar_sample_nucleotide = 4
rar_tar = 5
rar_tar_type = 6
rar_tar_dye = 7
rar_tar_chemistry = 8
rar_excl = 9
rar_note = 10
rar_baseline = 11
rar_plat_val = 12
rar_plat_base = 13
rar_n_log = 14
rar_n_included = 15
rar_stop_log = 16
rar_indiv_PCR_eff_x = 17
rar_indiv_PCR_eff_y = 18
rar_indiv_PCR_eff = 19
rar_R2 = 20
rar_PCR_eff = 21
rar_PCR_eff_err = 22
rar_TD0_fluor = 23
rar_TD0 = 24
rar_indiv_Ncopy = 25
rar_Ncopy = 26
rar_amplification = 27
rar_baseline_error = 28
rar_plateau = 29
rar_excl_eff = 30
rar_isUsedForMeanPCREff = 31
rar_nAmpli = 32
rar_vol = 33
rar_nCopyFact = 34
rar_dyeConc = 35
rar_primer_len_for = 36
rar_primer_conc_for = 37
rar_primer_len_rev = 38
rar_primer_conc_rev = 39
rar_probe1_len = 40
rar_probe1_conc = 41
rar_probe2_len = 42
rar_probe2_conc = 43
rar_amplicon_len = 44

if use_results:
    with open("data_vermeulen_methods_quant.csv", newline='') as tfile:
        table = list(csv.reader(tfile, delimiter=";"))
        startRow = 17 * (set_method) + 1
        for row in range(startRow + 1, startRow + 16):
            for col in range(3, len(table[row])):
                if len(table[startRow][col]) < 1:
                    continue
                if table[startRow][col] == "ALUsq(Eurogentec)":
                    continue
                if table[startRow][col] not in data:
                    data[table[startRow][col]] = {}
                    data[table[startRow][col]]["raw"] = {}
                if table[row][2] not in data[table[startRow][col]]["raw"]:
                    data[table[startRow][col]]["raw"][table[row][2]] = []
                try:
                    vFloat = float(table[row][col])
                except ValueError:
                    print("Error " + str(row) + " " + str(col) + " : " + table[startRow][col] + " - " + table[row][2] + " = " + table[row][col])
                else:
                    if math.isfinite(vFloat):
                        data[table[startRow][col]]["raw"][table[row][2]].append(vFloat)
                    else:
                        print("Error " + str(row) + " " + str(col) + " : " + table[startRow][col] + " - " + table[row][2] + " = " + table[row][col])

    if set_method not in [4, 5, 7]:
        with open("data_vermeulen_methods_cq.csv", newline='') as tfile:
            table = list(csv.reader(tfile, delimiter=";"))
            startRow = 17 * (set_method) + 1
            for row in range(startRow + 1, startRow + 16):
                for col in range(3, len(table[row])):
                    if len(table[startRow][col]) < 1:
                        continue
                    if table[startRow][col] == "ALUsq(Eurogentec)":
                       continue
                    if table[startRow][col] not in data:
                        data[table[startRow][col]] = {}
                    if "Cq" not in data[table[startRow][col]]:
                        data[table[startRow][col]]["Cq"] = {}
                    if table[row][2] not in data[table[startRow][col]]["Cq"]:
                        data[table[startRow][col]]["Cq"][table[row][2]] = []
                    try:
                        vFloat = float(table[row][col])
                    except ValueError:
                        print("Error Cq " + str(row) + " " + str(col) + " : " + table[startRow][col] + " - " + table[row][2] + " = " + table[row][col])
                    else:
                        if math.isfinite(vFloat):
                            data[table[startRow][col]]["Cq"][table[row][2]].append(vFloat)
                        else:
                            print("Error Cq " + str(row) + " " + str(col) + " : " + table[startRow][col] + " - " + table[row][2] + " = " + table[row][col])
    else:
        noCq = True

else:
    noCq = True
    col = rar_Ncopy + 2
    genes = ["AHCY", "AKR1C1", "ARHGEF7", "BIRC5", "CAMTA1", "CAMTA2", "CD44", "CDCA5", "CDH5", "CDKN3", 
             "CLSTN1", "CPSG3", "DDC", "DPYSL3", "ECEL1", "ELAVL4", "EPB41L3", "EPHA5", "EPN2", "FYN", 
             "GNB1", "HIVEP2", "HMBS", "HPRT1", "IGSF4", "INPP1", "MAP2K4", "MAP7", "MAPT", "MCM2", 
             "MRPL3", "MTSS1", "MYCN(4)", "NHLH2", "NM23A", "NRCAM", "NTRK1", "ODC1", "PAICS", "PDE4DIP", 
             "PIK3R1", "PLAGL1", "PLAT", "PMP22", "PRAME", "PRDM2", "PRKACB", "PRKCZ", "PTN", "PTPRF", 
             "PTPRH", "PTPRN2", "QPCT", "SCG2", "SDHA(1)", "SLC25A5", "SLC6A8", "SNAPC1", "TNFRSF", 
             "TYMS", "UBC(2)", "ULK2", "WSB1"]
    sampes = ["STD_15", "STD_150", "STD_1500", "STD_15000", "STD_150000"]
    with open("temp_vermeulen_resuls.csv", newline='') as tfile:
        table = list(csv.reader(tfile, delimiter="\t"))
        for row in range(1, len(table)):
            if table[row][0] != "biomarker_set":
                continue
            geneName = table[row][rar_tar + 2].split("_")
            if geneName[0] not in genes:
                continue
            if table[row][rar_sample + 2] not in sampes:
                continue
            if geneName[0] not in data:
                    data[geneName[0]] = {}
                    data[geneName[0]]["raw"] = {}
            if table[row][rar_sample + 2] not in data[geneName[0]]["raw"]:
                data[geneName[0]]["raw"][table[row][rar_sample + 2]] = []
            try:
                vFloat = float(table[row][col])
            except ValueError:
                print("Error " + str(geneName[0]) + " " + str(table[row][rar_sample + 2]) + " : " + table[row][col])
            else:
                if math.isfinite(vFloat):
                    data[geneName[0]]["raw"][table[row][rar_sample + 2]].append(vFloat)
                else:
                    print("Error " + str(geneName[0]) + " " + str(table[row][rar_sample + 2]) + " : " + table[row][col])

targets = sorted(list(data.keys()))

# Calculate everything
data = rdml.standardCurveStats(data, noCq)

#############################
###  Write out Bias tabe  ###
#############################

# create an Empty line
emptyLine = "\t\t"
shortEmptyLine = ""
for tar in range(0, len(targets)):
    if tar < len(targets) - 1:
        shortEmptyLine += "\t"
    else:
        shortEmptyLine += "\n"
emptyLine += shortEmptyLine

# Write target names
ww.write("\t\t")
for tar in range(0, len(targets)):
    ww.write(targets[tar])
    if tar < len(targets) - 1:
        ww.write("\t")
    else:
        ww.write("\n")

# Write the input data
for concStr in ["15", "150", "1500", "15000", "150000"]:
    for pos in range(0, 3):
        ww.write("\t" + concStr + "\t")
        for tar in range(0, len(targets)):
            ww.write("{0:0.3e}".format(data[targets[tar]]["raw"]["STD_" + concStr][pos]))
            if tar < len(targets) - 1:
                ww.write("\t")
            else:
                ww.write("\n")

# mean(max) write out
ww.write(emptyLine)
ww.write(emptyLine)
ww.write("\tmean(max)\t")
for tar in range(0, len(targets)):
    if tar < len(targets) - 1:
        ww.write("mean(max)\t")
    else:
        ww.write("mean(max)\n")
ww.write("\t150000\t")
for tar in range(0, len(targets)):
    ww.write("{0:0.3e}".format(data[targets[tar]]["mean_max"]))
    if tar < len(targets) - 1:
        ww.write("\t")
    else:
        ww.write("\n")

# Write the scaled data
ww.write(emptyLine)
ww.write(emptyLine)
for concStr in ["15", "150", "1500", "15000", "150000"]:
    for pos in range(0, 3):
        ww.write("\t" + str(float(concStr) / 150000) + "\t")
        for tar in range(0, len(targets)):
            ww.write("{0:0.6f}".format(data[targets[tar]]["scaled"]["STD_" + concStr][pos]))
            if tar < len(targets) - 1:
                ww.write("\t")
            else:
                ww.write("\n")

# Calculate scaled(min) and write out
ww.write(emptyLine)
ww.write("scaled(min)\t" + str(15 / 150000) + "\t")
for tar in range(0, len(targets)):
    ww.write("{0:0.6f}".format(data[targets[tar]]["scaled_mean_min"]))
    if tar < len(targets) - 1:
        ww.write("\t")
    else:
        ww.write("\n")
# Calculate scaled(max) and write out
ww.write("scaled(max)\t" + str(150000 / 150000) + "\t")
for tar in range(0, len(targets)):
    ww.write("{0:0.6f}".format(data[targets[tar]]["scaled_mean_max"]))
    if tar < len(targets) - 1:
        ww.write("\t")
    else:
        ww.write("\n")
# Calculate Bias as ratio mean and write out
ww.write("ratio\tBias as ratio mean\t")
for tar in range(0, len(targets)):
    ww.write("{0:0.2f}".format(data[targets[tar]]["bias_as_ratio"]))
    if tar < len(targets) - 1:
        ww.write("\t")
    else:
        ww.write("\n")

# Write the log scaled data
ww.write(emptyLine)
ww.write("\ttake log\t" + shortEmptyLine)
for concStr in ["15", "150", "1500", "15000", "150000"]:
    for pos in range(0, 3):
        ww.write("\t" + str(np.log10(float(concStr) / 150000)) + "\t")
        for tar in range(0, len(targets)):
            ww.write("{0:0.6f}".format(np.log10(data[targets[tar]]["scaled"]["STD_" + concStr][pos])))
            if tar < len(targets) - 1:
                ww.write("\t")
            else:
                ww.write("\n")           
ww.write("SSx\t30\t" + shortEmptyLine)
# Write target names
ww.write("\t\t")
for tar in range(0, len(targets)):
    ww.write(targets[tar])
    if tar < len(targets) - 1:
        ww.write("\t")
    else:
        ww.write("\n")

# Slope write out
ww.write("\tslope = bias\t")
for tar in range(0, len(targets)):
    ww.write("{0:0.6f}".format(data[targets[tar]]["slope_bias"]))
    if tar < len(targets) - 1:
        ww.write("\t")
    else:
        ww.write("\n")
# correlation_R write out
ww.write("\tcorrelation = R\t")
for tar in range(0, len(targets)):
    ww.write("{0:0.6f}".format(data[targets[tar]]["correlation_R"]))
    if tar < len(targets) - 1:
        ww.write("\t")
    else:
        ww.write("\n")

ww.write(emptyLine)
ww.write("\tSres\t")
for tar in range(0, len(targets)):
    ww.write("{0:0.6f}".format(data[targets[tar]]["Sres"]))
    if tar < len(targets) - 1:
        ww.write("\t")
    else:
        ww.write("\n")
ww.write("\tSxx\t")
for tar in range(0, len(targets)):
    ww.write("30")
    if tar < len(targets) - 1:
        ww.write("\t")
    else:
        ww.write("\n")
ww.write("\tse of slope\t")
for tar in range(0, len(targets)):
    ww.write("{0:0.6f}".format(data[targets[tar]]["se_of_slope"]))
    if tar < len(targets) - 1:
        ww.write("\t")
    else:
        ww.write("\n")
ww.write("expected slope\t1\t" + shortEmptyLine)
ww.write("\tt\t")
for tar in range(0, len(targets)):
    ww.write("{0:0.6f}".format(data[targets[tar]]["t"]))
    if tar < len(targets) - 1:
        ww.write("\t")
    else:
        ww.write("\n")
ww.write("\tp (slope = 1)\t")
for tar in range(0, len(targets)):
    ww.write("{0:0.4e}".format(data[targets[tar]]["p"]))
    if tar < len(targets) - 1:
        ww.write("\t")
    else:
        ww.write("\n")

# Derivation from regression
ww.write(emptyLine)
ww.write("\tDerivation from regression\t" + shortEmptyLine)
ww.write("\tis independent of scaling\t" + shortEmptyLine)
ww.write("\tncount\t")
for tar in range(0, len(targets)):
    if tar < len(targets) - 1:
        ww.write("15\t")
    else:
        ww.write("15\n")
ww.write("\tngroups\t")
for tar in range(0, len(targets)):
    if tar < len(targets) - 1:
        ww.write("5\t")
    else:
        ww.write("5\n")
ww.write(emptyLine)
ww.write("\tvar\t")
for tar in range(0, len(targets)):
    ww.write("{0:0.6f}".format(data[targets[tar]]["all_variance"]))
    if tar < len(targets) - 1:
        ww.write("\t")
    else:
        ww.write("\n")
ww.write("\tSStotal\t")
for tar in range(0, len(targets)):
    ww.write("{0:0.6f}".format(data[targets[tar]]["SStotal"]))
    if tar < len(targets) - 1:
        ww.write("\t")
    else:
        ww.write("\n")
ww.write(emptyLine)
for concStr in ["15", "150", "1500", "15000", "150000"]:
    if concStr == "15":
        ww.write("\tvariance per concentration\t")
    else:
        ww.write("\t\t")
    for tar in range(0, len(targets)):
        ww.write("{0:0.6f}".format(data[targets[tar]]["variance"]["STD_" + concStr]))
        if tar < len(targets) - 1:
            ww.write("\t")
        else:
            ww.write("\n")
for concStr in ["15", "150", "1500", "15000", "150000"]:
    if concStr == "15":
        ww.write("\tSS per concentration\t")
    else:
        ww.write("\t\t")
    for tar in range(0, len(targets)):
        ww.write("{0:0.6f}".format(data[targets[tar]]["variance"]["STD_" + concStr] * 2))
        if tar < len(targets) - 1:
            ww.write("\t")
        else:
            ww.write("\n")
ww.write("\tSSwithin\t")
for tar in range(0, len(targets)):
    ww.write("{0:0.6f}".format(data[targets[tar]]["SSwithin"]))
    if tar < len(targets) - 1:
        ww.write("\t")
    else:
        ww.write("\n")
ww.write("\tSSbetween\t")
for tar in range(0, len(targets)):
    ww.write("{0:0.6f}".format(data[targets[tar]]["SSbetween"]))
    if tar < len(targets) - 1:
        ww.write("\t")
    else:
        ww.write("\n")
ww.write(emptyLine)
ww.write("\tstexy\t")
for tar in range(0, len(targets)):
    ww.write("{0:0.6f}".format(data[targets[tar]]["stexy"]))
    if tar < len(targets) - 1:
        ww.write("\t")
    else:
        ww.write("\n")
ww.write("\tSSres\t")
for tar in range(0, len(targets)):
    ww.write("{0:0.6f}".format(data[targets[tar]]["SSres"]))
    if tar < len(targets) - 1:
        ww.write("\t")
    else:
        ww.write("\n")
ww.write("\tSSregression\t")
for tar in range(0, len(targets)):
    ww.write("{0:0.6f}".format(data[targets[tar]]["SSregression"]))
    if tar < len(targets) - 1:
        ww.write("\t")
    else:
        ww.write("\n")
ww.write("\tSSderivation\t")
for tar in range(0, len(targets)):
    ww.write("{0:0.6f}".format(data[targets[tar]]["SSderivation"]))
    if tar < len(targets) - 1:
        ww.write("\t")
    else:
        ww.write("\n")

ww.write(emptyLine)
ww.write("16\t\t")
for tar in range(0, len(targets)):
    if tar < len(targets) - 1:
        ww.write("SS\t")
    else:
        ww.write("SS\n")
ww.write("\tbetween\t")
for tar in range(0, len(targets)):
    ww.write("{0:0.6f}".format(data[targets[tar]]["SSbetween"]))
    if tar < len(targets) - 1:
        ww.write("\t")
    else:
        ww.write("\n")
ww.write("\tregression\t")
for tar in range(0, len(targets)):
    ww.write("{0:0.6f}".format(data[targets[tar]]["SSregression"]))
    if tar < len(targets) - 1:
        ww.write("\t")
    else:
        ww.write("\n")
ww.write("\tderivation\t")
for tar in range(0, len(targets)):
    ww.write("{0:0.6f}".format(data[targets[tar]]["SSderivation"]))
    if tar < len(targets) - 1:
        ww.write("\t")
    else:
        ww.write("\n")
ww.write("\twithin\t")
for tar in range(0, len(targets)):
    ww.write("{0:0.6f}".format(data[targets[tar]]["SSwithin"]))
    if tar < len(targets) - 1:
        ww.write("\t")
    else:
        ww.write("\n")
ww.write(emptyLine)

ww.write("\t\t")
for tar in range(0, len(targets)):
    if tar < len(targets) - 1:
        ww.write("df\t")
    else:
        ww.write("df\n")
ww.write("\tbetween\t")
for tar in range(0, len(targets)):
    if tar < len(targets) - 1:
        ww.write("4\t")
    else:
        ww.write("4\n")
ww.write("\tregression\t")
for tar in range(0, len(targets)):
    if tar < len(targets) - 1:
        ww.write("1\t")
    else:
        ww.write("1\n")
ww.write("\tderivation\t")
for tar in range(0, len(targets)):
    if tar < len(targets) - 1:
        ww.write("3\t")
    else:
        ww.write("3\n")
ww.write("\twithin\t")
for tar in range(0, len(targets)):
    if tar < len(targets) - 1:
        ww.write("10\t")
    else:
        ww.write("10\n")
ww.write(emptyLine)
ww.write("\t\t")
for tar in range(0, len(targets)):
    if tar < len(targets) - 1:
        ww.write("MS\t")
    else:
        ww.write("MS\n")
ww.write("\tbetween\t")
for tar in range(0, len(targets)):
    ww.write("{0:0.6f}".format(data[targets[tar]]["MSbetween"]))
    if tar < len(targets) - 1:
        ww.write("\t")
    else:
        ww.write("\n")
ww.write("\tregression\t")
for tar in range(0, len(targets)):
    ww.write("{0:0.6f}".format(data[targets[tar]]["MSregression"]))
    if tar < len(targets) - 1:
        ww.write("\t")
    else:
        ww.write("\n")
ww.write("\tLinearity: MS derivation\t")
for tar in range(0, len(targets)):
    ww.write("{0:0.6f}".format(data[targets[tar]]["MSderivation"]))
    if tar < len(targets) - 1:
        ww.write("\t")
    else:
        ww.write("\n")
ww.write("\tReproducibility: MS within\t")
for tar in range(0, len(targets)):
    ww.write("{0:0.6f}".format(data[targets[tar]]["MSwithin"]))
    if tar < len(targets) - 1:
        ww.write("\t")
    else:
        ww.write("\n")
ww.close()

##################################
###  Write out variation tabe  ###
##################################
if not noCq:
    # Write target names
    wv.write("\t\t")
    for tar in range(0, len(targets)):
        wv.write(targets[tar])
        if tar < len(targets) - 1:
            wv.write("\t")
        else:
            wv.write("\n")

    # Write the input data
    for concStr in ["15", "150", "1500", "15000", "150000"]:
        for pos in range(0, 3):
            wv.write("\t" + concStr + "\t")
            for tar in range(0, len(targets)):
                wv.write("{0:0.2f}".format(data[targets[tar]]["Cq"]["STD_" + concStr][pos]))
                if tar < len(targets) - 1:
                    wv.write("\t")
                else:
                    wv.write("\n")

    # mean(max) write out
    wv.write(emptyLine)
    wv.write("\tstandard curve slope\t")
    for tar in range(0, len(targets)):
        wv.write("{0:0.3e}".format(data[targets[tar]]["cq_slope_bias"]))
        if tar < len(targets) - 1:
            wv.write("\t")
        else:
            wv.write("\n")
    wv.write("\tstandard curve derived eff\t")
    for tar in range(0, len(targets)):
        wv.write("{0:0.3e}".format(data[targets[tar]]["cq_eff"]))
        if tar < len(targets) - 1:
            wv.write("\t")
        else:
            wv.write("\n")

    wv.write(emptyLine)
    wv.write("\tFq\t")
    for tar in range(0, len(targets)):
        wv.write("1")
        if tar < len(targets) - 1:
            wv.write("\t")
        else:
            wv.write("\n")
    wv.write("\tlog(input) - log(expected F0 with SC-Eff\t" + shortEmptyLine)
    # Write the input data
    for concStr in ["15", "150", "1500", "15000", "150000"]:
        for pos in range(0, 3):
            wv.write("\t" + "{0:0.4f}".format(np.log10(float(concStr))) + "\t")
            for tar in range(0, len(targets)):
                wv.write("{0:0.4f}".format(data[targets[tar]]["Cq_F0"]["STD_" + concStr][pos]))
                if tar < len(targets) - 1:
                    wv.write("\t")
                else:
                    wv.write("\n")
    wv.write(emptyLine)
    wv.write("\tslope = bias\t")
    for tar in range(0, len(targets)):
        wv.write("{0:0.4f}".format(data[targets[tar]]["cq_slope_check"]))
        if tar < len(targets) - 1:
            wv.write("\t")
        else:
            wv.write("\n")
    wv.write(emptyLine)
    wv.write("\tn per gloup\t")
    for tar in range(0, len(targets)):
        wv.write("3")
        if tar < len(targets) - 1:
            wv.write("\t")
        else:
            wv.write("\n")
    wv.write("\t\t")
    for tar in range(0, len(targets)):
        wv.write("3")
        if tar < len(targets) - 1:
            wv.write("\t")
        else:
            wv.write("\n")
    wv.write("\t\t")
    for tar in range(0, len(targets)):
        wv.write("3")
        if tar < len(targets) - 1:
            wv.write("\t")
        else:
            wv.write("\n")
    wv.write("\t\t")
    for tar in range(0, len(targets)):
        wv.write("3")
        if tar < len(targets) - 1:
            wv.write("\t")
        else:
            wv.write("\n")
    wv.write("\t\t")
    for tar in range(0, len(targets)):
        wv.write("3")
        if tar < len(targets) - 1:
            wv.write("\t")
        else:
            wv.write("\n")
    wv.write("\tngroup\t")
    for tar in range(0, len(targets)):
        wv.write("5")
        if tar < len(targets) - 1:
            wv.write("\t")
        else:
            wv.write("\n")
    wv.write(emptyLine)
    for concStr in ["15", "150", "1500", "15000", "150000"]:
        if concStr == "15":
            wv.write("\tSS per dilution\t")
        else:
            wv.write("\t\t")
        for tar in range(0, len(targets)):
            wv.write("{0:0.6f}".format(data[targets[tar]]["cq_ss_per_dil"]["STD_" + concStr]))
            if tar < len(targets) - 1:
                wv.write("\t")
            else:
                wv.write("\n")
    wv.write("\tSSwithin\t")
    for tar in range(0, len(targets)):
        wv.write("{0:0.6f}".format(data[targets[tar]]["cq_SSwithin"]))
        if tar < len(targets) - 1:
            wv.write("\t")
        else:
            wv.write("\n")
    wv.write("\tdf-within\t")
    for tar in range(0, len(targets)):
        wv.write("10")
        if tar < len(targets) - 1:
            wv.write("\t")
        else:
            wv.write("\n")
    wv.write("\tMS within Cq-standard\t")
    for tar in range(0, len(targets)):
        wv.write("{0:0.6f}".format(data[targets[tar]]["cq_ms_within_cq"]))
        if tar < len(targets) - 1:
            wv.write("\t")
        else:
            wv.write("\n")
    wv.write(emptyLine)
    wv.write(emptyLine)
    wv.write(emptyLine)


    # Write target names
    wv.write("\t\t")
    for tar in range(0, len(targets)):
        wv.write(targets[tar])
        if tar < len(targets) - 1:
            wv.write("\t")
        else:
            wv.write("\n")

    # Write the input data
    for concStr in ["15", "150", "1500", "15000", "150000"]:
        for pos in range(0, 3):
            wv.write("observed F0\t" + concStr + "\t")
            for tar in range(0, len(targets)):
                wv.write("{0:0.3e}".format(data[targets[tar]]["raw"]["STD_" + concStr][pos]))
                if tar < len(targets) - 1:
                    wv.write("\t")
                else:
                    wv.write("\n")
    wv.write(emptyLine)
    for concStr in ["15", "150", "1500", "15000", "150000"]:
        for pos in range(0, 3):
            wv.write("log(observed F0)\t" + concStr + "\t")
            for tar in range(0, len(targets)):
                wv.write("{0:0.6f}".format(np.log10(data[targets[tar]]["raw"]["STD_" + concStr][pos])))
                if tar < len(targets) - 1:
                    wv.write("\t")
                else:
                    wv.write("\n")
    wv.write(emptyLine)

    wv.write(emptyLine)
    for concStr in ["15", "150", "1500", "15000", "150000"]:
        if concStr == "15":
            wv.write("\tSS per dilution\t")
        else:
            wv.write("\t\t")
        for tar in range(0, len(targets)):
            wv.write("{0:0.6f}".format(data[targets[tar]]["variance"]["STD_" + concStr] * 2))
            if tar < len(targets) - 1:
                wv.write("\t")
            else:
                wv.write("\n")
    wv.write("\tSSwithin\t")
    for tar in range(0, len(targets)):
        wv.write("{0:0.6f}".format(data[targets[tar]]["SSwithin"]))
        if tar < len(targets) - 1:
            wv.write("\t")
        else:
            wv.write("\n")
    wv.write("\tdf-within\t")
    for tar in range(0, len(targets)):
        wv.write("10")
        if tar < len(targets) - 1:
            wv.write("\t")
        else:
            wv.write("\n")
    wv.write("\tMS within LinRegPCR\t")
    for tar in range(0, len(targets)):
        wv.write("{0:0.6f}".format(data[targets[tar]]["cq_ms_within_linregpcr"]))
        if tar < len(targets) - 1:
            wv.write("\t")
        else:
            wv.write("\n")

    wv.write(emptyLine)
    # Write target names
    wv.write("\t\t")
    for tar in range(0, len(targets)):
        wv.write(targets[tar])
        if tar < len(targets) - 1:
            wv.write("\t")
        else:
            wv.write("\n")
    wv.write("\tratio = measure for added or reduced variability\t")
    for tar in range(0, len(targets)):
        wv.write("{0:0.6f}".format(data[targets[tar]]["cq_ms_ratio"]))
        if tar < len(targets) - 1:
            wv.write("\t")
        else:
            wv.write("\n")

    wv.close()

####################################
###  Write out differences tabe  ###
####################################

# Write target names
wd.write("\t\t\t")
for tar in range(0, len(targets)):
    wd.write(targets[tar])
    if tar < len(targets) - 1:
        wd.write("\t")
    else:
        wd.write("\n")

# Write the input data
for concStr in ["15", "150", "1500", "15000", "150000"]:
    for pos in range(0, 3):
        wd.write("\t\t" + concStr + "\t")
        for tar in range(0, len(targets)):
            wd.write("{0:0.3e}".format(data[targets[tar]]["raw"]["STD_" + concStr][pos]))
            if tar < len(targets) - 1:
                wd.write("\t")
            else:
                wd.write("\n")
wd.write("\t" + emptyLine)
wd.write("\t" + emptyLine)
for concStr in ["15", "150", "1500", "15000", "150000"]:
    for pos in range(0, 3):
        wd.write("\t\t" + str(np.log10(float(concStr))) + "\t")
        for tar in range(0, len(targets)):
            wd.write("{0:0.6f}".format(np.log10(data[targets[tar]]["raw"]["STD_" + concStr][pos])))
            if tar < len(targets) - 1:
                wd.write("\t")
            else:
                wd.write("\n")

wd.write("\t" + emptyLine)
wd.write("\t\tsteyx\t")
for tar in range(0, len(targets)):
    wd.write("{0:0.6f}".format(data[targets[tar]]["stexy"]))
    if tar < len(targets) - 1:
        wd.write("\t")
    else:
        wd.write("\n")
wd.write("\t" + emptyLine)
for concStr in ["15", "150", "1500", "15000", "150000"]:
    wd.write("\t\t" + str(np.log10(float(concStr))) + "\t")
    for tar in range(0, len(targets)):
        wd.write("{0:0.6f}".format(data[targets[tar]]["mean_log"]["STD_" + concStr]))
        if tar < len(targets) - 1:
            wd.write("\t")
        else:
            wd.write("\n")

dif_alpha = 0.95
dif_detect_num = 0
dif_detect_sum = 0.0
dif_SSx = np.var(np.log10(data[targets[tar]]["all_std"]), ddof=1) * (data[targets[tar]]["var_n"] - 2.0)
dif_x_mean = np.mean(np.log10(data[targets[tar]]["std"]["val"]))
dif_n = 3
dif_t = -1.0 * scp.t.ppf((1.0 - dif_alpha) / 4.0, dif_n - 1)

wd.write("x-mean\t" + "{0:0.6f}".format(dif_x_mean) + "\t\t" + shortEmptyLine)
wd.write("SSx\t" + "{0:0.6f}".format(dif_SSx) + "\t\t" + shortEmptyLine)
wd.write("n\t" + "{0:0.6f}".format(dif_n) + "\t\t" + shortEmptyLine)
wd.write("alpha\t" + "{0:0.6f}".format(dif_alpha) + "\t\t" + shortEmptyLine)
wd.write("t\t" + "{0:0.6f}".format(dif_t) + "\t\t" + shortEmptyLine)

for concStr in ["15", "150", "1500", "15000", "150000"]:
    wd.write("\t\t\t")
    for tar in range(0, len(targets)):
        wd.write("{0:0.6f}".format(data[targets[tar]]["dif_se_yfit"]["STD_" + concStr]))
        if tar < len(targets) - 1:
            wd.write("\t")
        else:
            wd.write("\n")
wd.write("\t" + emptyLine)
wd.write("\t\tlog(upper)\t" + shortEmptyLine)
for concStr in ["15", "150", "1500", "15000", "150000"]:
    wd.write("\t\t\t")
    for tar in range(0, len(targets)):
        wd.write("{0:0.6f}".format(data[targets[tar]]["dif_log_upper"]["STD_" + concStr]))
        if tar < len(targets) - 1:
            wd.write("\t")
        else:
            wd.write("\n")
wd.write("\t\tlog(lower)\t" + shortEmptyLine)
for concStr in ["15", "150", "1500", "15000", "150000"]:
    wd.write("\t\t\t")
    for tar in range(0, len(targets)):
        wd.write("{0:0.6f}".format(data[targets[tar]]["dif_log_lower"]["STD_" + concStr]))
        if tar < len(targets) - 1:
            wd.write("\t")
        else:
            wd.write("\n")

wd.write("\t" + emptyLine)
wd.write("\t\tupper\t" + shortEmptyLine)
for concStr in ["15", "150", "1500", "15000", "150000"]:
    wd.write("\t\t\t")
    for tar in range(0, len(targets)):
        wd.write("{0:0.3e}".format(np.power(10.0, data[targets[tar]]["dif_log_upper"]["STD_" + concStr])))
        if tar < len(targets) - 1:
            wd.write("\t")
        else:
            wd.write("\n")
wd.write("\t\tlower\t" + shortEmptyLine)
for concStr in ["15", "150", "1500", "15000", "150000"]:
    wd.write("\t\t\t")
    for tar in range(0, len(targets)):
        wd.write("{0:0.3e}".format(np.power(10.0, data[targets[tar]]["dif_log_lower"]["STD_" + concStr])))
        if tar < len(targets) - 1:
            wd.write("\t")
        else:
            wd.write("\n")
wd.write("\t\tmean\t" + shortEmptyLine)
for concStr in ["15", "150", "1500", "15000", "150000"]:
    wd.write("\t\t\t")
    for tar in range(0, len(targets)):
        wd.write("{0:0.3e}".format(np.power(10.0, data[targets[tar]]["mean_log"]["STD_" + concStr])))
        if tar < len(targets) - 1:
            wd.write("\t")
        else:
            wd.write("\n")
wd.write("\t\tfold up\t" + shortEmptyLine)
for concStr in ["15", "150", "1500", "15000", "150000"]:
    wd.write("\t\t\t")
    for tar in range(0, len(targets)):
        wd.write("{0:0.6f}".format(data[targets[tar]]["dif_fold_up"]["STD_" + concStr]))
        if tar < len(targets) - 1:
            wd.write("\t")
        else:
            wd.write("\n")
wd.write("\t\tfold down\t" + shortEmptyLine)
for concStr in ["15", "150", "1500", "15000", "150000"]:
    wd.write("\t\t\t")
    for tar in range(0, len(targets)):
        wd.write("{0:0.6f}".format(data[targets[tar]]["dif_fold_down"]["STD_" + concStr]))
        if tar < len(targets) - 1:
            wd.write("\t")
        else:
            wd.write("\n")

wd.write("\t" + emptyLine)
# Write target names
wd.write("\t\t\t")
for tar in range(0, len(targets)):
    wd.write(targets[tar])
    if tar < len(targets) - 1:
        wd.write("\t")
    else:
        wd.write("\n")
wd.write("\t\tdetectable difference\t")
for tar in range(0, len(targets)):
    wd.write("{0:0.6f}".format(data[targets[tar]]["dif_detectable_diff"]))
    if tar < len(targets) - 1:
        wd.write("\t")
    else:
        wd.write("\n")

wd.close()
sys.exit(0)
