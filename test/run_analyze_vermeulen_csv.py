#!/usr/bin/python

from __future__ import division
from __future__ import print_function

import os
import json
import sys
import time
import math
import csv
import numpy as np
import scipy
import scipy.stats

use_results = True  # True - results from paper, False - use output form run_analyze_vermeulen_rdml.py
set_method = 1  # 1 - LinRegPCR, 2 - 5PSM, 3 - FPK-PCR, 4 - LRE-Emax, 5 - LRE-E100, 6 - Cy0, 7 - MAK2, 8 - DART, 9 - 4PLM, 10 - PCR-Miner

out_file_bias = "temp_vermeulen_report_bias.csv"
out_file_var = "temp_vermeulen_report_variation.csv"
out_file_diff = "temp_vermeulen_report_difference.csv"

ww = open(out_file_bias, "w")
wv = open(out_file_var, "w")
wd = open(out_file_diff, "w")

data = {}

if use_results:
    with open("data_vermeulen_methods_quant.csv", newline='') as tfile:
        table = list(csv.reader(tfile, delimiter=";"))
        startRow = 17 * (set_method - 1) + 1
        for row in range(startRow + 1, startRow + 16):
            for col in range(3, len(table[row])):
                if table[startRow][col] not in data:
                    data[table[startRow][col]] = {}
                if table[row][2] not in data[table[startRow][col]]:
                    data[table[startRow][col]][table[row][2]] = {}
                    data[table[startRow][col]][table[row][2]]["raw"] = []
                try:
                    vFloat = float(table[row][col])
                except ValueError:
                    print("Error: " + table[startRow][col] + " - " + table[row][2] + " = " + table[row][col])
                else:
                    if math.isfinite(vFloat):
                        data[table[startRow][col]][table[row][2]]["raw"].append(vFloat)
                    else:
                        print("Error: " + table[startRow][col] + " - " + table[row][2] + " = " + table[row][col])
    with open("data_vermeulen_methods_cq.csv", newline='') as tfile:
        table = list(csv.reader(tfile, delimiter=";"))
        startRow = 17 * (set_method - 1) + 1
        for row in range(startRow + 1, startRow + 16):
            for col in range(3, len(table[row])):
                if table[startRow][col] not in data:
                    data[table[startRow][col]] = {}
                if table[row][2] not in data[table[startRow][col]]:
                    data[table[startRow][col]][table[row][2]] = {}
                if "Cq" not in data[table[startRow][col]][table[row][2]]:
                    data[table[startRow][col]][table[row][2]]["Cq"] = []
                try:
                    vFloat = float(table[row][col])
                except ValueError:
                    print("Error Cq: " + table[startRow][col] + " - " + table[row][2] + " = " + table[row][col])
                else:
                    if math.isfinite(vFloat):
                        data[table[startRow][col]][table[row][2]]["Cq"].append(vFloat)
                    else:
                        print("Error Cq: " + table[startRow][col] + " - " + table[row][2] + " = " + table[row][col])

targets = sorted(list(data.keys()))

# Calculate everything
for tar in range(0, len(targets)):
    data[targets[tar]]["mean_max"] = np.mean(data[targets[tar]]["STD_150000"]["raw"])
    log_x = []
    log_y = []
    concSSsum = 0.0
    for conc in ["15", "150", "1500", "15000", "150000"]:
        data[targets[tar]]["STD_" + conc]["scaled"] = []
        for pos in range(0, 3):
            data[targets[tar]]["STD_" + conc]["scaled"].append(data[targets[tar]]["STD_" + conc]["raw"][pos] / data[targets[tar]]["mean_max"])

        data[targets[tar]]["STD_" + conc]["variance"] = np.var(np.log10(data[targets[tar]]["STD_" + conc]["scaled"]), ddof=1)
        data[targets[tar]]["STD_" + conc]["dif_log_mean"] = np.mean(np.log10(data[targets[tar]]["STD_" + conc]["raw"]))
        concSSsum += data[targets[tar]]["STD_" + conc]["variance"]
        for pos in range(0, 3):
            log_x.append(np.log10(float(conc) / 150000))
            log_y.append(np.log10(data[targets[tar]]["STD_" + conc]["scaled"][pos]))
    slope, intercept, r_val, p_val, std_err = scipy.stats.linregress(log_x, log_y)
    data[targets[tar]]["slope_bias"] = slope
    data[targets[tar]]["correlation_R"] = r_val
    data[targets[tar]]["scaled_mean_min"] = np.mean(data[targets[tar]]["STD_15"]["scaled"])
    data[targets[tar]]["scaled_mean_max"] = np.mean(data[targets[tar]]["STD_150000"]["scaled"])
    data[targets[tar]]["bias_as_ratio"] = data[targets[tar]]["scaled_mean_max"] / data[targets[tar]]["scaled_mean_min"]
    data[targets[tar]]["variance"] = np.var(log_y, ddof=1)
    data[targets[tar]]["SStotal"] = data[targets[tar]]["variance"] * (15 - 1)
    data[targets[tar]]["SSwithin"] = concSSsum * 2
    data[targets[tar]]["SSbetween"] = data[targets[tar]]["SStotal"] - data[targets[tar]]["SSwithin"]
    # Calc stexy
    fit = np.polyfit(log_x, log_y, deg=1)
    nnn = len(log_x)
    mmm = fit[0]
    ccc = fit[1]
    y_pred = mmm * np.asarray(log_x) + ccc
    data[targets[tar]]["stexy"] = (((np.asarray(log_y) - y_pred)**2).sum() / (nnn - 2))**0.5
    data[targets[tar]]["SSres"] = np.power(data[targets[tar]]["stexy"], 2) * (15-2)
    data[targets[tar]]["Sres"] = np.power(data[targets[tar]]["SSres"] / (15 - 2), 0.5)
    data[targets[tar]]["se_of_slope"] = data[targets[tar]]["Sres"] / np.power(30, 0.5)
    data[targets[tar]]["t"] = np.abs(1 - data[targets[tar]]["slope_bias"]) / data[targets[tar]]["se_of_slope"]
    data[targets[tar]]["p"] = scipy.stats.t.sf(data[targets[tar]]["t"], df=(15-2)) * 2
    data[targets[tar]]["SSregression"] = data[targets[tar]]["SStotal"] - data[targets[tar]]["SSres"]
    data[targets[tar]]["SSderivation"] = data[targets[tar]]["SSres"] - data[targets[tar]]["SSwithin"]
    data[targets[tar]]["MSbetween"] = data[targets[tar]]["SSbetween"] / 4
    data[targets[tar]]["MSregression"] = data[targets[tar]]["SSregression"] / 1
    data[targets[tar]]["MSderivation"] = data[targets[tar]]["SSderivation"] / 3
    data[targets[tar]]["MSwithin"] = data[targets[tar]]["SSwithin"] / 10

    # Cq variance calc
    cq_x = []
    cq_y = []
    cq_check_y = []
    for conc in ["15", "150", "1500", "15000", "150000"]:
        for pos in range(0, 3):
            cq_x.append(np.log10(float(conc) / 150000))
            cq_y.append(data[targets[tar]]["STD_" + conc]["Cq"][pos])
    slope, intercept, r_val, p_val, std_err = scipy.stats.linregress(cq_x, cq_y)
    data[targets[tar]]["cq_slope_bias"] = slope
    data[targets[tar]]["cq_eff"] = np.power(10, -1.0 / data[targets[tar]]["cq_slope_bias"])
    data[targets[tar]]["cq_SSwithin"] = 0.0
    for conc in ["15", "150", "1500", "15000", "150000"]:
        if "Cq_F0" not in data[targets[tar]]["STD_" + conc]:
            data[targets[tar]]["STD_" + conc]["Cq_F0"] = []
        for pos in range(0, 3):
            data[targets[tar]]["STD_" + conc]["Cq_F0"].append(np.log10(1.0 / np.power(data[targets[tar]]["cq_eff"], data[targets[tar]]["STD_" + conc]["Cq"][pos])))
            cq_check_y.append(np.log10(1.0 / np.power(data[targets[tar]]["cq_eff"], data[targets[tar]]["STD_" + conc]["Cq"][pos])))
        data[targets[tar]]["STD_" + conc]["cq_ss_per_dil"] = np.var(data[targets[tar]]["STD_" + conc]["Cq_F0"], ddof=1) * (3 - 1)
        data[targets[tar]]["cq_SSwithin"] += data[targets[tar]]["STD_" + conc]["cq_ss_per_dil"]
        
    slope, intercept, r_val, p_val, std_err = scipy.stats.linregress(cq_x, cq_check_y)
    data[targets[tar]]["cq_slope_check"] = slope
    data[targets[tar]]["cq_ms_within_cq"] = data[targets[tar]]["cq_SSwithin"] / 10.0
    data[targets[tar]]["cq_ms_within_linregpcr"] = data[targets[tar]]["SSwithin"] / 10.0
    data[targets[tar]]["cq_ms_ratio"] = data[targets[tar]]["cq_ms_within_linregpcr"] / data[targets[tar]]["cq_ms_within_cq"]

    # difference calc
    dif_conc = [15.0, 150.0, 1500.0, 15000.0, 150000.0]
    dif_conc2 = [15.0, 15.0, 15.0, 150.0, 150.0, 150.0, 1500.0, 1500.0, 1500.0, 15000.0, 15000.0, 15000.0, 150000.0, 150000.0, 150000.0]
    dif_x_mean = np.mean(np.log10(np.asarray(dif_conc)))
    dif_SSx = np.var(np.log10(np.asarray(dif_conc2)), ddof=1) * (15.0 - 2.0)
    dif_n = 3
    dif_alpha = 0.95
    dif_t = -1.0 * scipy.stats.t.ppf((1.0 - dif_alpha) / 4.0, dif_n - 1)
    dif_detect_num = 0
    dif_detect_sum = 0.0
    for conc in ["15", "150", "1500", "15000", "150000"]:
        data[targets[tar]]["STD_" + conc]["scaled"] = []
        for pos in range(0, 3):
            data[targets[tar]]["STD_" + conc]["scaled"].append(data[targets[tar]]["STD_" + conc]["raw"][pos] / data[targets[tar]]["mean_max"])

        data[targets[tar]]["STD_" + conc]["variance"] = np.var(np.log10(data[targets[tar]]["STD_" + conc]["scaled"]), ddof=1)
        data[targets[tar]]["STD_" + conc]["dif_se_yfit"] = data[targets[tar]]["stexy"] * np.power(1.0 / dif_n + np.power(np.log10(float(conc)) - dif_x_mean, 2.0) / dif_SSx , 0.5)
        data[targets[tar]]["STD_" + conc]["dif_log_upper"] = data[targets[tar]]["STD_" + conc]["dif_log_mean"] + dif_t * data[targets[tar]]["STD_" + conc]["dif_se_yfit"]
        data[targets[tar]]["STD_" + conc]["dif_log_lower"] = data[targets[tar]]["STD_" + conc]["dif_log_mean"] - dif_t * data[targets[tar]]["STD_" + conc]["dif_se_yfit"]
        data[targets[tar]]["STD_" + conc]["dif_fold_up"] = np.power(10.0, data[targets[tar]]["STD_" + conc]["dif_log_upper"]) / np.power(10.0, data[targets[tar]]["STD_" + conc]["dif_log_mean"])
        dif_detect_sum += np.log(data[targets[tar]]["STD_" + conc]["dif_fold_up"])
        dif_detect_num += 1
        data[targets[tar]]["STD_" + conc]["dif_fold_down"] = np.power(10.0, data[targets[tar]]["STD_" + conc]["dif_log_mean"]) / np.power(10.0, data[targets[tar]]["STD_" + conc]["dif_log_lower"])

    data[targets[tar]]["dif_detectable_diff"] = np.exp(dif_detect_sum / dif_detect_num)

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
for conc in ["15", "150", "1500", "15000", "150000"]:
    for pos in range(0, 3):
        ww.write("\t" + conc + "\t")
        for tar in range(0, len(targets)):
            ww.write("{0:0.3e}".format(data[targets[tar]]["STD_" + conc]["raw"][pos]))
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
for conc in ["15", "150", "1500", "15000", "150000"]:
    for pos in range(0, 3):
        ww.write("\t" + str(float(conc) / 150000) + "\t")
        for tar in range(0, len(targets)):
            ww.write("{0:0.6f}".format(data[targets[tar]]["STD_" + conc]["scaled"][pos]))
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
for conc in ["15", "150", "1500", "15000", "150000"]:
    for pos in range(0, 3):
        ww.write("\t" + str(np.log10(float(conc) / 150000)) + "\t")
        for tar in range(0, len(targets)):
            ww.write("{0:0.6f}".format(np.log10(data[targets[tar]]["STD_" + conc]["scaled"][pos])))
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
    ww.write("{0:0.6f}".format(data[targets[tar]]["variance"]))
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
for conc in ["15", "150", "1500", "15000", "150000"]:
    if conc == "15":
        ww.write("\tvariance per concentration\t")
    else:
        ww.write("\t\t")
    for tar in range(0, len(targets)):
        ww.write("{0:0.6f}".format(data[targets[tar]]["STD_" + conc]["variance"]))
        if tar < len(targets) - 1:
            ww.write("\t")
        else:
            ww.write("\n")
for conc in ["15", "150", "1500", "15000", "150000"]:
    if conc == "15":
        ww.write("\tSS per concentration\t")
    else:
        ww.write("\t\t")
    for tar in range(0, len(targets)):
        ww.write("{0:0.6f}".format(data[targets[tar]]["STD_" + conc]["variance"] * 2))
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

# Write target names
wv.write("\t\t")
for tar in range(0, len(targets)):
    wv.write(targets[tar])
    if tar < len(targets) - 1:
        wv.write("\t")
    else:
        wv.write("\n")

# Write the input data
for conc in ["15", "150", "1500", "15000", "150000"]:
    for pos in range(0, 3):
        wv.write("\t" + conc + "\t")
        for tar in range(0, len(targets)):
            wv.write("{0:0.2f}".format(data[targets[tar]]["STD_" + conc]["Cq"][pos]))
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
for conc in ["15", "150", "1500", "15000", "150000"]:
    for pos in range(0, 3):
        wv.write("\t" + "{0:0.4f}".format(np.log10(float(conc))) + "\t")
        for tar in range(0, len(targets)):
            wv.write("{0:0.4f}".format(data[targets[tar]]["STD_" + conc]["Cq_F0"][pos]))
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
for conc in ["15", "150", "1500", "15000", "150000"]:
    if conc == "15":
        wv.write("\tSS per dilution\t")
    else:
        wv.write("\t\t")
    for tar in range(0, len(targets)):
        wv.write("{0:0.6f}".format(data[targets[tar]]["STD_" + conc]["cq_ss_per_dil"]))
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
for conc in ["15", "150", "1500", "15000", "150000"]:
    for pos in range(0, 3):
        wv.write("observed F0\t" + conc + "\t")
        for tar in range(0, len(targets)):
            wv.write("{0:0.3e}".format(data[targets[tar]]["STD_" + conc]["raw"][pos]))
            if tar < len(targets) - 1:
                wv.write("\t")
            else:
                wv.write("\n")
wv.write(emptyLine)
for conc in ["15", "150", "1500", "15000", "150000"]:
    for pos in range(0, 3):
        wv.write("log(observed F0)\t" + conc + "\t")
        for tar in range(0, len(targets)):
            wv.write("{0:0.6f}".format(np.log10(data[targets[tar]]["STD_" + conc]["raw"][pos])))
            if tar < len(targets) - 1:
                wv.write("\t")
            else:
                wv.write("\n")
wv.write(emptyLine)

wv.write(emptyLine)
for conc in ["15", "150", "1500", "15000", "150000"]:
    if conc == "15":
        wv.write("\tSS per dilution\t")
    else:
        wv.write("\t\t")
    for tar in range(0, len(targets)):
        wv.write("{0:0.6f}".format(data[targets[tar]]["STD_" + conc]["variance"] * 2))
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
for conc in ["15", "150", "1500", "15000", "150000"]:
    for pos in range(0, 3):
        wd.write("\t\t" + conc + "\t")
        for tar in range(0, len(targets)):
            wd.write("{0:0.3e}".format(data[targets[tar]]["STD_" + conc]["raw"][pos]))
            if tar < len(targets) - 1:
                wd.write("\t")
            else:
                wd.write("\n")
wd.write("\t" + emptyLine)
wd.write("\t" + emptyLine)
for conc in ["15", "150", "1500", "15000", "150000"]:
    for pos in range(0, 3):
        wd.write("\t\t" + str(np.log10(float(conc))) + "\t")
        for tar in range(0, len(targets)):
            wd.write("{0:0.6f}".format(np.log10(data[targets[tar]]["STD_" + conc]["raw"][pos])))
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
for conc in ["15", "150", "1500", "15000", "150000"]:
    wd.write("\t\t" + str(np.log10(float(conc))) + "\t")
    for tar in range(0, len(targets)):
        wd.write("{0:0.6f}".format(data[targets[tar]]["STD_" + conc]["dif_log_mean"]))
        if tar < len(targets) - 1:
            wd.write("\t")
        else:
            wd.write("\n")
wd.write("x-mean\t" + "{0:0.6f}".format(dif_x_mean) + "\t\t" + shortEmptyLine)
wd.write("SSx\t" + "{0:0.6f}".format(dif_SSx) + "\t\t" + shortEmptyLine)
wd.write("n\t" + "{0:0.6f}".format(dif_n) + "\t\t" + shortEmptyLine)
wd.write("alpha\t" + "{0:0.6f}".format(dif_alpha) + "\t\t" + shortEmptyLine)
wd.write("t\t" + "{0:0.6f}".format(dif_t) + "\t\t" + shortEmptyLine)

for conc in ["15", "150", "1500", "15000", "150000"]:
    wd.write("\t\t\t")
    for tar in range(0, len(targets)):
        wd.write("{0:0.6f}".format(data[targets[tar]]["STD_" + conc]["dif_se_yfit"]))
        if tar < len(targets) - 1:
            wd.write("\t")
        else:
            wd.write("\n")
wd.write("\t" + emptyLine)
wd.write("\t\tlog(upper)\t" + shortEmptyLine)
for conc in ["15", "150", "1500", "15000", "150000"]:
    wd.write("\t\t\t")
    for tar in range(0, len(targets)):
        wd.write("{0:0.6f}".format(data[targets[tar]]["STD_" + conc]["dif_log_upper"]))
        if tar < len(targets) - 1:
            wd.write("\t")
        else:
            wd.write("\n")
wd.write("\t\tlog(lower)\t" + shortEmptyLine)
for conc in ["15", "150", "1500", "15000", "150000"]:
    wd.write("\t\t\t")
    for tar in range(0, len(targets)):
        wd.write("{0:0.6f}".format(data[targets[tar]]["STD_" + conc]["dif_log_lower"]))
        if tar < len(targets) - 1:
            wd.write("\t")
        else:
            wd.write("\n")

wd.write("\t" + emptyLine)
wd.write("\t\tupper\t" + shortEmptyLine)
for conc in ["15", "150", "1500", "15000", "150000"]:
    wd.write("\t\t\t")
    for tar in range(0, len(targets)):
        wd.write("{0:0.3e}".format(np.power(10.0, data[targets[tar]]["STD_" + conc]["dif_log_upper"])))
        if tar < len(targets) - 1:
            wd.write("\t")
        else:
            wd.write("\n")
wd.write("\t\tlower\t" + shortEmptyLine)
for conc in ["15", "150", "1500", "15000", "150000"]:
    wd.write("\t\t\t")
    for tar in range(0, len(targets)):
        wd.write("{0:0.3e}".format(np.power(10.0, data[targets[tar]]["STD_" + conc]["dif_log_lower"])))
        if tar < len(targets) - 1:
            wd.write("\t")
        else:
            wd.write("\n")
wd.write("\t\tmean\t" + shortEmptyLine)
for conc in ["15", "150", "1500", "15000", "150000"]:
    wd.write("\t\t\t")
    for tar in range(0, len(targets)):
        wd.write("{0:0.3e}".format(np.power(10.0, data[targets[tar]]["STD_" + conc]["dif_log_mean"])))
        if tar < len(targets) - 1:
            wd.write("\t")
        else:
            wd.write("\n")
wd.write("\t\tfold up\t" + shortEmptyLine)
for conc in ["15", "150", "1500", "15000", "150000"]:
    wd.write("\t\t\t")
    for tar in range(0, len(targets)):
        wd.write("{0:0.6f}".format(data[targets[tar]]["STD_" + conc]["dif_fold_up"]))
        if tar < len(targets) - 1:
            wd.write("\t")
        else:
            wd.write("\n")
wd.write("\t\tfold down\t" + shortEmptyLine)
for conc in ["15", "150", "1500", "15000", "150000"]:
    wd.write("\t\t\t")
    for tar in range(0, len(targets)):
        wd.write("{0:0.6f}".format(data[targets[tar]]["STD_" + conc]["dif_fold_down"]))
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
