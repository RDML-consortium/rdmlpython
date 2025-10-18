#!/usr/bin/python

from __future__ import division
from __future__ import print_function

import os
import json
import sys
import time
import math
import numpy as np

parent_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
library_dir = os.path.abspath(os.path.join(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir), os.pardir))
sys.path.append(library_dir)

import rdmlpython as rdml

printExpRun = True

vermeulen_file = os.path.join(parent_dir, "experiments/vermeulen/data_vermeulen_raw.rdml")
vermeulen_std_5 = os.path.join(parent_dir, "experiments/vermeulen/data_vermeulen_std_150000_15.rdml")
vermeulen_std_4 = os.path.join(parent_dir, "experiments/vermeulen/data_vermeulen_std_150000_150.rdml")
vermeulen_std_3 = os.path.join(parent_dir, "experiments/vermeulen/data_vermeulen_std_150000_1500.rdml")
vermeulen_std_2 = os.path.join(parent_dir, "experiments/vermeulen/data_vermeulen_std_150000_15000.rdml")
vermeulen_std_1 = os.path.join(parent_dir, "experiments/vermeulen/data_vermeulen_std_150000_150000.rdml")
rdml_vol_file = os.path.join(parent_dir, "experiments/untergasser/volume_machine.rdml")
rdml_amp_file = os.path.join(parent_dir, "experiments/untergasser/amplicon_primer_mix.rdml")
rdml_dil_file = os.path.join(parent_dir, "experiments/untergasser/large_DNA_dilutions.rdml")
rdml_probes_file = os.path.join(parent_dir, "experiments/untergasser/probes.rdml")
rdml_add_sybr_file = os.path.join(parent_dir, "experiments/untergasser/sybr_conc.rdml")
rdml_ownmix_a_file = os.path.join(parent_dir, "experiments/untergasser/own_mix_a.rdml")
rdml_ownmix_b_file = os.path.join(parent_dir, "experiments/untergasser/own_mix_b.rdml")
rdml_ownmix_c_file = os.path.join(parent_dir, "experiments/untergasser/own_mix_c.rdml")
rdml_eva_a_file = os.path.join(parent_dir, "experiments/untergasser/eva_green_a.rdml")
rdml_eva_b_file = os.path.join(parent_dir, "experiments/untergasser/eva_green_b.rdml")
rdml_eva_c_file = os.path.join(parent_dir, "experiments/untergasser/eva_green_c.rdml")
out_vol_file = "temp_volume_resuls.csv"
out_amp_file = "temp_amplicon_primer_resuls.csv"
out_dil_file = "temp_large_DNA_dilutions_resuls.csv"
out_probes_file = "temp_probes_resuls.csv"
out_add_sybr_file = "temp_sybr_conc_resuls.csv"
out_ownmix_a_file = "temp_ownmix_a_resuls.csv"
out_ownmix_b_file = "temp_ownmix_b_resuls.csv"
out_ownmix_c_file = "temp_ownmix_c_resuls.csv"
out_eva_a_file = "temp_eva_green_a_resuls.csv"
out_eva_b_file = "temp_eva_green_b_resuls.csv"
out_eva_c_file = "temp_eva_green_c_resuls.csv"
out_file = "temp_vermeulen_resuls.csv"
in_json = "stored_results_test.json"
out_json = "temp_test.json"

vermeulenNcopyCol = 24

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




def saveCSV(arr, fil):
    outF = open(fil, "w")
    for row in range(0, len(arr)):
        for col in range(0, len(arr[row])):
            outF.write(str(arr[row][col]) + "\t")
        outF.write("\n")
    outF.write("\n")
    outF.close()


def saveNum(store, ele):
    if ele not in store:
        return -1.0
    else:
        return store[ele]


def stringNice(text, orgVal, storeVal, dist, dir, prec):
    aa = ""
    bb = ""
    diff = 0
    if dir == "+":
        diff = orgVal - storeVal
    else:
        diff = storeVal - orgVal
    if diff - dist > 0.0:
        aa = '\033[42m'
        bb = '\033[0m'
    elif diff + dist < 0.0:
        aa = '\033[41m'
        bb = '\033[0m'
    if prec == 0:
        return(aa + text + str(orgVal) + " (" + str(storeVal) + ")" + bb)
    return ""


def printNice(text, orgVal, storeVal, dist, dir, prec):
    print(stringNice(text, orgVal, storeVal, dist, dir, prec))


def stringNice2(orgStore, storeStore, ele, form, dist, dir, org):
    orgVal = float(saveNum(orgStore, ele))
    storeVal = float(saveNum(storeStore, ele))
    aa = ""
    bb = ""
    diff = 0
    if dir == "+":
        diff = orgVal - storeVal
    else:
        diff = storeVal - orgVal
    if diff - dist > 0.0:
        aa = '\033[42m'
        bb = '\033[0m'
    elif diff + dist < 0.0:
        aa = '\033[41m'
        bb = '\033[0m'
    if org:
        return(aa + form.format(orgVal) + bb)
    else:
        return(aa + form.format(storeVal) + bb)

def colorDiff(currStore, saveStore, ele, form):
    orgVal = form.format(saveNum(currStore, ele))
    storeVal = form.format(saveNum(saveStore, ele))
    aa = ""
    bb = ""
    if orgVal != storeVal:
        aa = '\033[41m'
        bb = '\033[0m'
    return(aa + orgVal + " (" + storeVal + ")" + bb)


def printPrimerTable(tars, mix, cur, sav):
    primConc = ["_100nM", "", "_750nM", "_100nM", "", "_750nM"]
    printConc = ["100nM", "250nM", "750nM", "100nM", "250nM", "750nM"]
    datOut = ["_mean", "_mean", "_mean", "_cv", "_cv", "_cv"]
    printDat = ["mean", "mean", "mean", "CV", "CV", "CV"]
    prec = [1.0, 1.0, 1.0, 0.01, 0.01, 0.01]
    prin = [0,0,0,5,5,5]
    res = "copies".ljust(20)
    for col in printDat:
        res += col.ljust(20)
    print(res)
    res = "exp: 1500".ljust(20)
    for col in printConc:
        res += col.ljust(20)
    print(res)
    colCurSum = np.zeros([len(primConc), len(tars)], dtype=np.float64)
    colSavSum = np.zeros([len(primConc), len(tars)], dtype=np.float64)
    colCurSum[:] = np.nan
    colSavSum[:] = np.nan
    roNum = 0
    for row in tars:
        res = row.ljust(20)
        count = 0
        for col in primConc:
            keyVal = row.replace(" ", "_") + col + mix + datOut[count]
            if keyVal in cur and keyVal in sav:
                colCurSum[count][roNum] = cur[keyVal]
                if keyVal not in sav:
                    continue
                colSavSum[count][roNum] = saveNum(sav, keyVal)
                direct = cur[keyVal] - saveNum(sav, keyVal)
                diff = abs(cur[keyVal] - saveNum(sav, keyVal))
                if diff > prec[count]:
                    if direct > 0.0:
                        res +=  '\033[44m'
                    else:
                        res +=  '\033[41m'
                if prin[count] == 0:
                    res += "{:8.0f}".format(cur[keyVal]) + " (" + "{:8.0f}".format(saveNum(sav, keyVal)) + ") "
                elif prin[count] == 1:
                    res += "{:8.1f}".format(cur[keyVal]) + " (" + "{:8.1f}".format(saveNum(sav, keyVal)) + ") "
                elif prin[count] == 2:
                    res += "{:8.2f}".format(cur[keyVal]) + " (" + "{:8.2f}".format(saveNum(sav, keyVal)) + ") "
                elif prin[count] == 3:
                    res += "{:8.3f}".format(cur[keyVal]) + " (" + "{:8.3f}".format(saveNum(sav, keyVal)) + ") "
                else:
                    res += "{:8.5f}".format(cur[keyVal]) + " (" + "{:8.5f}".format(saveNum(sav, keyVal)) + ") "
                if diff > prec[count]:
                    res +=  '\033[0m'
            else:
                res += "".ljust(20)
            count += 1
        print(res)
        roNum += 1
    res = "Mean".ljust(20)
    count = 0
    for col in primConc:
        keyVal = row.replace(" ", "_") + col + mix + datOut[count]
        meanCur = np.nanmean(colCurSum[count])
        if keyVal not in colSavSum:
            continue
        meanSav = np.nanmean(colSavSum[count])
        direct = meanCur - meanSav
        diff = abs(meanCur - meanSav)
        if diff > prec[count]:
            if direct > 0.0:
                res +=  '\033[44m'
            else:
                res +=  '\033[41m'
        if prin[count] == 0:
            res += "{:8.0f}".format(meanCur) + " (" + "{:8.0f}".format(meanSav) + ") "
        elif prin[count] == 1:
            res += "{:8.1f}".format(meanCur) + " (" + "{:8.1f}".format(meanSav) + ") "
        elif prin[count] == 2:
            res += "{:8.2f}".format(meanCur) + " (" + "{:8.2f}".format(meanSav) + ") "
        elif prin[count] == 3:
            res += "{:8.3f}".format(meanCur) + " (" + "{:8.3f}".format(meanSav) + ") "
        else:
            res += "{:8.5f}".format(meanCur) + " (" + "{:8.5f}".format(meanSav) + ") "
        if diff > prec[count]:
            res +=  '\033[0m'
        count += 1
    print(res)
    res = "CV".ljust(20)
    count = 0
    for col in primConc:
        keyVal = row.replace(" ", "_") + col + mix + datOut[count]
        if keyVal not in colSavSum:
            continue
        cvCur = np.nanstd(colCurSum[count]) / np.nanmean(colCurSum[count])
        cvSav = np.nanstd(colSavSum[count]) / np.nanmean(colSavSum[count])
        direct = cvCur - cvSav
        diff = abs(cvCur - cvSav)
        if diff > prec[count]:
            if direct > 0.0:
                res +=  '\033[44m'
            else:
                res +=  '\033[41m'
        res += "{:8.4f}".format(cvCur) + " (" + "{:8.4f}".format(cvSav) + ") "
        if diff > prec[count]:
            res +=  '\033[0m'
        count += 1
    print(res)


def printCVTable(tars, mix, cur, sav):
    print("".ljust(20) + "mean CV".ljust(20))
    for row in tars:
        res = row.ljust(20)
        keyVal = "dilution_" + row.replace(" ", "_") + mix + "_cv"
        if keyVal in cur:
            direct = cur[keyVal] - saveNum(sav, keyVal)
            diff = abs(cur[keyVal] - saveNum(sav, keyVal))
            if diff > 0.01:
                if direct > 0.0:
                    res +=  '\033[44m'
                else:
                    res +=  '\033[41m'
            res += "{:8.5f}".format(cur[keyVal]) + " (" + "{:8.5f}".format(saveNum(sav, keyVal)) + ") "
            if diff > 0.01:
                res +=  '\033[0m'
        else:
            res += "".ljust(20)
        print(res)


def printCVProbes(tars, mix, cur, sav):
    print("".ljust(30) + "mean CV".ljust(20))
    for row in tars:
        res = row.ljust(30)
        keyVal = "probes_" + row.replace(" ", "_") + mix + "_cv"
        if keyVal in cur:
            direct = cur[keyVal] - saveNum(sav, keyVal)
            diff = abs(cur[keyVal] - saveNum(sav, keyVal))
            if diff > 0.01:
                if direct > 0.0:
                    res +=  '\033[44m'
                else:
                    res +=  '\033[41m'
            res += "{:8.5f}".format(cur[keyVal]) + " (" + "{:8.5f}".format(saveNum(sav, keyVal)) + ") "
            if diff > 0.01:
                res +=  '\033[0m'
        else:
            res += "".ljust(20)
        print(res)


def printTable(spc, pref, rows, cols, prec, prin, cur, sav):
    res = "".ljust(spc)
    for col in cols:
        res += col.ljust(20)
    print(res)
    colCurSum = np.zeros([len(cols), len(rows)], dtype=np.float64)
    colSavSum = np.zeros([len(cols), len(rows)], dtype=np.float64)
    colCurSum[:] = np.nan
    colSavSum[:] = np.nan
    roNum = 0
    for row in rows:
        res = row.ljust(spc)
        count = 0
        for col in cols:
            keyVal = pref + row.replace(" ", "_") + "_" + col
            colCurSum[count][roNum] = cur[keyVal]
            colSavSum[count][roNum] = saveNum(sav, keyVal)
            direct = cur[keyVal] - saveNum(sav, keyVal)
            diff = abs(cur[keyVal] - saveNum(sav, keyVal))
            if diff > prec[count]:
                if direct > 0.0:
                    res +=  '\033[44m'
                else:
                    res +=  '\033[41m'
            if prin[count] == 0:
                res += "{:8.0f}".format(cur[keyVal]) + " (" + "{:8.0f}".format(saveNum(sav, keyVal)) + ") "
            elif prin[count] == 1:
                res += "{:8.1f}".format(cur[keyVal]) + " (" + "{:8.1f}".format(saveNum(sav, keyVal)) + ") "
            elif prin[count] == 2:
                res += "{:8.2f}".format(cur[keyVal]) + " (" + "{:8.2f}".format(saveNum(sav, keyVal)) + ") "
            elif prin[count] == 3:
                res += "{:8.3f}".format(cur[keyVal]) + " (" + "{:8.3f}".format(saveNum(sav, keyVal)) + ") "
            else:
                res += "{:8.5f}".format(cur[keyVal]) + " (" + "{:8.5f}".format(saveNum(sav, keyVal)) + ") "
            if diff > prec[count]:
                res +=  '\033[0m'
            count += 1
        print(res)
        roNum += 1
    res = "Mean".ljust(spc)
    count = 0
    for col in cols:
        meanCur = np.nanmean(colCurSum[count])
        meanSav = np.nanmean(colSavSum[count])
        direct = meanCur - meanSav
        diff = abs(meanCur - meanSav)
        if diff > prec[count]:
            if direct > 0.0:
                res +=  '\033[44m'
            else:
                res +=  '\033[41m'
        if prin[count] == 0:
            res += "{:8.0f}".format(meanCur) + " (" + "{:8.0f}".format(meanSav) + ") "
        elif prin[count] == 1:
            res += "{:8.1f}".format(meanCur) + " (" + "{:8.1f}".format(meanSav) + ") "
        elif prin[count] == 2:
            res += "{:8.2f}".format(meanCur) + " (" + "{:8.2f}".format(meanSav) + ") "
        elif prin[count] == 3:
            res += "{:8.3f}".format(meanCur) + " (" + "{:8.3f}".format(meanSav) + ") "
        else:
            res += "{:8.5f}".format(meanCur) + " (" + "{:8.5f}".format(meanSav) + ") "
        if diff > prec[count]:
            res +=  '\033[0m'
        count += 1
    print(res)
    res = "CV".ljust(spc)
    count = 0
    for col in cols:
        cvCur = np.nanstd(colCurSum[count]) / np.nanmean(colCurSum[count])
        cvSav = np.nanstd(colSavSum[count]) / np.nanmean(colSavSum[count])
        direct = cvCur - cvSav
        diff = abs(cvCur - cvSav)
        if diff > prec[count]:
            if direct > 0.0:
                res +=  '\033[44m'
            else:
                res +=  '\033[41m'
        res += "{:8.4f}".format(cvCur) + " (" + "{:8.4f}".format(cvSav) + ") "
        if diff > prec[count]:
            res +=  '\033[0m'
        count += 1
    print(res)
    if "curve_pcr_eff" in cols:
        res = "Diff cur dil".ljust(spc)
        std_ar = ["curve_pcr_eff", "dilution_pcr_eff"]
        sum_eff_cur = 0.0
        sum_eff_sav = 0.0
        num_eff = 0.0
        for row in rows:
            keyCurv = pref + row.replace(" ", "_") + "_curve_pcr_eff"
            keyDil = pref + row.replace(" ", "_") + "_dilution_pcr_eff"
            if keyCurv in cur and keyDil in cur:
                if keyCurv in sav and keyDil in sav:
                    sum_eff_cur += abs(cur[keyDil] - cur[keyCurv])
                    sum_eff_sav += abs(sav[keyDil] - sav[keyCurv])
                    num_eff += 1.0
        if num_eff > 0.0:
            eff_diff_cur = sum_eff_cur / num_eff
            eff_diff_sav = sum_eff_sav / num_eff
            res += "{:8.3f}".format(eff_diff_cur) + " (" + "{:8.3f}".format(eff_diff_sav) + ") "
            print(res)


laDa = {}
curDa = {}

with open(in_json, 'r') as openfile:
    laDa = json.load(openfile)

print("\n##############################\n### Remove Temporary Files ###\n##############################")
print('rm ' + os.path.join(parent_dir, "test/temp_*"))
os.system('rm ' + os.path.join(parent_dir, "test/temp_*"))

print("\n#######################\n### Technical Tests ###\n#######################")

print("Test: Fill Matrix Gaps:")
print("This RDML file has gaps in the matrix which need to be filled automatically.")
print("----------------------")

os.system('python3 ' + os.path.join(library_dir, "rdmlpython/rdml.py") +  
          ' -lrp ' + os.path.join(parent_dir, "experiments/technical/matrix_gaps.rdml") +
          ' -e "Experiment 1"' +
          ' -r "Run 1"' +
          ' --pcrEfficiencyExl 0.05' +
          ' --excludeNoPlateau' +
          ' --excludeEfficiency "outlier"' +
          ' --ignoreExclusion' +
          ' --saveRaw ' + os.path.join(parent_dir, "test/temp_matrix_gaps_raw_data.tsv") +
          ' --saveRawFixed ' + os.path.join(parent_dir, "test/temp_matrix_gaps_raw_fixed_data.tsv") +
          ' --saveBaslineCorr ' + os.path.join(parent_dir, "test/temp_matrix_gaps_baseline_corrected.tsv") +
          ' --saveResults ' + os.path.join(parent_dir, "test/temp_matrix_gaps_results.tsv"))

os.system('python3 ' + os.path.join(parent_dir, "test/diff_table.py") + ' ' + 
          os.path.join(parent_dir, "test/temp_matrix_gaps_raw_data.tsv") + ' ' + 
          os.path.join(parent_dir, "test/matrix_gaps_raw_data.tsv") + ' ' + 
          '"Matrix Gaps - test raw data" 20 N Y')
os.system('python3 ' + os.path.join(parent_dir, "test/diff_table.py") + ' ' + 
          os.path.join(parent_dir, "test/temp_matrix_gaps_raw_fixed_data.tsv") + ' ' + 
          os.path.join(parent_dir, "test/matrix_gaps_raw_fixed_data.tsv") + ' ' + 
          '"Matrix Gaps - test fixed raw data" 20 N Y')
os.system('python3 ' + os.path.join(parent_dir, "test/diff_table.py") + ' ' + 
          os.path.join(parent_dir, "test/temp_matrix_gaps_baseline_corrected.tsv") + ' ' + 
          os.path.join(parent_dir, "test/matrix_gaps_baseline_corrected.tsv") + ' ' + 
          '"Matrix Gaps - test baseline corrected data" 20 N Y')
os.system('python3 ' + os.path.join(parent_dir, "test/diff_table.py") + ' ' + 
          os.path.join(parent_dir, "test/temp_matrix_gaps_results.tsv") + ' ' + 
          os.path.join(parent_dir, "test/matrix_gaps_results.tsv") + ' ' + 
          '"Matrix Gaps - test results" 20 N Y')


print("\n################################\n### Test Machine and Volumes ###\n################################")
print("This test uses the dilution series with different volumes and machines. In principle all should have ")
print("the same TD0 and Ncopy should be as expected.")
print("----------------------")

rd = rdml.Rdml(rdml_vol_file)

expList = rd.experiments()
if len(expList) < 1:
    print("No experiments found!")
    sys.exit(0)
ww = open(out_vol_file, "w")
startLine = 0
linRegRes = {}
colTar = []
quant = []
for exp in expList:
    runList = exp.runs()
    if len(runList) < 1:
        print("No runs found!")
        sys.exit(0)
    for run in runList:
        if printExpRun:
            print("Experiment: " + exp["id"] + " Run: " + run["id"])
        res = run.webAppLinRegPCR(pcrEfficiencyExl=0.05, updateTargetEfficiency=True, updateRDML=True, excludeNoPlateau=True, excludeEfficiency="outlier", excludeInstableBaseline=True)
        resTab = json.loads(res["LinRegPCR_Result_Table"])
        for tabRow in range(0, len(resTab)):
            if startLine == 0:
                ww.write("Experiment\tRun\t")
                startLine = 1
            else:
                ww.write(exp["id"] + "\t" + run["id"] + "\t")
            for tabCol in range(0, len(resTab[tabRow])):
                outCell = str(resTab[tabRow][tabCol]).replace("\t", ";")
                if tabCol < len(resTab[tabRow]) - 1:
                    ww.write(outCell + "\t")
                else:
                    ww.write(outCell + "\n")
            if startLine == 1:
                if resTab[tabRow][rar_sample_type] in ["unkn", "std"]:
                    if float(resTab[tabRow][rar_Ncopy]) > 5.0:
                        if exp["id"] not in linRegRes:
                            linRegRes[exp["id"]] = {}
                        if run["id"] not in linRegRes[exp["id"]]:
                            linRegRes[exp["id"]][run["id"]] = {}
                        if resTab[tabRow][rar_tar] not in linRegRes[exp["id"]][run["id"]]:
                            linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]] = {}
                        if resTab[tabRow][rar_tar] not in colTar:
                            colTar.append(resTab[tabRow][rar_tar])
                        if resTab[tabRow][rar_sample] not in linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]:  # Sample
                            linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]][resTab[tabRow][rar_sample]] = {}
                            linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]][resTab[tabRow][rar_sample]]["Ncopy"] = []
                            linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]][resTab[tabRow][rar_sample]]["TD0"] = []
                        
                        linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]][resTab[tabRow][rar_sample]]["Ncopy"].append(resTab[tabRow][rar_Ncopy])  # Ncopy
                        linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]][resTab[tabRow][rar_sample]]["TD0"].append(resTab[tabRow][rar_TD0])  # TD0

        for tar in colTar:
            for sam in linRegRes[exp["id"]][run["id"]][tar]:
                linRegRes[exp["id"]][run["id"]][tar][sam]["Ncopy mean"] = np.mean(linRegRes[exp["id"]][run["id"]][tar][sam]["Ncopy"])
                linRegRes[exp["id"]][run["id"]][tar][sam]["TD0 mean"] = np.mean(linRegRes[exp["id"]][run["id"]][tar][sam]["TD0"])

ww.close()
rd.save("temp_volume.rdml")

curDa["Test_Vol_AMC_384_05_TD0"] = linRegRes["AMC Amsterdam NL - 384 Wells"]["Run 1"]["FSTL1"]["M1 gDNA 0.5ng/ul - 5ul"]["TD0 mean"]
curDa["Test_Vol_AMC_384_10_TD0"] = linRegRes["AMC Amsterdam NL - 384 Wells"]["Run 1"]["FSTL1"]["M1 gDNA 0.5ng/ul STD"]["TD0 mean"]
curDa["Test_Vol_AMC_384_20_TD0"] = linRegRes["AMC Amsterdam NL - 384 Wells"]["Run 1"]["FSTL1"]["M1 gDNA 0.5ng/ul - 20ul"]["TD0 mean"]
curDa["Test_Vol_EMBL_384_05_TD0"] = linRegRes["EMBL Heidelberg DE - 384 Wells"]["Run 1"]["FSTL1"]["M1 gDNA 0.5ng/ul - 5ul"]["TD0 mean"]
curDa["Test_Vol_EMBL_384_10_TD0"] = linRegRes["EMBL Heidelberg DE - 384 Wells"]["Run 1"]["FSTL1-STD"]["M1 gDNA 0.5ng/ul STD"]["TD0 mean"]
curDa["Test_Vol_EMBL_384_20_TD0"] = linRegRes["EMBL Heidelberg DE - 384 Wells"]["Run 1"]["FSTL1"]["M1 gDNA 0.5ng/ul - 20ul"]["TD0 mean"]
curDa["Test_Vol_EMBL_96_10_TD0"] = linRegRes["EMBL Heidelberg DE - 96 Wells"]["Run 1"]["FSTL1"]["M1 gDNA 0.5ng/ul - 10ul"]["TD0 mean"]
curDa["Test_Vol_EMBL_96_20_TD0"] = linRegRes["EMBL Heidelberg DE - 96 Wells"]["Run 1"]["FSTL1-STD"]["M1 gDNA 0.5ng/ul STD"]["TD0 mean"]
curDa["Test_Vol_EMBL_96_40_TD0"] = linRegRes["EMBL Heidelberg DE - 96 Wells"]["Run 1"]["FSTL1"]["M1 gDNA 0.5ng/ul - 40ul"]["TD0 mean"]

print("\nAMC Amsterdam NL - 384 Wells")
print("  5ul - 0.5ng/ul - TD0: " + colorDiff(curDa, laDa, "Test_Vol_AMC_384_05_TD0", "{:6.2f}"))
print(" 10ul - 0.5ng/ul - TD0: " + colorDiff(curDa, laDa, "Test_Vol_AMC_384_10_TD0", "{:6.2f}"))
print(" 20ul - 0.5ng/ul - TD0: " + colorDiff(curDa, laDa, "Test_Vol_AMC_384_20_TD0", "{:6.2f}"))
print("\nEMBL Heidelberg DE - 384 Wells")
print("  5ul - 0.5ng/ul - TD0: " + colorDiff(curDa, laDa, "Test_Vol_EMBL_384_05_TD0", "{:6.2f}"))
print(" 10ul - 0.5ng/ul - TD0: " + colorDiff(curDa, laDa, "Test_Vol_EMBL_384_10_TD0", "{:6.2f}"))
print(" 20ul - 0.5ng/ul - TD0: " + colorDiff(curDa, laDa, "Test_Vol_EMBL_384_20_TD0", "{:6.2f}"))
print("\nEMBL Heidelberg DE - 96 Wells")
print(" 10ul - 0.5ng/ul - TD0: " + colorDiff(curDa, laDa, "Test_Vol_EMBL_96_10_TD0", "{:6.2f}"))
print(" 20ul - 0.5ng/ul - TD0: " + colorDiff(curDa, laDa, "Test_Vol_EMBL_96_20_TD0", "{:6.2f}"))
print(" 40ul - 0.5ng/ul - TD0: " + colorDiff(curDa, laDa, "Test_Vol_EMBL_96_40_TD0", "{:6.2f}"))

curDa["Test_Vol_AMC_384_05"] = linRegRes["AMC Amsterdam NL - 384 Wells"]["Run 1"]["FSTL1"]["M1 gDNA 0.5ng/ul - 5ul"]["Ncopy mean"]
curDa["Test_Vol_AMC_384_10"] = linRegRes["AMC Amsterdam NL - 384 Wells"]["Run 1"]["FSTL1"]["M1 gDNA 0.5ng/ul STD"]["Ncopy mean"]
curDa["Test_Vol_AMC_384_20"] = linRegRes["AMC Amsterdam NL - 384 Wells"]["Run 1"]["FSTL1"]["M1 gDNA 0.5ng/ul - 20ul"]["Ncopy mean"]
curDa["Test_Vol_EMBL_384_05"] = linRegRes["EMBL Heidelberg DE - 384 Wells"]["Run 1"]["FSTL1"]["M1 gDNA 0.5ng/ul - 5ul"]["Ncopy mean"]
curDa["Test_Vol_EMBL_384_10"] = linRegRes["EMBL Heidelberg DE - 384 Wells"]["Run 1"]["FSTL1-STD"]["M1 gDNA 0.5ng/ul STD"]["Ncopy mean"]
curDa["Test_Vol_EMBL_384_20"] = linRegRes["EMBL Heidelberg DE - 384 Wells"]["Run 1"]["FSTL1"]["M1 gDNA 0.5ng/ul - 20ul"]["Ncopy mean"]
curDa["Test_Vol_EMBL_96_10"] = linRegRes["EMBL Heidelberg DE - 96 Wells"]["Run 1"]["FSTL1"]["M1 gDNA 0.5ng/ul - 10ul"]["Ncopy mean"]
curDa["Test_Vol_EMBL_96_20"] = linRegRes["EMBL Heidelberg DE - 96 Wells"]["Run 1"]["FSTL1-STD"]["M1 gDNA 0.5ng/ul STD"]["Ncopy mean"]
curDa["Test_Vol_EMBL_96_40"] = linRegRes["EMBL Heidelberg DE - 96 Wells"]["Run 1"]["FSTL1"]["M1 gDNA 0.5ng/ul - 40ul"]["Ncopy mean"]

print("\nAMC Amsterdam NL - 384 Wells")
print("  5ul - 0.5ng/ul -  750 copies: " + colorDiff(curDa, laDa, "Test_Vol_AMC_384_05", "{:6.1f}"))
print(" 10ul - 0.5ng/ul - 1500 copies: " + colorDiff(curDa, laDa, "Test_Vol_AMC_384_10", "{:6.1f}"))
print(" 20ul - 0.5ng/ul - 3000 copies: " + colorDiff(curDa, laDa, "Test_Vol_AMC_384_20", "{:6.1f}"))
print("\nEMBL Heidelberg DE - 384 Wells")
print("  5ul - 0.5ng/ul -  750 copies: " + colorDiff(curDa, laDa, "Test_Vol_EMBL_384_05", "{:6.1f}"))
print(" 10ul - 0.5ng/ul - 1500 copies: " + colorDiff(curDa, laDa, "Test_Vol_EMBL_384_10", "{:6.1f}"))
print(" 20ul - 0.5ng/ul - 3000 copies: " + colorDiff(curDa, laDa, "Test_Vol_EMBL_384_20", "{:6.1f}"))
print("\nEMBL Heidelberg DE - 96 Wells")
print(" 10ul - 0.5ng/ul - 1500 copies: " + colorDiff(curDa, laDa, "Test_Vol_EMBL_96_10", "{:6.1f}"))
print(" 20ul - 0.5ng/ul - 3000 copies: " + colorDiff(curDa, laDa, "Test_Vol_EMBL_96_20", "{:6.1f}"))
print(" 40ul - 0.5ng/ul - 6000 copies: " + colorDiff(curDa, laDa, "Test_Vol_EMBL_96_40", "{:6.1f}"))

curDa["Test_Vol_AMC_384_c0.2"] = linRegRes["AMC Amsterdam NL - 384 Wells"]["Run 1"]["FSTL1"]["gDNA 0.2ng/ul"]["Ncopy mean"]
curDa["Test_Vol_AMC_384_c0.5"] = linRegRes["AMC Amsterdam NL - 384 Wells"]["Run 1"]["FSTL1"]["M1 gDNA 0.5ng/ul STD"]["Ncopy mean"]
curDa["Test_Vol_AMC_384_c1"] = linRegRes["AMC Amsterdam NL - 384 Wells"]["Run 1"]["FSTL1"]["gDNA 1ng/ul"]["Ncopy mean"]
curDa["Test_Vol_AMC_384_c2"] = linRegRes["AMC Amsterdam NL - 384 Wells"]["Run 1"]["FSTL1"]["gDNA 2ng/ul"]["Ncopy mean"]
curDa["Test_Vol_AMC_384_c4"] = linRegRes["AMC Amsterdam NL - 384 Wells"]["Run 1"]["FSTL1"]["gDNA 4ng/ul"]["Ncopy mean"]
curDa["Test_Vol_AMC_384_c8"] = linRegRes["AMC Amsterdam NL - 384 Wells"]["Run 1"]["FSTL1"]["gDNA 8ng/ul"]["Ncopy mean"]
curDa["Test_Vol_EMBL_384_c0.2"] = linRegRes["EMBL Heidelberg DE - 384 Wells"]["Run 1"]["FSTL1"]["gDNA 0.2ng/ul"]["Ncopy mean"]
curDa["Test_Vol_EMBL_384_c0.5"] = linRegRes["EMBL Heidelberg DE - 384 Wells"]["Run 1"]["FSTL1-STD"]["M1 gDNA 0.5ng/ul STD"]["Ncopy mean"]
curDa["Test_Vol_EMBL_384_c1"] = linRegRes["EMBL Heidelberg DE - 384 Wells"]["Run 1"]["FSTL1"]["gDNA 1ng/ul"]["Ncopy mean"]
curDa["Test_Vol_EMBL_384_c2"] = linRegRes["EMBL Heidelberg DE - 384 Wells"]["Run 1"]["FSTL1"]["gDNA 2ng/ul"]["Ncopy mean"]
curDa["Test_Vol_EMBL_384_c4"] = linRegRes["EMBL Heidelberg DE - 384 Wells"]["Run 1"]["FSTL1"]["gDNA 4ng/ul"]["Ncopy mean"]
curDa["Test_Vol_EMBL_96_c0.2"] = linRegRes["EMBL Heidelberg DE - 96 Wells"]["Run 1"]["FSTL1"]["gDNA 0.2ng/ul"]["Ncopy mean"]
curDa["Test_Vol_EMBL_96_c0.5"] = linRegRes["EMBL Heidelberg DE - 96 Wells"]["Run 1"]["FSTL1-STD"]["M1 gDNA 0.5ng/ul STD"]["Ncopy mean"]
curDa["Test_Vol_EMBL_96_c1"] = linRegRes["EMBL Heidelberg DE - 96 Wells"]["Run 1"]["FSTL1"]["gDNA 1ng/ul"]["Ncopy mean"]
curDa["Test_Vol_EMBL_96_c2"] = linRegRes["EMBL Heidelberg DE - 96 Wells"]["Run 1"]["FSTL1"]["gDNA 2ng/ul"]["Ncopy mean"]
curDa["Test_Vol_EMBL_96_c4"] = linRegRes["EMBL Heidelberg DE - 96 Wells"]["Run 1"]["FSTL1"]["gDNA 4ng/ul"]["Ncopy mean"]

print("\nAMC Amsterdam NL - 384 Wells")
print(" 10ul - 0.2ng/ul -   600 copies: " + colorDiff(curDa, laDa, "Test_Vol_AMC_384_c0.2", "{:8.1f}"))
print(" 10ul - 0.5ng/ul -  1500 copies: " + colorDiff(curDa, laDa, "Test_Vol_AMC_384_c0.5", "{:8.1f}"))
print(" 10ul - 1.0ng/ul -  3000 copies: " + colorDiff(curDa, laDa, "Test_Vol_AMC_384_c1", "{:8.1f}"))
print(" 10ul - 2.0ng/ul -  6000 copies: " + colorDiff(curDa, laDa, "Test_Vol_AMC_384_c2", "{:8.1f}"))
print(" 10ul - 4.0ng/ul - 12000 copies: " + colorDiff(curDa, laDa, "Test_Vol_AMC_384_c4", "{:8.1f}"))
print(" 10ul - 8.0ng/ul - 24000 copies: " + colorDiff(curDa, laDa, "Test_Vol_AMC_384_c8", "{:8.1f}"))
print("\nEMBL Heidelberg DE - 384 Wells")
print(" 10ul - 0.2ng/ul -   600 copies: " + colorDiff(curDa, laDa, "Test_Vol_EMBL_384_c0.2", "{:8.1f}"))
print(" 10ul - 0.5ng/ul -  1500 copies: " + colorDiff(curDa, laDa, "Test_Vol_EMBL_384_c0.2", "{:8.1f}"))
print(" 10ul - 1.0ng/ul -  3000 copies: " + colorDiff(curDa, laDa, "Test_Vol_EMBL_384_c1", "{:8.1f}"))
print(" 10ul - 2.0ng/ul -  6000 copies: " + colorDiff(curDa, laDa, "Test_Vol_EMBL_384_c2", "{:8.1f}"))
print(" 10ul - 4.0ng/ul - 12000 copies: " + colorDiff(curDa, laDa, "Test_Vol_EMBL_384_c4", "{:8.1f}"))
print("\nEMBL Heidelberg DE - 96 Wells")
print(" 20ul - 0.2ng/ul -  1200 copies: " + colorDiff(curDa, laDa, "Test_Vol_EMBL_96_c0.2", "{:8.1f}"))
print(" 20ul - 0.5ng/ul -  3000 copies: " + colorDiff(curDa, laDa, "Test_Vol_EMBL_96_c0.2", "{:8.1f}"))
print(" 20ul - 1.0ng/ul -  6000 copies: " + colorDiff(curDa, laDa, "Test_Vol_EMBL_96_c1", "{:8.1f}"))
print(" 20ul - 2.0ng/ul - 12000 copies: " + colorDiff(curDa, laDa, "Test_Vol_EMBL_96_c2", "{:8.1f}"))
print(" 20ul - 4.0ng/ul - 24000 copies: " + colorDiff(curDa, laDa, "Test_Vol_EMBL_96_c4", "{:8.1f}"))

# Time the test
startTime = time.time()
rd = rdml.Rdml(rdml_amp_file)

expList = rd.experiments()
if len(expList) < 1:
    print("No experiments found in " + rdml_amp_file + "!")
    sys.exit(0)
ww = open(out_amp_file, "w")
startLine = 0
reactionDataTrue = 0
reactionDataFalse = 0
linRegRes = {}
for exp in expList:
    if exp["id"] == "Primer Test - EMBL":
        print("\n###################\n### Primer Test ###\n###################")
        print("This test uses the primer test data. In principle all targets should have a high PCReff and a ")
        print("Ncopy of 1200. The TD0 should decrease with amplicon length.")
        print("----------------------")

    if exp["id"] == "Primer Conc - AMC":
        print("\n#############################\n### Test Primer Dilutions ###\n#############################")
    runList = exp.runs()
    if len(runList) < 1:
        print("No runs found!")
        sys.exit(0)
    for run in runList:
        if printExpRun:
            print("Experiment: " + exp["id"] + " Run: " + run["id"])
        res = run.webAppLinRegPCR(pcrEfficiencyExl=0.05, updateTargetEfficiency=False, updateRDML=True, excludeNoPlateau=True, excludeEfficiency="outlier", excludeInstableBaseline=True)
        resTab = json.loads(res["LinRegPCR_Result_Table"])
        for tabRow in range(startLine, len(resTab)):
            if startLine == 0:
                ww.write("Experiment\tRun\t")
                startLine = 1
            else:
                ww.write(exp["id"] + "\t" + run["id"] + "\t")
            for tabCol in range(0, len(resTab[tabRow])):
                outCell = str(resTab[tabRow][tabCol]).replace("\t", ";")
                if tabCol < len(resTab[tabRow]) - 1:
                    ww.write(outCell + "\t")
                else:
                    ww.write(outCell + "\n")
            if startLine == 1:
                if resTab[tabRow][rar_sample_type] in ["unkn", "std"]:
                    if float(resTab[tabRow][rar_Ncopy]) > 5.0:
                        if exp["id"] not in linRegRes:
                            linRegRes[exp["id"]] = {}
                        if run["id"] not in linRegRes[exp["id"]]:
                            linRegRes[exp["id"]][run["id"]] = {}
                        if resTab[tabRow][rar_tar] not in linRegRes[exp["id"]][run["id"]]:
                            linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]] = {}
                            linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["PCR Eff indiv"] = []
                            linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["Ncopy indiv"] = []
                            linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["Ncopy"] = []
                            linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["TD0"] = []

                        linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["PCR Eff indiv"].append(resTab[tabRow][rar_indiv_PCR_eff])
                        linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["Ncopy indiv"].append(resTab[tabRow][rar_indiv_Ncopy])
                        linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["Ncopy"].append(resTab[tabRow][rar_Ncopy])
                        linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["TD0"].append(resTab[tabRow][rar_TD0])
                        linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["PCReff"] = resTab[tabRow][rar_PCR_eff]

                        reactionDataTrue += 1
                    else:
                        reactionDataFalse += 1

    if exp["id"] == "Primer Test - EMBL":
        prim_test_key = "Primer_Test_"
        prim_test_tar = ["FSTL_1_A_047", "FSTL_1_B_042", "FSTL_1_C_040", "FSTL_1_D_105", "FSTL_1_E_097", "FSTL_1_F_109", "FSTL_1_G_211",
                         "FSTL_1_H_201", "FSTL_1_I_204", "FSTL_1_K_219", "FSTL_1_*_259", "FSTL_1_*_259 AMC", "FSTL_1_L_412", "FSTL_1_M_398",
                         "FSTL_1_N_386", "FSTL_1_O_417", "FSTL_1_P_820", "FSTL_1_R_790", "FSTL_1_S_778", "FSTL_1_T_770"]
        for tar in prim_test_tar:
            curDa[prim_test_key + "PCReff_" + tar] = linRegRes["Primer Test - EMBL"]["Primer Test"][tar]["PCReff"]
            curDa[prim_test_key + "TD0_" + tar] = np.mean(linRegRes["Primer Test - EMBL"]["Primer Test"][tar]["TD0"])
            curDa[prim_test_key + "Ncopy_" + tar] = np.mean(linRegRes["Primer Test - EMBL"]["Primer Test"][tar]["Ncopy"])
        print("Target                TD0    (TD0   )    PCReff (PCReff)     Ncopy  (Ncopy  )")
        for tar in prim_test_tar:
            print(tar.ljust(22) + colorDiff(curDa, laDa, prim_test_key + "TD0_" + tar, "{:6.3f}") + "    "
                                + colorDiff(curDa, laDa, prim_test_key + "PCReff_" + tar, "{:6.3f}") + "     "
                                + colorDiff(curDa, laDa, prim_test_key + "Ncopy_" + tar, "{:6.0f}"))

    if exp["id"] == "Primer Conc - AMC":
        quant = exp.quantify()
        tars = ["FSTL_1_A_047", "FSTL_1_B_042", "FSTL_1_D_105", "FSTL_1_E_097", "FSTL_1_F_109", "FSTL_1_H_201",
                "FSTL_1_I_204", "FSTL_1_K_219", "FSTL_1_*_259", "FSTL_1_L_412", "FSTL_1_M_398", "FSTL_1_O_417", "FSTL_1_P_820"]
        tar_roche = []
        tar_sensi = []
        tar_lc    = []
        for tar in tars:
            tar_roche.append(tar + " 100nM")
            tar_sensi.append(tar + " 100nM SF")
            tar_lc.append(tar + " 100nM LC")
        for tar in tars:
            tar_roche.append(tar)
            tar_sensi.append(tar + " SF")
            tar_lc.append(tar + " LC")
        for tar in tars:
            tar_roche.append(tar + " 750nM")
            tar_sensi.append(tar + " 750nM SF")
            tar_lc.append(tar + " 750nM LC")
        for ele in ["FSTL_1_O_417 100nM", "FSTL_1_O_417", "FSTL_1_O_417 750nM", "FSTL_1_M_398 100nM", "FSTL_1_P_820 100nM","FSTL_1_M_398", "FSTL_1_P_820","FSTL_1_M_398 750nM", "FSTL_1_P_820 750nM"]:
            tar_roche.remove(ele)
        for ele in ["FSTL_1_D_105 100nM SF", "FSTL_1_M_398 100nM SF", "FSTL_1_O_417 100nM SF", 
                    "FSTL_1_D_105 SF", "FSTL_1_M_398 SF", "FSTL_1_O_417 SF", "FSTL_1_D_105 750nM SF", 
                    "FSTL_1_M_398 750nM SF", "FSTL_1_O_417 750nM SF"]:
            tar_sensi.remove(ele)
        for ele in ["FSTL_1_A_047 100nM LC", "FSTL_1_B_042 100nM LC", "FSTL_1_D_105 100nM LC", "FSTL_1_E_097 100nM LC", 
                    "FSTL_1_*_259 100nM LC", "FSTL_1_L_412 100nM LC", "FSTL_1_A_047 LC", "FSTL_1_B_042 LC", 
                    "FSTL_1_D_105 LC", "FSTL_1_E_097 LC", "FSTL_1_*_259 LC", "FSTL_1_L_412 LC", "FSTL_1_A_047 750nM LC", 
                    "FSTL_1_B_042 750nM LC", "FSTL_1_D_105 750nM LC", "FSTL_1_E_097 750nM LC", "FSTL_1_*_259 750nM LC", 
                    "FSTL_1_L_412 750nM LC"]:
            tar_lc.remove(ele)


        for tar in tar_roche:
            if tar not in quant["tec_data"]["genomic DNA 4ng"]:
                print("Error: unknown target " + tar)
            keyVal = tar.replace(" ", "_") + "_mean"
            curDa[keyVal] = quant["tec_data"]["genomic DNA 4ng"][tar]["Ncopy_mean"]
            keyVal = tar.replace(" ", "_") + "_cv"
            curDa[keyVal] = quant["tec_data"]["genomic DNA 4ng"][tar]["Ncopy_cv"]
        for tar in tar_sensi:
            if tar not in quant["tec_data"]["genomic DNA 4ng"]:
                print("Error: unknown target " + tar)
            keyVal = tar.replace(" ", "_") + "_mean"
            curDa[keyVal] = quant["tec_data"]["genomic DNA 4ng"][tar]["Ncopy_mean"]
            keyVal = tar.replace(" ", "_") + "_cv"
            curDa[keyVal] = quant["tec_data"]["genomic DNA 4ng"][tar]["Ncopy_cv"]   
        for tar in tar_lc:
            if tar not in quant["tec_data"]["genomic DNA 4ng"]:
                print("Error: unknown target " + tar)
            keyVal = tar.replace(" ", "_") + "_mean"
            curDa[keyVal] = quant["tec_data"]["genomic DNA 4ng"][tar]["Ncopy_mean"]
            keyVal = tar.replace(" ", "_") + "_cv"
            curDa[keyVal] = quant["tec_data"]["genomic DNA 4ng"][tar]["Ncopy_cv"]   


        print("\n####################\n### Ncopy and CV ###\n####################")
        print("This test uses the primer dilution data. In principle all targets should have a Ncopy of 1200. The CV ")
        print("of the seven reactions should be as low as possible.")
        print("----------------------")

        print("Roche:")
        printPrimerTable(tars, "", curDa, laDa)
        print("\nSensi:")
        printPrimerTable(tars, "_SF", curDa, laDa)
        print("\nLC-Green:")
        printPrimerTable(tars, "_LC", curDa, laDa)

        print("\n#####################\n### indiv PCR eff ###\n#####################")
        print("This test uses the primer dilution data. The individual PCR efficiency and the mean ")
        print("of the seven reactions is calculated and the CV of the efficiencies.")
        print("----------------------")
        print("\nRoche PCR efficiencies:                                                                           Mean                      CV")
        doTD0 = False
        xe = exp["id"]
        xr = "P1 - Primer Conc - Roche"
        for tar in tar_roche:
            if tar not in linRegRes[xe][xr]:
                continue
            keyValmean = tar.replace(" ", "_") + "_PCR_eff_mean"
            keyValstd = tar.replace(" ", "_") + "_PCR_eff_cv"
            lin = tar.ljust(22)
            linTD0 = tar.ljust(22)
            linNcopy = tar.ljust(22)
            cur_mean = np.mean(linRegRes[xe][xr][tar]["PCR Eff indiv"])
            curDa[keyValmean] = cur_mean
            cur_cv = np.std(linRegRes[xe][xr][tar]["PCR Eff indiv"]) / cur_mean
            curDa[keyValstd] = cur_cv
            for pos in range(0, 7):
                if pos < len(linRegRes[xe][xr][tar]["TD0"]):
                    linTD0 +=  "{:10.2f}".format(linRegRes[xe][xr][tar]["TD0"][pos])
                else:
                    linTD0 += "          "
                if pos < len(linRegRes[xe][xr][tar]["Ncopy indiv"]):
                    linNcopy +=  "{:10.2f}".format(linRegRes[xe][xr][tar]["Ncopy indiv"][pos])
                else:
                    linNcopy += "          "
                if pos < len(linRegRes[xe][xr][tar]["PCR Eff indiv"]):
                    lin +=  "{:10.6f}".format(linRegRes[xe][xr][tar]["PCR Eff indiv"][pos])
                else:
                    lin += "          "
            if keyValmean in laDa:
                if abs(cur_mean - laDa[keyValmean]) > 0.00001:
                    if cur_mean - laDa[keyValmean] > 0.0:
                        lin +=  "\033[41m"
                    else:
                        lin +=  "\033[44m"
            lin +=  "  {:10.6f}".format(cur_mean)
            if keyValmean in laDa:
                lin +=  " ({:10.6f})".format(laDa[keyValmean])
                if abs(cur_mean - laDa[keyValmean]) > 0.00001:
                    lin +=  "\033[0m"
            else:
                lin += " (         )"
            if keyValstd in laDa:
                if abs(cur_cv - laDa[keyValstd]) > 0.00001:
                    if cur_cv - laDa[keyValstd] > 0.0:
                        lin +=  "\033[41m"
                    else:
                        lin +=  "\033[44m"
            lin +=  "  {:10.6f}".format(cur_cv)
            if keyValstd in laDa:
                lin +=  " ({:10.6f})".format(laDa[keyValstd])
                if abs(cur_cv - laDa[keyValstd]) > 0.00001:
                    lin +=  "\033[0m"
            else:
                lin += " (         )"
            if doTD0:
                print(linTD0)
                print(linNcopy)
            print(lin)

        print("\nSensi PCR efficiencies:                                                                           Mean                      CV")
        xr = "P5 - Primer Conc - SensiFast"
        for tar in tar_sensi:
            if tar not in linRegRes[xe][xr]:
                continue
            keyValmean = tar.replace(" ", "_") + "_PCR_eff_mean"
            keyValstd = tar.replace(" ", "_") + "_PCR_eff_cv"
            lin = tar.ljust(22)
            cur_mean = np.mean(linRegRes[xe][xr][tar]["PCR Eff indiv"])
            curDa[keyValmean] = cur_mean
            cur_cv = np.std(linRegRes[xe][xr][tar]["PCR Eff indiv"]) / cur_mean
            curDa[keyValstd] = cur_cv
            for pos in range(0, 7):
                if pos < len(linRegRes[xe][xr][tar]["PCR Eff indiv"]):
                    lin +=  "{:10.6f}".format(linRegRes[xe][xr][tar]["PCR Eff indiv"][pos])
                else:
                    lin += "          "
            if keyValmean in laDa:
                if abs(cur_mean - laDa[keyValmean]) > 0.00001:
                    if cur_mean - laDa[keyValmean] > 0.0:
                        lin +=  "\033[41m"
                    else:
                        lin +=  "\033[44m"
            lin +=  "  {:10.6f}".format(cur_mean)
            if keyValmean in laDa:
                lin +=  " ({:10.6f})".format(laDa[keyValmean])
                if abs(cur_mean - laDa[keyValmean]) > 0.00001:
                    lin +=  "\033[0m"
            else:
                lin += " (         )"
            if keyValstd in laDa:
                if abs(cur_cv - laDa[keyValstd]) > 0.00001:
                    if cur_cv - laDa[keyValstd] > 0.0:
                        lin +=  "\033[41m"
                    else:
                        lin +=  "\033[44m"
            lin +=  "  {:10.6f}".format(cur_cv)
            if keyValstd in laDa:
                lin +=  " ({:10.6f})".format(laDa[keyValstd])
                if abs(cur_cv - laDa[keyValstd]) > 0.00001:
                    lin +=  "\033[0m"
            else:
                lin += " (         )"
            print(lin)

        print("\nLC-Green PCR efficiencies:                                                                        Mean                      CV")
        xr = "P2 - Primer Conc - LC Green"
        for tar in tar_lc:
            if tar not in linRegRes[xe][xr]:
                continue
            keyValmean = tar.replace(" ", "_") + "_PCR_eff_mean"
            keyValstd = tar.replace(" ", "_") + "_PCR_eff_cv"
            lin = tar.ljust(22)
            cur_mean = np.mean(linRegRes[xe][xr][tar]["PCR Eff indiv"])
            curDa[keyValmean] = cur_mean
            cur_cv = np.std(linRegRes[xe][xr][tar]["PCR Eff indiv"]) / cur_mean
            curDa[keyValstd] = cur_cv
            for pos in range(0, 7):
                if pos < len(linRegRes[xe][xr][tar]["PCR Eff indiv"]):
                    lin +=  "{:10.6f}".format(linRegRes[xe][xr][tar]["PCR Eff indiv"][pos])
                else:
                    lin += "          "
            if keyValmean in laDa:
                if abs(cur_mean - laDa[keyValmean]) > 0.00001:
                    if cur_mean - laDa[keyValmean] > 0.0:
                        lin +=  "\033[41m"
                    else:
                        lin +=  "\033[44m"
            lin +=  "  {:10.6f}".format(cur_mean)
            if keyValmean in laDa:
                lin +=  " ({:10.6f})".format(laDa[keyValmean])
                if abs(cur_mean - laDa[keyValmean]) > 0.00001:
                    lin +=  "\033[0m"
            else:
                lin += " (         )"
            if keyValstd in laDa:
                if abs(cur_cv - laDa[keyValstd]) > 0.00001:
                    if cur_cv - laDa[keyValstd] > 0.0:
                        lin +=  "\033[41m"
                    else:
                        lin +=  "\033[44m"
            lin +=  "  {:10.6f}".format(cur_cv)
            if keyValstd in laDa:
                lin +=  " ({:10.6f})".format(laDa[keyValstd])
                if abs(cur_cv - laDa[keyValstd]) > 0.00001:
                    lin +=  "\033[0m"
            else:
                lin += " (         )"
            print(lin)

    if exp["id"] == "DNA Dilution - AMC":
        quant = exp.quantify()
        head_var = ["curve_pcr_eff", "dilution_pcr_eff", "expected_max", "mean_max", "expected_min", "mean_min", 
                    "expected_ratio", "mean_ratio", "slope_bias", "correlation_R", "linearity", 
                    "reproducibility", "detectable_diff"]
        head_a = ["curve_pcr_eff", "dilution_pcr_eff", "expected_max", "mean_max", "expected_min", "mean_min", 
                    "expected_ratio", "mean_ratio"]
        head_b = ["slope_bias", "correlation_R", "linearity", "reproducibility", "detectable_diff"]
        prec_a = [0.001, 0.001, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        prec_b = [0.0001, 0.00001, 0.00001, 0.00001, 0.0001]
        print_a = [3, 3, 0, 0, 0, 0, 0, 0]
        print_b = [5, 5, 5, 5, 5]
        tar_roche = ["FSTL_1_A_047", "FSTL_1_B_042", "FSTL_1_D_105", "FSTL_1_E_097", "FSTL_1_F_109", "FSTL_1_H_201", 
                     "FSTL_1_I_204", "FSTL_1_K_219", "FSTL_1_*_259", "FSTL_1_L_412"] # , "FSTL_1_M_398", "FSTL_1_P_820"]
        tar_sensi = ["FSTL_1_A_047 SF", "FSTL_1_D_105 SF", "FSTL_1_E_097 SF", "FSTL_1_F_109 SF", "FSTL_1_H_201 SF",
                     "FSTL_1_I_204 SF", "FSTL_1_K_219 SF", "FSTL_1_*_259 SF", "FSTL_1_L_412 SF", "FSTL_1_M_398 SF",
                     "FSTL_1_O_417 SF", "FSTL_1_P_820 SF", "FSTL_1_T_770 SF"]
        tar_lc    = ["FSTL_1_F_109 LC", "FSTL_1_H_201 LC", "FSTL_1_I_204 LC", "FSTL_1_K_219 LC", "FSTL_1_M_398 LC",
                     "FSTL_1_O_417 LC", "FSTL_1_P_820 LC"]

        for tar in tar_roche:
            for var in head_var:
                keyVal = tar.replace(" ", "_") + "_" + var
                curDa[keyVal] = quant["dil_std"][tar][var]
        for tar in tar_sensi:
            for var in head_var:
                keyVal = tar.replace(" ", "_") + "_" + var
                curDa[keyVal] = quant["dil_std"][tar][var]
        for tar in tar_lc:
            for var in head_var:
                keyVal = tar.replace(" ", "_") + "_" + var
                curDa[keyVal] = quant["dil_std"][tar][var]

        samDna = ["genomic DNA 20ng", "genomic DNA 7ng", "genomic DNA 2ng", "genomic DNA 0.7ng", "genomic DNA 0.2ng"]
        for tar in tar_roche:
            num = 0.0
            sum = 0
            for sam in samDna:
                if sam not in quant["tec_data"]:
                    continue
                if tar not in quant["tec_data"][sam]:
                    continue
                num += quant["tec_data"][sam][tar]["Ncopy_cv"]
                sum += 1
            keyVal = "dilution_" + tar.replace(" ", "_") + "_cv"
            if sum > 0:
                curDa[keyVal] = num / float(sum)
        for tar in tar_sensi:
            num = 0.0
            sum = 0
            for sam in samDna:
                if sam not in quant["tec_data"]:
                    continue
                if tar not in quant["tec_data"][sam]:
                    continue
                num += quant["tec_data"][sam][tar]["Ncopy_cv"]
                sum += 1
            keyVal = "dilution_" + tar.replace(" ", "_") + "_cv"
            if sum > 0:
                curDa[keyVal] = num / float(sum)
        for tar in tar_lc:
            num = 0.0
            sum = 0
            for sam in samDna:
                if sam not in quant["tec_data"]:
                    continue
                if tar not in quant["tec_data"][sam]:
                    continue
                num += quant["tec_data"][sam][tar]["Ncopy_cv"]
                sum += 1
            keyVal = "dilution_" + tar.replace(" ", "_") + "_cv"
            if sum > 0:
                curDa[keyVal] = num / float(sum)

        print("\n##########################\n### Test DNA Dilutions ###\n##########################")
        print("This test uses the DNA dilution data. The mean of the CV on all Ncopy of one target ")
        print("is calculated. It should be as low as possible.")
        print("----------------------")

        print("Roche:")
        printCVTable(tar_roche, "", curDa, laDa)
        print("\nSensi:")
        printCVTable(tar_roche, "_SF", curDa, laDa)
        print("\nLC-Green:")
        printCVTable(tar_roche, "_LC", curDa, laDa)

        print("\n#########################\n### Test DNA Standard ###\n#########################")
        print("This test uses the DNA dilution data. The parameters of a dilution standard are calculated. ")
        print("The efficiencies of the curve and the dilution should be similar. The Ncopy match the expectations.")
        print("----------------------")

        print("\n\nRoche:")
        printTable(20, "", tar_roche, head_a, prec_a, print_a, curDa, laDa)
        print("\nSensi:")
        printTable(20, "", tar_sensi, head_a, prec_a, print_a, curDa, laDa)
        print("\nLC-Green:")
        printTable(20, "", tar_lc, head_a, prec_a, print_a, curDa, laDa)
        print("\n\nRoche:")
        printTable(20, "", tar_roche, head_b, prec_b, print_b, curDa, laDa)
        print("\nSensi:")
        printTable(20, "", tar_sensi, head_b, prec_b, print_b, curDa, laDa)
        print("\nLC-Green:")
        printTable(20, "", tar_lc, head_b, prec_b, print_b, curDa, laDa)

curDa["AMP_reactionDataFalse"] = reactionDataFalse
curDa["AMP_reactionDataSum"] = reactionDataFalse + reactionDataTrue

printNice("Failing Reactions: ", curDa["AMP_reactionDataFalse"], laDa["AMP_reactionDataFalse"], 1, "-" , 0)
printNice("Sum Reactions: ", curDa["AMP_reactionDataSum"], laDa["AMP_reactionDataSum"], 1, "+", 0)

ww.close()
rd.save("temp_test_amplicon_primer_linregpcr.rdml")
endTime = time.time()
runTime = endTime - startTime
runMin = math.floor(runTime / 60.0)
runSec = runTime - runMin* 60.0
print("Runtime: " + str(runMin) + ":" + "{:.3f}".format(runSec))

# exit(0)

print("\n############################\n### Test Large Dilutions ###\n############################")
print("This test uses the large dilution data. It calculates the TD0 and the Ncopy. ")
print("The calculated values should match the expected values.")
print("----------------------")
startTime = time.time()
rd = rdml.Rdml(rdml_dil_file)

expList = rd.experiments()
if len(expList) < 1:
    print("No experiments found in " + rdml_dil_file + "!")
    sys.exit(0)
ww = open(out_dil_file, "w")
startLine = 0
reactionDataTrue = 0
reactionDataFalse = 0
linRegRes = {}
for exp in expList:
    runList = exp.runs()
    if len(runList) < 1:
        print("No runs found!")
        sys.exit(0)
    for run in runList:
        if printExpRun:
            print("Experiment: " + exp["id"] + " Run: " + run["id"])
        res = run.webAppLinRegPCR(pcrEfficiencyExl=0.05, updateTargetEfficiency=False, updateRDML=True, excludeNoPlateau=True, excludeEfficiency="outlier", excludeInstableBaseline=True)
        resTab = json.loads(res["LinRegPCR_Result_Table"])
        for tabRow in range(startLine, len(resTab)):
            if startLine == 0:
                ww.write("Experiment\tRun\t")
                startLine = 1
            else:
                ww.write(exp["id"] + "\t" + run["id"] + "\t")
            for tabCol in range(0, len(resTab[tabRow])):
                outCell = str(resTab[tabRow][tabCol]).replace("\t", ";")
                if tabCol < len(resTab[tabRow]) - 1:
                    ww.write(outCell + "\t")
                else:
                    ww.write(outCell + "\n")
            if startLine == 1:
                if resTab[tabRow][rar_sample_type] in ["unkn", "std"]:
                    if float(resTab[tabRow][rar_Ncopy]) > 5.0:
                        if exp["id"] not in linRegRes:
                            linRegRes[exp["id"]] = {}
                        if run["id"] not in linRegRes[exp["id"]]:
                            linRegRes[exp["id"]][run["id"]] = {}
                        if resTab[tabRow][rar_tar] not in linRegRes[exp["id"]][run["id"]]:
                            linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]] = {}
                        if resTab[tabRow][rar_sample] not in linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]:
                            linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]][resTab[tabRow][rar_sample]] = {}
                            linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]][resTab[tabRow][rar_sample]]["Ncopy"] = []
                            linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]][resTab[tabRow][rar_sample]]["TD0"] = []
                        linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]][resTab[tabRow][rar_sample]]["Ncopy"].append(resTab[tabRow][rar_Ncopy])
                        linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]][resTab[tabRow][rar_sample]]["TD0"].append(resTab[tabRow][rar_TD0])
                        linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["PCReff"] = resTab[tabRow][rar_PCR_eff]

                        reactionDataTrue += 1
                    else:
                        reactionDataFalse += 1

    if exp["id"] == "DNA Dilutions - EMBL":
        dil_test_key = "Large_DNA_Dil_A_"
        dil_test_tar = ["FSTL_1_F_109", "FSTL_1_K_219", "FSTL_1_H_201", "FSTL_1_L_412"]
        dil_test_sam = ["Dil BSG - 10e8", "Dil BSG - 10e7", "Dil BSG - 10e6", "Dil BSG - 10e5", "Dil BSG - 10e4", "Dil BSG - 10e3", "Dil BSG - 10e2",
                        "Dil 1:10 - 10e8", "Dil 1:10 - 10e7", "Dil 1:10 - 10e6", "Dil 1:10 - 10e5", "Dil 1:10 - 10e4", "Dil 1:10 - 10e3", "Dil 1:10 - 10e2"]
        dil_test_sam_A = ["Dil 1:10 - 10e8", "Dil 1:10 - 10e7", "Dil 1:10 - 10e6", "Dil 1:10 - 10e5", "Dil 1:10 - 10e4", "Dil 1:10 - 10e3", "Dil 1:10 - 10e2"]
        dil_test_sam_B = ["Dil BSG - 10e8", "Dil BSG - 10e7", "Dil BSG - 10e6", "Dil BSG - 10e5", "Dil BSG - 10e4", "Dil BSG - 10e3", "Dil BSG - 10e2"]
        dil_test_sam_C = ["10e8", "10e7", "10e6", "10e5", "10e4", "10e3", "10e2"]
        for tar in dil_test_tar:
            for sam in dil_test_sam:
                if sam not in linRegRes["DNA Dilutions - EMBL"]["P10 - DNA Dilutions"][resTab[tabRow][rar_tar]]:
                    continue
                curDa[dil_test_key + "PCReff_" + tar] = linRegRes["DNA Dilutions - EMBL"]["P10 - DNA Dilutions"][tar]["PCReff"]
                curDa[dil_test_key + "TD0_" + tar + "_" + sam] = np.mean(linRegRes["DNA Dilutions - EMBL"]["P10 - DNA Dilutions"][tar][sam]["TD0"])
                curDa[dil_test_key + "Ncopy_" + tar + "_" + sam] = np.mean(linRegRes["DNA Dilutions - EMBL"]["P10 - DNA Dilutions"][tar][sam]["Ncopy"])

        # PCR eff
        out_str = " ".ljust(17)
        for tar in dil_test_tar:
            out_str +=  tar.ljust(17)
        print(out_str)
        out_str = "PCR eff".ljust(17)
        for tar in dil_test_tar:
            out_str +=  colorDiff(curDa, laDa, dil_test_key + "PCReff_" + tar, "{:6.3f}")  + "  "
        print(out_str)
        print("\n")
        # TD0
        out_str = " ".ljust(10)
        for tar in dil_test_tar:
            out_str +=  tar.ljust(17)
            out_str +=  tar.ljust(17)
        print(out_str)
        out_str = " ".ljust(10)
        for pos in range(0,8):
            if pos % 2 == 0:
                out_str +=  "Dil 1:10".ljust(17)
            else:
                out_str +=  "Dil BSG".ljust(17)
        print(out_str)
        for pos in range(0,7): 
            out_str = dil_test_sam_C[pos].ljust(10)
            for tar in dil_test_tar:
                out_str +=  colorDiff(curDa, laDa, dil_test_key + "TD0_" + tar + "_" + dil_test_sam_A[pos], "{:6.3f}") + "  "
                out_str +=  colorDiff(curDa, laDa, dil_test_key + "TD0_" + tar + "_" + dil_test_sam_B[pos], "{:6.3f}") + "  "
            print(out_str)
        print("\n")
        # Ncopy
        out_str = " ".ljust(10)
        for tar in dil_test_tar:
            out_str +=  tar.ljust(23)
            out_str +=  tar.ljust(23)
        print(out_str)
        out_str = " ".ljust(10)
        for pos in range(0,8):
            if pos % 2 == 0:
                out_str +=  "Dil 1:10".ljust(23)
            else:
                out_str +=  "Dil BSG".ljust(23)
        print(out_str)
        for pos in range(0,7): 
            out_str = dil_test_sam_C[pos].ljust(10)
            for tar in dil_test_tar:
                out_str +=  colorDiff(curDa, laDa, dil_test_key + "Ncopy_" + tar + "_" + dil_test_sam_A[pos], "{:9.0f}") + "  "
                out_str +=  colorDiff(curDa, laDa, dil_test_key + "Ncopy_" + tar + "_" + dil_test_sam_B[pos], "{:9.0f}") + "  "
            print(out_str)
        print("\n")

    if exp["id"] == "dPCR - DKFZ":
        dil_test_key = "Large_DNA_Dil_B_"
        dil_test_tar = ["FSTL_1_F_109", "FSTL_1_K_219", "FSTL_1_H_201", "FSTL_1_L_412"]
        dil_test_sam = ["Freeze 10e8", "Freeze 10e7", "Freeze 10e6", "Freeze 10e5", "Freeze 10e4", "Freeze 10e3", "Freeze 10e2",
                        "Dil 1:10 - 10e8", "Dil 1:10 - 10e7", "Dil 1:10 - 10e6", "Dil 1:10 - 10e5", "Dil 1:10 - 10e4", "Dil 1:10 - 10e3", "Dil 1:10 - 10e2"]
        dil_test_sam_A = ["Dil 1:10 - 10e8", "Dil 1:10 - 10e7", "Dil 1:10 - 10e6", "Dil 1:10 - 10e5", "Dil 1:10 - 10e4", "Dil 1:10 - 10e3", "Dil 1:10 - 10e2"]
        dil_test_sam_B = ["Freeze 10e8", "Freeze 10e7", "Freeze 10e6", "Freeze 10e5", "Freeze 10e4", "Freeze 10e3", "Freeze 10e2"]
        dil_test_sam_C = ["10e8", "10e7", "10e6", "10e5", "10e4", "10e3", "10e2"]
        for tar in dil_test_tar:
            for sam in dil_test_sam:
                if sam not in linRegRes["dPCR - DKFZ"]["p11 - dPCR"][resTab[tabRow][rar_tar]]:
                    continue
                curDa[dil_test_key + "PCReff_" + tar] = linRegRes["dPCR - DKFZ"]["p11 - dPCR"][tar]["PCReff"]
                curDa[dil_test_key + "TD0_" + tar + "_" + sam] = np.mean(linRegRes["dPCR - DKFZ"]["p11 - dPCR"][tar][sam]["TD0"])
                curDa[dil_test_key + "Ncopy_" + tar + "_" + sam] = np.mean(linRegRes["dPCR - DKFZ"]["p11 - dPCR"][tar][sam]["Ncopy"])
        # PCR eff
        out_str = " ".ljust(17)
        for tar in dil_test_tar:
            out_str +=  tar.ljust(17)
        print(out_str)
        out_str = "PCR eff".ljust(17)
        for tar in dil_test_tar:
            out_str +=  colorDiff(curDa, laDa, dil_test_key + "PCReff_" + tar, "{:6.3f}")  + "  "
        print(out_str)
        print("\n")
        # TD0
        out_str = " ".ljust(10)
        for tar in dil_test_tar:
            out_str +=  tar.ljust(17)
            out_str +=  tar.ljust(17)
        print(out_str)
        out_str = " ".ljust(10)
        for pos in range(0,8):
            if pos % 2 == 0:
                out_str +=  "Dil 1:10".ljust(17)
            else:
                out_str +=  "Freeze".ljust(17)
        print(out_str)
        for pos in range(0,7): 
            out_str = dil_test_sam_C[pos].ljust(10)
            for tar in dil_test_tar:
                out_str +=  colorDiff(curDa, laDa, dil_test_key + "TD0_" + tar + "_" + dil_test_sam_A[pos], "{:6.3f}") + "  "
                out_str +=  colorDiff(curDa, laDa, dil_test_key + "TD0_" + tar + "_" + dil_test_sam_B[pos], "{:6.3f}") + "  "
            print(out_str)
        print("\n")
        # Ncopy
        out_str = " ".ljust(10)
        for tar in dil_test_tar:
            out_str +=  tar.ljust(23)
            out_str +=  tar.ljust(23)
        print(out_str)
        out_str = " ".ljust(10)
        for pos in range(0,8):
            if pos % 2 == 0:
                out_str +=  "Dil 1:10".ljust(23)
            else:
                out_str +=  "Freeze".ljust(23)
        print(out_str)
        for pos in range(0,7): 
            out_str = dil_test_sam_C[pos].ljust(10)
            for tar in dil_test_tar:
                out_str +=  colorDiff(curDa, laDa, dil_test_key + "Ncopy_" + tar + "_" + dil_test_sam_A[pos], "{:9.0f}") + "  "
                out_str +=  colorDiff(curDa, laDa, dil_test_key + "Ncopy_" + tar + "_" + dil_test_sam_B[pos], "{:9.0f}") + "  "
            print(out_str)
        print("\n")


ww.close()
rd.save("temp_test_large_DNA_dilutions_linregpcr.rdml")
endTime = time.time()
runTime = endTime - startTime
runMin = math.floor(runTime / 60.0)
runSec = runTime - runMin* 60.0
print("Runtime: " + str(runMin) + ":" + "{:.3f}".format(runSec))

print("\n###################\n### Test Probes ###\n###################")
print("This test uses DNA dilutions but quantifies by probes. First calculate results.")
print("----------------------")

# Time the test
startTime = time.time()
rd = rdml.Rdml(rdml_probes_file)

expList = rd.experiments()
if len(expList) < 1:
    print("No experiments found!")
    sys.exit(0)
ww = open(out_probes_file, "w")
startLine = 0
linRegRes = {}
reactionDataTrue = 0
reactionDataFalse = 0
colTar = []
quant = []
for exp in expList:
    runList = exp.runs()
    if len(runList) < 1:
        print("No runs found!")
        sys.exit(0)
    for run in runList:
        if printExpRun:
            print("Experiment: " + exp["id"] + " Run: " + run["id"])
        res = run.webAppLinRegPCR(pcrEfficiencyExl=0.05, updateTargetEfficiency=True, updateRDML=True, excludeNoPlateau=True, excludeEfficiency="outlier", excludeInstableBaseline=True)
        resTab = json.loads(res["LinRegPCR_Result_Table"])
        for tabRow in range(0, len(resTab)):
            if startLine == 0:
                ww.write("Experiment\tRun\t")
                startLine = 1
            else:
                ww.write(exp["id"] + "\t" + run["id"] + "\t")
            for tabCol in range(0, len(resTab[tabRow])):
                outCell = str(resTab[tabRow][tabCol]).replace("\t", ";")
                if tabCol < len(resTab[tabRow]) - 1:
                    ww.write(outCell + "\t")
                else:
                    ww.write(outCell + "\n")
            if startLine == 1:
                if resTab[tabRow][rar_sample_type] in ["unkn", "std"]:
                    if float(resTab[tabRow][rar_Ncopy]) > 0.0:
                        if exp["id"] not in linRegRes:
                            linRegRes[exp["id"]] = {}
                        if run["id"] not in linRegRes[exp["id"]]:
                            linRegRes[exp["id"]][run["id"]] = {}
                        if resTab[tabRow][rar_tar] not in linRegRes[exp["id"]][run["id"]]:
                            linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]] = {}
                            linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["Ncopy"] = []
                            linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["TD0"] = []

                        linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["Ncopy"].append(resTab[tabRow][rar_Ncopy])
                        linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["TD0"].append(resTab[tabRow][rar_TD0])
                        linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["PCReff"] = resTab[tabRow][rar_PCR_eff]

                        reactionDataTrue += 1
                    else:
                        reactionDataFalse += 1

    if exp["id"] == "Probes - AMC":
        quant = exp.quantify()
        head_var = ["curve_pcr_eff", "dilution_pcr_eff", "expected_max", "mean_max", "expected_min", "mean_min", 
                    "expected_ratio", "mean_ratio", "slope_bias", "correlation_R", "linearity", 
                    "reproducibility", "detectable_diff"]
        head_a = ["curve_pcr_eff", "dilution_pcr_eff", "expected_max", "mean_max", "expected_min", "mean_min", 
                    "expected_ratio", "mean_ratio"]
        head_b = ["slope_bias", "correlation_R", "linearity", "reproducibility", "detectable_diff"]
        prec_a = [0.001, 0.001, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        prec_b = [0.0001, 0.00001, 0.00001, 0.00001, 0.0001]
        print_a = [3, 3, 0, 0, 0, 0, 0, 0]
        print_b = [5, 5, 5, 5, 5]
        tar_roche = ["FSTL_1_D_105 250nM 150nM", "FSTL_1_F_109 250nM 150nM", "FSTL_1_K_219 250nM 150nM", "FSTL_1_O_417 250nM 150nM"]
        tar_idt = ["FSTL_1_D_105 250nM 150nM IDT", "FSTL_1_F_109 250nM 150nM IDT", "FSTL_1_K_219 250nM 150nM IDT", "FSTL_1_O_417 250nM 150nM IDT"]

        for tar in tar_roche:
            for var in head_var:
                keyVal = tar.replace(" ", "_") + "_" + var
                curDa[keyVal] = quant["dil_std"][tar][var]
        for tar in tar_idt:
            for var in head_var:
                keyVal = tar.replace(" ", "_") + "_" + var
                curDa[keyVal] = quant["dil_std"][tar][var]

        print("\n#####################################\n### Test Probes with DNA Standard ###\n#####################################")
        print("This test uses the DNA dilution but quantifies by probes. The parameters of a dilution ")
        print("standard are calculated. The efficiencies of the curve and the dilution should be similar. ")
        print("The Ncopy match the expectations.")
        print("----------------------")

        print("\n\nRoche:")
        printTable(30, "", tar_roche, head_a, prec_a, print_a, curDa, laDa)
        print("\nSensi:")
        printTable(30, "", tar_idt, head_a, prec_a, print_a, curDa, laDa)
        print("\n\nRoche:")
        printTable(30, "", tar_roche, head_b, prec_b, print_b, curDa, laDa)
        print("\nSensi:")
        printTable(30, "", tar_idt, head_b, prec_b, print_b, curDa, laDa)

        print("\n#############################################\n### Test Probes Primer Conc vs Probe Conc ###\n#############################################")
        print("This test uses the 4ng DNA dilution data but quantifies by different primer and probe concentraion.")
        print("The TD0 should be good, PCR efficiencies are calculated of 3 replicates and have big noise. Ncopy ")
        print("is noisy due to the PCR efficiencies.")
        print("----------------------")

        probe_tar = ["FSTL_1_D_105 ", "FSTL_1_F_109 ", "FSTL_1_K_219 ", "FSTL_1_O_417 "]
        probe_conc = ["500nM", "250nM", "150nM", "75nM", "50nM", "25nM"]
        primer_conc = ["100nM", "250nM", "750nM"]
        res = ""
        for prim in primer_conc:
            res += (prim + " TD0").ljust(15)
            for col in probe_conc:
                res += col.ljust(16)
            res += "\n"
            for tar in probe_tar:
                res += tar.ljust(15)
                count = 0
                for probeC in probe_conc:
                    tar_id = tar + prim + " " + probeC
                    keyVal = "probe_prim_prob_" + tar_id.replace(" ", "_") + "_TD0"
                    curDa[keyVal] = -1.0
                    if tar_id in linRegRes[exp["id"]]["P8 - Primer conc - Probe conc"]:
                        curDa[keyVal] = np.mean(linRegRes[exp["id"]]["P8 - Primer conc - Probe conc"][tar_id]["TD0"])
                    res += colorDiff(curDa, laDa, keyVal, "{:6.2f}") + " "
                res += "\n"
            res += "\n"
        print(res)
        res = ""
        for prim in primer_conc:
            res += (prim + " PCReff").ljust(15)
            for col in probe_conc:
                res += col.ljust(18)
            res += "\n"
            for tar in probe_tar:
                res += tar.ljust(15)
                count = 0
                for probeC in probe_conc:
                    tar_id = tar + prim + " " + probeC
                    keyVal = "probe_prim_prob_" + tar_id.replace(" ", "_") + "_PCReff"
                    curDa[keyVal] = -1.0
                    if tar_id in linRegRes[exp["id"]]["P8 - Primer conc - Probe conc"]:
                        curDa[keyVal] = linRegRes[exp["id"]]["P8 - Primer conc - Probe conc"][tar_id]["PCReff"]
                    res += colorDiff(curDa, laDa, keyVal, "{:7.4f}") + " "
                res += "\n"
            res += "\n"
        print(res)
        res = ""
        for prim in primer_conc:
            res += (prim + " Ncopy").ljust(17)
            for col in probe_conc:
                res += col.ljust(20)
            res += "\n"
            for tar in probe_tar:
                res += tar.ljust(15)
                count = 0
                for probeC in probe_conc:
                    tar_id = tar + prim + " " + probeC
                    keyVal = "probe_prim_prob_" + tar_id.replace(" ", "_") + "_Ncopy"
                    curDa[keyVal] = -1.0
                    if tar_id in linRegRes[exp["id"]]["P8 - Primer conc - Probe conc"]:
                        curDa[keyVal] = np.mean(linRegRes[exp["id"]]["P8 - Primer conc - Probe conc"][tar_id]["Ncopy"])
                    res += colorDiff(curDa, laDa, keyVal, "{:8.1f}") + " "
                res += "\n"
            res += "\n"
        print(res)

ww.close()
rd.save("temp_probes.rdml")


print("\n#######################################\n### Test Own SYBR Mix - DMSO and Mg ###\n#######################################")
print("This test uses the FastStart Taq DNA Polymerase from Roche (Cat. No. 04 738 357 001) and adds defined amounts of SYBR I. ")
print("5ng gDNA were used 1500 copies should be expected. As PCR efficiency is based on 3 reactions, Ncopy will have big noise. ")
print("----------------------")

rd = rdml.Rdml(rdml_ownmix_a_file)

expList = rd.experiments()
if len(expList) < 1:
    print("No experiments found!")
    sys.exit(0)
ww = open(out_ownmix_a_file, "w")
startLine = 0
linRegRes = {}
colTar = []
quant = []
for exp in expList:
    runList = exp.runs()
    if len(runList) < 1:
        print("No runs found!")
        sys.exit(0)
    for run in runList:
        if printExpRun:
            print("Experiment: " + exp["id"] + " Run: " + run["id"])
        res = run.webAppLinRegPCR(pcrEfficiencyExl=0.05, updateTargetEfficiency=True, updateRDML=True, excludeNoPlateau=True, excludeEfficiency="outlier", excludeInstableBaseline=True)
        resTab = json.loads(res["LinRegPCR_Result_Table"])
        for tabRow in range(0, len(resTab)):
            if startLine == 0:
                ww.write("Experiment\tRun\t")
                startLine = 1
            else:
                ww.write(exp["id"] + "\t" + run["id"] + "\t")
            for tabCol in range(0, len(resTab[tabRow])):
                outCell = str(resTab[tabRow][tabCol]).replace("\t", ";")
                if tabCol < len(resTab[tabRow]) - 1:
                    ww.write(outCell + "\t")
                else:
                    ww.write(outCell + "\n")
            if startLine == 1:
                if resTab[tabRow][rar_sample_type] in ["unkn", "std"]:
                    if float(resTab[tabRow][rar_Ncopy]) > 5.0:
                        if exp["id"] not in linRegRes:
                            linRegRes[exp["id"]] = {}
                        if run["id"] not in linRegRes[exp["id"]]:
                            linRegRes[exp["id"]][run["id"]] = {}
                        if resTab[tabRow][rar_tar] not in colTar:
                            colTar.append(resTab[tabRow][rar_tar])
                        if resTab[tabRow][rar_tar] not in linRegRes[exp["id"]][run["id"]]:
                            linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]] = {}
                            linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["Ncopy"] = []
                            linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["TD0"] = []
                            linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["PCReff"] = resTab[tabRow][rar_PCR_eff]
                        
                        linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["Ncopy"].append(resTab[tabRow][rar_Ncopy])  # Ncopy
                        linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["TD0"].append(resTab[tabRow][rar_TD0])  # TD0

        for tar in colTar:
            if tar in linRegRes[exp["id"]][run["id"]]:
                linRegRes[exp["id"]][run["id"]][tar]["Ncopy mean"] = np.mean(linRegRes[exp["id"]][run["id"]][tar]["Ncopy"])
                linRegRes[exp["id"]][run["id"]][tar]["TD0 mean"] = np.mean(linRegRes[exp["id"]][run["id"]][tar]["TD0"])    

for dmso in ["0.0", "2.5", "5.0"]:
    for mg in ["3.0", "5.0"]:
        curDa["Test_OWNMIX_A_" + dmso + "_" + mg + "_TD0"] = linRegRes["Own Mix A"]["Run 1"]["DMSO " + dmso + " - Mg " + mg]["TD0 mean"]
        curDa["Test_OWNMIX_A_" + dmso + "_" + mg + "_PCReff"] = linRegRes["Own Mix A"]["Run 1"]["DMSO " + dmso + " - Mg " + mg]["PCReff"]
        curDa["Test_OWNMIX_A_" + dmso + "_" + mg + "_Ncopy"] = linRegRes["Own Mix A"]["Run 1"]["DMSO " + dmso + " - Mg " + mg]["Ncopy mean"]

res = "\n                    TD0               PCReff              Ncopy\n"
for dmso in ["0.0", "2.5", "5.0"]:
    for mg in ["3.0", "5.0"]:
        res += ("DMSO " + dmso + " - Mg " + mg).ljust(20)
        res += colorDiff(curDa, laDa, "Test_OWNMIX_A_" + dmso + "_" + mg + "_TD0", "{:6.2f}") + "   "
        res += colorDiff(curDa, laDa, "Test_OWNMIX_A_" + dmso + "_" + mg + "_PCReff", "{:7.4f}") + "   "
        res += colorDiff(curDa, laDa, "Test_OWNMIX_A_" + dmso + "_" + mg + "_Ncopy", "{:7.1f}")
        res += "\n"
print(res)

ww.close()
rd.save("temp_ownmix_a.rdml")


print("\n########################################\n### Test Own SYBR Mix - SYBR and BSA ###\n########################################")
print("This test uses the FastStart Taq DNA Polymerase from Roche (Cat. No. 04 738 357 001) and adds defined amounts of SYBR I. ")
print("5ng gDNA were used 1500 copies should be expected. As PCR efficiency is based on 3 reactions, Ncopy will have big noise. ")
print("----------------------")

rd = rdml.Rdml(rdml_ownmix_b_file)

expList = rd.experiments()
if len(expList) < 1:
    print("No experiments found!")
    sys.exit(0)
ww = open(out_ownmix_b_file, "w")
startLine = 0
linRegRes = {}
colTar = []
quant = []
for exp in expList:
    runList = exp.runs()
    if len(runList) < 1:
        print("No runs found!")
        sys.exit(0)
    for run in runList:
        if printExpRun:
            print("Experiment: " + exp["id"] + " Run: " + run["id"])
        res = run.webAppLinRegPCR(pcrEfficiencyExl=0.05, updateTargetEfficiency=True, updateRDML=True, excludeNoPlateau=True, excludeEfficiency="outlier", excludeInstableBaseline=True)
        resTab = json.loads(res["LinRegPCR_Result_Table"])
        for tabRow in range(0, len(resTab)):
            if startLine == 0:
                ww.write("Experiment\tRun\t")
                startLine = 1
            else:
                ww.write(exp["id"] + "\t" + run["id"] + "\t")
            for tabCol in range(0, len(resTab[tabRow])):
                outCell = str(resTab[tabRow][tabCol]).replace("\t", ";")
                if tabCol < len(resTab[tabRow]) - 1:
                    ww.write(outCell + "\t")
                else:
                    ww.write(outCell + "\n")
            if startLine == 1:
                if resTab[tabRow][rar_sample_type] in ["unkn", "std"]:
                    if exp["id"] not in linRegRes:
                        linRegRes[exp["id"]] = {}
                    if run["id"] not in linRegRes[exp["id"]]:
                        linRegRes[exp["id"]][run["id"]] = {}
                    if resTab[tabRow][rar_tar] not in colTar:
                        colTar.append(resTab[tabRow][rar_tar])
                    if resTab[tabRow][rar_tar] not in linRegRes[exp["id"]][run["id"]]:
                        linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]] = {}
                        linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["Ncopy"] = []
                        linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["TD0"] = []
                        linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["PCReff"] = resTab[tabRow][rar_PCR_eff]
                    
                    linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["Ncopy"].append(resTab[tabRow][rar_Ncopy])  # Ncopy
                    linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["TD0"].append(resTab[tabRow][rar_TD0])  # TD0

        for tar in colTar:
            if tar in linRegRes[exp["id"]][run["id"]]:
                linRegRes[exp["id"]][run["id"]][tar]["Ncopy mean"] = np.mean(linRegRes[exp["id"]][run["id"]][tar]["Ncopy"])
                linRegRes[exp["id"]][run["id"]][tar]["TD0 mean"] = np.mean(linRegRes[exp["id"]][run["id"]][tar]["TD0"])    

for dmso in ["SYBR", "No_BSA", "Sigma"]:
    for conc in ["100", "200", "400", "800"]:
        useConc = " BSA"
        startStr = dmso
        if dmso == "No_BSA":
            useConc = ""
            startStr = "SYBR"
        curDa["Test_OWNMIX_B_" + dmso + "_" + conc + "_TD0"] = linRegRes["Own Mix B"]["Run 1"][startStr + " 1:" + conc + useConc]["TD0 mean"]
        curDa["Test_OWNMIX_B_" + dmso + "_" + conc + "_PCReff"] = linRegRes["Own Mix B"]["Run 1"][startStr + " 1:" + conc + useConc]["PCReff"]
        curDa["Test_OWNMIX_B_" + dmso + "_" + conc + "_Ncopy"] = linRegRes["Own Mix B"]["Run 1"][startStr + " 1:" + conc + useConc]["Ncopy mean"]

res = "\n                    TD0               PCReff              Ncopy\n"
for conc in ["100", "200", "400", "800"]:
    for dmso in ["SYBR", "No_BSA", "Sigma"]:
        useConc = " BSA"
        startStr = dmso
        if dmso == "No_BSA":
            useConc = ""
            startStr = "SYBR"
        res += (startStr.ljust(5) + " 1:" + conc + useConc).ljust(20)
        res += colorDiff(curDa, laDa, "Test_OWNMIX_B_" + dmso + "_" + conc + "_TD0", "{:6.2f}") + "   "
        res += colorDiff(curDa, laDa, "Test_OWNMIX_B_" + dmso + "_" + conc + "_PCReff", "{:7.4f}") + "   "
        res += colorDiff(curDa, laDa, "Test_OWNMIX_B_" + dmso + "_" + conc + "_Ncopy", "{:7.1f}")
        res += "\n"
print(res)

ww.close()
rd.save("temp_ownmix_b.rdml")


print("\n#######################################\n### Test Own SYBR Mix - Primer Conc ###\n#######################################")
print("This test uses the FastStart Taq DNA Polymerase from Roche (Cat. No. 04 738 357 001) and adds defined amounts of SYBR I. ")
print("5ng gDNA were used 1500 copies should be expected. This experiment had many failing technical replicates, Ncopy will ")
print("have big noise. ")
print("----------------------")

rd = rdml.Rdml(rdml_ownmix_c_file)

expList = rd.experiments()
if len(expList) < 1:
    print("No experiments found!")
    sys.exit(0)
ww = open(out_ownmix_c_file, "w")
startLine = 0
linRegRes = {}
colTar = []
quant = []
for exp in expList:
    runList = exp.runs()
    if len(runList) < 1:
        print("No runs found!")
        sys.exit(0)
    for run in runList:
        if printExpRun:
            print("Experiment: " + exp["id"] + " Run: " + run["id"])
        res = run.webAppLinRegPCR(pcrEfficiencyExl=0.05, updateTargetEfficiency=True, updateRDML=True, excludeNoPlateau=False, excludeEfficiency="outlier", excludeInstableBaseline=False)
        resTab = json.loads(res["LinRegPCR_Result_Table"])
        for tabRow in range(0, len(resTab)):
            if startLine == 0:
                ww.write("Experiment\tRun\t")
                startLine = 1
            else:
                ww.write(exp["id"] + "\t" + run["id"] + "\t")
            for tabCol in range(0, len(resTab[tabRow])):
                outCell = str(resTab[tabRow][tabCol]).replace("\t", ";")
                if tabCol < len(resTab[tabRow]) - 1:
                    ww.write(outCell + "\t")
                else:
                    ww.write(outCell + "\n")
            if startLine == 1:
                if resTab[tabRow][rar_sample_type] in ["unkn", "std"]:
                    if float(resTab[tabRow][rar_Ncopy]) > 5.0:
                        if exp["id"] not in linRegRes:
                            linRegRes[exp["id"]] = {}
                        if run["id"] not in linRegRes[exp["id"]]:
                            linRegRes[exp["id"]][run["id"]] = {}
                        if resTab[tabRow][rar_tar] not in colTar:
                            colTar.append(resTab[tabRow][rar_tar])
                        if resTab[tabRow][rar_tar] not in linRegRes[exp["id"]][run["id"]]:
                            linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]] = {}
                            linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["Ncopy"] = []
                            linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["TD0"] = []
                            linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["PCReff"] = resTab[tabRow][rar_PCR_eff]

                        linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["Ncopy"].append(resTab[tabRow][rar_Ncopy])  # Ncopy
                        linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["TD0"].append(resTab[tabRow][rar_TD0])  # TD0

        for tar in colTar:
            if tar in linRegRes[exp["id"]][run["id"]]:
                linRegRes[exp["id"]][run["id"]][tar]["Ncopy mean"] = np.mean(linRegRes[exp["id"]][run["id"]][tar]["Ncopy"])
                linRegRes[exp["id"]][run["id"]][tar]["TD0 mean"] = np.mean(linRegRes[exp["id"]][run["id"]][tar]["TD0"])    
    
tars = ["FSTL_1_A_047", "FSTL_1_B_042", "FSTL_1_D_105", "FSTL_1_E_097", "FSTL_1_F_109",
        "FSTL_1_H_201", "FSTL_1_I_204", "FSTL_1_K_219", "FSTL_1_*_259", "FSTL_1_L_412"]

for tar in tars:
    for conc in ["100nM", "250nM", "750nM"]:
        useConc = " " + conc
        if conc == "250nM":
            useConc = ""
        curDa["Test_OWNMIX_C_" + tar + "_" + conc + "_TD0"] = linRegRes["Own Mix C"]["Run 1"][tar + useConc]["TD0 mean"]
        curDa["Test_OWNMIX_C_" + tar + "_" + conc + "_PCReff"] = linRegRes["Own Mix C"]["Run 1"][tar + useConc]["PCReff"]
        curDa["Test_OWNMIX_C_" + tar + "_" + conc + "_Ncopy"] = linRegRes["Own Mix C"]["Run 1"][tar + useConc]["Ncopy mean"]

res = "\n                    TD0               PCReff              Ncopy\n"
for conc in ["100nM", "250nM", "750nM"]:
    for tar in tars:
        res += (tar + " " + conc).ljust(20)
        res += colorDiff(curDa, laDa, "Test_OWNMIX_C_" + tar + "_" + conc + "_TD0", "{:6.2f}") + "   "
        res += colorDiff(curDa, laDa, "Test_OWNMIX_C_" + tar + "_" + conc + "_PCReff", "{:7.4f}") + "   "
        res += colorDiff(curDa, laDa, "Test_OWNMIX_C_" + tar + "_" + conc + "_Ncopy", "{:7.1f}")
        res += "\n"
    res += "\n"
print(res)

ww.close()
rd.save("temp_ownmix_c.rdml")


print("\n###################################\n### Test Added Dye to Probe Mix ###\n###################################")
print("This test uses the probe mixes from Roche and IDT and add defined amounts of SYBR I. ")
print("7ng DNA were used 2600 copies should be expected. As PCR efficiency is based on ")
print("3 reactions, Ncopy will have big noise. ")
print("----------------------")

rd = rdml.Rdml(rdml_add_sybr_file)

expList = rd.experiments()
if len(expList) < 1:
    print("No experiments found!")
    sys.exit(0)
ww = open(out_add_sybr_file, "w")
startLine = 0
linRegRes = {}
colTar = []
quant = []
for exp in expList:
    runList = exp.runs()
    if len(runList) < 1:
        print("No runs found!")
        sys.exit(0)
    for run in runList:
        if printExpRun:
            print("Experiment: " + exp["id"] + " Run: " + run["id"])
        res = run.webAppLinRegPCR(pcrEfficiencyExl=0.05, updateTargetEfficiency=True, updateRDML=True, excludeNoPlateau=True, excludeEfficiency="outlier", excludeInstableBaseline=True)
        resTab = json.loads(res["LinRegPCR_Result_Table"])
        for tabRow in range(0, len(resTab)):
            if startLine == 0:
                ww.write("Experiment\tRun\t")
                startLine = 1
            else:
                ww.write(exp["id"] + "\t" + run["id"] + "\t")
            for tabCol in range(0, len(resTab[tabRow])):
                outCell = str(resTab[tabRow][tabCol]).replace("\t", ";")
                if tabCol < len(resTab[tabRow]) - 1:
                    ww.write(outCell + "\t")
                else:
                    ww.write(outCell + "\n")
            if startLine == 1:
                if resTab[tabRow][rar_sample_type] in ["unkn", "std"]:
                    if exp["id"] not in linRegRes:
                        linRegRes[exp["id"]] = {}
                    if run["id"] not in linRegRes[exp["id"]]:
                        linRegRes[exp["id"]][run["id"]] = {}
                    if resTab[tabRow][rar_tar] not in colTar:
                        colTar.append(resTab[tabRow][rar_tar])
                    if resTab[tabRow][rar_tar] not in linRegRes[exp["id"]][run["id"]]:
                        linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]] = {}
                        linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["Ncopy"] = []
                        linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["TD0"] = []
                        linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["PCReff"] = resTab[tabRow][rar_PCR_eff]
                    
                    linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["Ncopy"].append(resTab[tabRow][rar_Ncopy])  # Ncopy
                    linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["TD0"].append(resTab[tabRow][rar_TD0])  # TD0

        for tar in colTar:
            if tar in linRegRes[exp["id"]][run["id"]]:
                linRegRes[exp["id"]][run["id"]][tar]["Ncopy mean"] = np.mean(linRegRes[exp["id"]][run["id"]][tar]["Ncopy"])
                linRegRes[exp["id"]][run["id"]][tar]["TD0 mean"] = np.mean(linRegRes[exp["id"]][run["id"]][tar]["TD0"])    

    
for conc in ["200", "400", "800", "1600"]:
    for uMix in ["ROCHE", "IDT"]:
        addMix = ""
        if uMix == "IDT":
            addMix = " IDT"
        curDa["Test_ADD_DYE_" + uMix + "_" + conc + "_TD0"] = linRegRes["Probe Mix SYBR"]["P4 - SYBR conc"]["FSTL_1_*_259 250nM 1:" + conc + addMix]["TD0 mean"]
        curDa["Test_ADD_DYE_" + uMix + "_" + conc + "_PCReff"] = linRegRes["Probe Mix SYBR"]["P4 - SYBR conc"]["FSTL_1_*_259 250nM 1:" + conc + addMix]["PCReff"]
        curDa["Test_ADD_DYE_" + uMix + "_" + conc + "_Ncopy"] = linRegRes["Probe Mix SYBR"]["P4 - SYBR conc"]["FSTL_1_*_259 250nM 1:" + conc + addMix]["Ncopy mean"]

res = "\n             TD0                                PCReff                                 Ncopy\n"
res +=  "             Roche            IDT               Roche              IDT                 Roche            IDT\n"
for conc in ["200", "400", "800", "1600"]:
    res += "SYBR 1:" + conc.ljust(6)
    res += colorDiff(curDa, laDa, "Test_ADD_DYE_ROCHE_" + conc + "_TD0", "{:6.2f}") + "  "
    res += colorDiff(curDa, laDa, "Test_ADD_DYE_IDT_" + conc + "_TD0", "{:6.2f}") + "   "
    res += colorDiff(curDa, laDa, "Test_ADD_DYE_ROCHE_" + conc + "_PCReff", "{:7.4f}") + "  "
    res += colorDiff(curDa, laDa, "Test_ADD_DYE_IDT_" + conc + "_PCReff", "{:7.4f}") + "   "
    res += colorDiff(curDa, laDa, "Test_ADD_DYE_ROCHE_" + conc + "_Ncopy", "{:6.1f}") + "  "
    res += colorDiff(curDa, laDa, "Test_ADD_DYE_IDT_" + conc + "_Ncopy", "{:6.1f}") + " "
    res += "\n"
print(res)

add_dye_test_tar = ["FSTL_1_A_047", "FSTL_1_B_042", "FSTL_1_C_040", "FSTL_1_D_105", "FSTL_1_E_097", "FSTL_1_F_109",
                    "FSTL_1_H_201", "FSTL_1_I_204", "FSTL_1_K_219", "FSTL_1_*_259", "FSTL_1_L_412", "FSTL_1_M_398",
                    "FSTL_1_O_417", "FSTL_1_P_820"]
for tar in add_dye_test_tar:
    for prim in ["100nM", "250nM", "750nM"]:
        for conc in ["500", "1600"]:
            curDa["Test_ADD_DYE_" + tar + "_" + prim + "_" + conc + "_TD0"] = linRegRes["Probe Mix SYBR"]["P7 - SYBR Primer"][tar + " " + prim + " 1:" + conc]["TD0 mean"]

res = "\nTD0             1:500                                               1:1600\n"
res +=  "                100nM            250nM            750nM             100nM            250nM            750nM\n"
for tar in add_dye_test_tar:
    res += tar.ljust(16)
    for conc in ["500", "1600"]:
        for prim in ["100nM", "250nM", "750nM"]:
            res += colorDiff(curDa, laDa, "Test_ADD_DYE_" + tar + "_" + prim + "_" + conc + "_TD0", "{:6.2f}") + "  "
        res += " "
    res += "\n"
print(res)

ww.close()
rd.save("temp_add_sybr.rdml")


print("\n########################################\n### Test Eva Green - Dye Conc vs Mix ###\n########################################")
print("This test uses the Roche Probe Mastermix and the IDT Probe Mix and adds defined amounts of Eva Green. ")
print("7ng DNA were used 2100 copies should be expected. The difference is mainly due to PCR efficiency ")
print("calculation. ")
print("----------------------")

rd = rdml.Rdml(rdml_eva_a_file)

expList = rd.experiments()
if len(expList) < 1:
    print("No experiments found!")
    sys.exit(0)
ww = open(out_eva_a_file, "w")
startLine = 0
linRegRes = {}
colTar = []
quant = []
for exp in expList:
    runList = exp.runs()
    if len(runList) < 1:
        print("No runs found!")
        sys.exit(0)
    for run in runList:
        if printExpRun:
            print("Experiment: " + exp["id"] + " Run: " + run["id"])
        res = run.webAppLinRegPCR(pcrEfficiencyExl=0.05, updateTargetEfficiency=True, updateRDML=True, excludeNoPlateau=False, excludeEfficiency="outlier", excludeInstableBaseline=False)
        resTab = json.loads(res["LinRegPCR_Result_Table"])
        for tabRow in range(0, len(resTab)):
            if startLine == 0:
                ww.write("Experiment\tRun\t")
                startLine = 1
            else:
                ww.write(exp["id"] + "\t" + run["id"] + "\t")
            for tabCol in range(0, len(resTab[tabRow])):
                outCell = str(resTab[tabRow][tabCol]).replace("\t", ";")
                if tabCol < len(resTab[tabRow]) - 1:
                    ww.write(outCell + "\t")
                else:
                    ww.write(outCell + "\n")
            if startLine == 1:
                if resTab[tabRow][rar_sample_type] in ["unkn", "std"]:
                    if float(resTab[tabRow][rar_Ncopy]) > 5.0:
                        if exp["id"] not in linRegRes:
                            linRegRes[exp["id"]] = {}
                        if run["id"] not in linRegRes[exp["id"]]:
                            linRegRes[exp["id"]][run["id"]] = {}
                        if resTab[tabRow][rar_tar] not in colTar:
                            colTar.append(resTab[tabRow][rar_tar])
                        if resTab[tabRow][rar_tar] not in linRegRes[exp["id"]][run["id"]]:
                            linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]] = {}
                            linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["Ncopy"] = []
                            linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["TD0"] = []
                            linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["PCReff"] = resTab[tabRow][rar_PCR_eff]

                        linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["Ncopy"].append(resTab[tabRow][rar_Ncopy])  # Ncopy
                        linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["TD0"].append(resTab[tabRow][rar_TD0])  # TD0

        for tar in colTar:
            if tar in linRegRes[exp["id"]][run["id"]]:
                linRegRes[exp["id"]][run["id"]][tar]["Ncopy mean"] = np.mean(linRegRes[exp["id"]][run["id"]][tar]["Ncopy"])
                linRegRes[exp["id"]][run["id"]][tar]["TD0 mean"] = np.mean(linRegRes[exp["id"]][run["id"]][tar]["TD0"])    
    
tars = ["FSTL_1_F_109", "FSTL_1_H_201"]
concs = ["L", "M", "H"]
cSt = ["0.125", "0.250", "1.250"]

for tar in tars:
    for conc in concs:
        for mix in ["R", "I"]:
            curDa["Test_EVAGREEN_A_" + tar + "_" + conc + "_" + mix + "_TD0"] = linRegRes["Eva Green A"]["Run 1"][tar + " " + conc + " " + mix]["TD0 mean"]
            curDa["Test_EVAGREEN_A_" + tar + "_" + conc + "_" + mix + "_PCReff"] = linRegRes["Eva Green A"]["Run 1"][tar + " " + conc + " " + mix]["PCReff"]
            curDa["Test_EVAGREEN_A_" + tar + "_" + conc + "_" + mix + "_Ncopy"] = linRegRes["Eva Green A"]["Run 1"][tar + " " + conc + " " + mix]["Ncopy mean"]

res = "\n                                   TD0                                PCReff                                 Ncopy\n"
res +=  "                                   Roche            IDT               Roche              IDT                 Roche              IDT\n"
for tar in tars:
    for cPos in range(0,3):
        res += ("Eva Green " + cSt[cPos] + "uM - " +  tar).ljust(35)
        res += colorDiff(curDa, laDa, "Test_EVAGREEN_A_" + tar + "_" + concs[cPos] + "_R_TD0", "{:6.2f}") + "  "
        res += colorDiff(curDa, laDa, "Test_EVAGREEN_A_" + tar + "_" + concs[cPos] + "_I_TD0", "{:6.2f}") + "   "
        res += colorDiff(curDa, laDa, "Test_EVAGREEN_A_" + tar + "_" + concs[cPos] + "_R_PCReff", "{:7.4f}") + "  "
        res += colorDiff(curDa, laDa, "Test_EVAGREEN_A_" + tar + "_" + concs[cPos] + "_I_PCReff", "{:7.4f}") + "   "
        res += colorDiff(curDa, laDa, "Test_EVAGREEN_A_" + tar + "_" + concs[cPos] + "_R_Ncopy", "{:7.1f}") + "  "
        res += colorDiff(curDa, laDa, "Test_EVAGREEN_A_" + tar + "_" + concs[cPos] + "_I_Ncopy", "{:7.1f}")
        res += "\n"
    res += "\n"
print(res)

ww.close()
rd.save("temp_eva_a.rdml")


print("\n######################################################\n### Test Eva Green - Primer Conc vs Eva Green Conc ###\n######################################################")
print("This test uses the IDT Probe Mix and adds defined amounts of Eva Green. It tests if Eva Green gets ")
print("limiting in low concentartions. ")
print("----------------------")

rd = rdml.Rdml(rdml_eva_b_file)

expList = rd.experiments()
if len(expList) < 1:
    print("No experiments found!")
    sys.exit(0)
ww = open(out_eva_b_file, "w")
startLine = 0
linRegRes = {}
quant = []
for exp in expList:
    runList = exp.runs()
    if len(runList) < 1:
        print("No runs found!")
        sys.exit(0)
    for run in runList:
        if printExpRun:
            print("Experiment: " + exp["id"] + " Run: " + run["id"])
        res = run.webAppLinRegPCR(pcrEfficiencyExl=0.05, updateTargetEfficiency=True, updateRDML=True, excludeNoPlateau=False, excludeEfficiency="outlier", excludeInstableBaseline=False)
        resTab = json.loads(res["LinRegPCR_Result_Table"])
        for tabRow in range(0, len(resTab)):
            if startLine == 0:
                ww.write("Experiment\tRun\t")
                startLine = 1
            else:
                ww.write(exp["id"] + "\t" + run["id"] + "\t")
            for tabCol in range(0, len(resTab[tabRow])):
                outCell = str(resTab[tabRow][tabCol]).replace("\t", ";")
                if tabCol < len(resTab[tabRow]) - 1:
                    ww.write(outCell + "\t")
                else:
                    ww.write(outCell + "\n")
            if startLine == 1:
                if resTab[tabRow][rar_sample_type] in ["unkn", "std"]:
                    if float(resTab[tabRow][rar_Ncopy]) > 5.0:
                        if exp["id"] not in linRegRes:
                            linRegRes[exp["id"]] = {}
                        if run["id"] not in linRegRes[exp["id"]]:
                            linRegRes[exp["id"]][run["id"]] = {}
                        if resTab[tabRow][rar_tar] not in linRegRes[exp["id"]][run["id"]]:
                            linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]] = {}
                        if resTab[tabRow][rar_sample] not in linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]:  # Sample
                            linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]][resTab[tabRow][rar_sample]] = {}
                            linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]][resTab[tabRow][rar_sample]]["Ncopy"] = []
                            linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]][resTab[tabRow][rar_sample]]["TD0"] = []
                            linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]][resTab[tabRow][rar_sample]]["PCReff"] = resTab[tabRow][rar_PCR_eff]

                        linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]][resTab[tabRow][rar_sample]]["Ncopy"].append(resTab[tabRow][rar_Ncopy])  # Ncopy
                        linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]][resTab[tabRow][rar_sample]]["TD0"].append(resTab[tabRow][rar_TD0])  # TD0

        for tar in linRegRes[exp["id"]][run["id"]]:
            for sam in linRegRes[exp["id"]][run["id"]][tar]:
                linRegRes[exp["id"]][run["id"]][tar][sam]["Ncopy mean"] = np.mean(linRegRes[exp["id"]][run["id"]][tar][sam]["Ncopy"])
                linRegRes[exp["id"]][run["id"]][tar][sam]["TD0 mean"] = np.mean(linRegRes[exp["id"]][run["id"]][tar][sam]["TD0"])
    quant = exp.quantify()

['genomic DNA 7ng', 'genomic DNA 20ng', 'genomic DNA 2ng', 'genomic DNA 0.7ng']

tars = ['FSTL_1_A_047', 'FSTL_1_B_042', 'FSTL_1_D_105', 'FSTL_1_E_097', 'FSTL_1_F_109', 'FSTL_1_H_201', 
        'FSTL_1_I_204', 'FSTL_1_K_219', 'FSTL_1_L_412', 'FSTL_1_M_398']

tars = ['FSTL_1_A_047', 'FSTL_1_B_042', 'FSTL_1_D_105', 'FSTL_1_E_097', 'FSTL_1_F_109', 'FSTL_1_H_201']
prims = [" 100nM", ""]
priSt = ["100nM", "250nM"]
concs = [" Low", ""]
cSt = ["0.125", "1.250"]

for tar in tars:
    for cPos in range(0, 2):
        for pPos in range(0, 2):
            if 'genomic DNA 7ng' in linRegRes["Eva Green B"]["Run 1"][tar + prims[pPos] + concs[cPos]]:
                curDa["Test_EVAGREEN_B_" + tar + "_7ng_" + priSt[pPos] + "_" + cSt[cPos] + "_TD0"] = linRegRes["Eva Green B"]["Run 1"][tar + prims[pPos] + concs[cPos]]['genomic DNA 7ng']["TD0 mean"]
                curDa["Test_EVAGREEN_B_" + tar + "_7ng_" + priSt[pPos] + "_" + cSt[cPos] + "_PCReff"] = linRegRes["Eva Green B"]["Run 1"][tar + prims[pPos] + concs[cPos]]['genomic DNA 7ng']["PCReff"]
                curDa["Test_EVAGREEN_B_" + tar + "_7ng_" + priSt[pPos] + "_" + cSt[cPos] + "_Ncopy"] = linRegRes["Eva Green B"]["Run 1"][tar + prims[pPos] + concs[cPos]]['genomic DNA 7ng']["Ncopy mean"]
            else:
                print("--" + tar + prims[pPos] + concs[cPos] + "---")

res = "\n                                   TD0                               PCReff                                Ncopy\n"
res +=  "Eva Green                          0.125            1.250            0.125              1.250              0.125              1.250             \n"
for pPos in range(0, 2):
    for tar in tars:
        res += ("Primer Conc " + priSt[pPos] + " - " +  tar).ljust(35)
        res += colorDiff(curDa, laDa, "Test_EVAGREEN_B_" + tar + "_7ng_" + priSt[pPos] + "_0.125_TD0", "{:6.2f}") + "  "
        res += colorDiff(curDa, laDa, "Test_EVAGREEN_B_" + tar + "_7ng_" + priSt[pPos] + "_1.250_TD0", "{:6.2f}") + "  "
        res += colorDiff(curDa, laDa, "Test_EVAGREEN_B_" + tar + "_7ng_" + priSt[pPos] + "_0.125_PCReff", "{:7.4f}") + "  "
        res += colorDiff(curDa, laDa, "Test_EVAGREEN_B_" + tar + "_7ng_" + priSt[pPos] + "_1.250_PCReff", "{:7.4f}") + "  "
        res += colorDiff(curDa, laDa, "Test_EVAGREEN_B_" + tar + "_7ng_" + priSt[pPos] + "_0.125_Ncopy", "{:7.1f}") + "  "
        res += colorDiff(curDa, laDa, "Test_EVAGREEN_B_" + tar + "_7ng_" + priSt[pPos] + "_1.250_Ncopy", "{:7.1f}") + "  "
        res += "\n"
    res += "\n"
print(res)



head_var = ["curve_pcr_eff", "dilution_pcr_eff", "expected_max", "mean_max", "expected_min", "mean_min", 
            "expected_ratio", "mean_ratio", "slope_bias", "correlation_R", "linearity", 
            "reproducibility", "detectable_diff"]
head_a = ["curve_pcr_eff", "dilution_pcr_eff", "expected_max", "mean_max", "expected_min", "mean_min", 
            "expected_ratio", "mean_ratio"]
head_b = ["slope_bias", "correlation_R", "linearity", "reproducibility", "detectable_diff"]
prec_a = [0.001, 0.001, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
prec_b = [0.0001, 0.00001, 0.00001, 0.00001, 0.0001]
print_a = [3, 3, 0, 0, 0, 0, 0, 0]
print_b = [5, 5, 5, 5, 5]
tars = ['FSTL_1_A_047', 'FSTL_1_B_042', 'FSTL_1_D_105', 'FSTL_1_E_097', 'FSTL_1_F_109', 'FSTL_1_H_201', 
        'FSTL_1_I_204', 'FSTL_1_K_219', 'FSTL_1_L_412', 'FSTL_1_M_398']

for tar in tars:
    for var in head_var:
        keyVal = "Test_EVAGREEN_B_DIL_" + tar.replace(" ", "_") + "_" + var
        curDa[keyVal] = quant["dil_std"][tar][var]

print("\n##########################################\n### Test Eva Green - Different Primers ###\n##########################################")
print("This test uses the IDT Probe Mix and adds defined amounts of Eva Green. It tests all primer pais in ")
print("normal Eva Green concentartions. ")
print("----------------------")


printTable(30, "Test_EVAGREEN_B_DIL_", tars, head_a, prec_a, print_a, curDa, laDa)
print(" ")
printTable(30, "Test_EVAGREEN_B_DIL_", tars, head_b, prec_b, print_b, curDa, laDa)


ww.close()
rd.save("temp_eva_a.rdml")


print("\n######################\n### Test Vermeulen ###\n######################")
print("This test uses the Vermeulen data. It first calculates everything.")
print("----------------------")
print("  Testing ony standards.")
# The standards only
verm_res = {}
vermeulen_stds = {"vermeulen_std_5": vermeulen_std_5,
                  "vermeulen_std_4": vermeulen_std_4,
                  "vermeulen_std_3": vermeulen_std_3,
                  "vermeulen_std_2": vermeulen_std_2,
                  "vermeulen_std_1": vermeulen_std_1}
for std_file in vermeulen_stds:
    if std_file not in verm_res:
        verm_res[std_file] = {}
    rdd = rdml.Rdml(vermeulen_stds[std_file])
    expList = rdd.experiments()
    ww = open("temp_" + std_file + ".csv", "w")
    startLine = 0
    reactionDataTrue = 0
    reactionDataFalse = 0
    for exp in expList:
        runList = exp.runs()
        for run in runList:
            res = run.webAppLinRegPCR(pcrEfficiencyExl=0.05, updateTargetEfficiency=False, updateRDML=True, excludeNoPlateau=True, excludeEfficiency="outlier", excludeInstableBaseline=True)
            resTab = json.loads(res["LinRegPCR_Result_Table"])
            for tabRow in range(startLine, len(resTab)):
                if startLine == 0:
                    ww.write("Experiment\tRun\t")
                    startLine = 1
                else:
                    ww.write(exp["id"] + "\t" + run["id"] + "\t")
                for tabCol in range(0, len(resTab[tabRow])):
                    outCell = str(resTab[tabRow][tabCol]).replace("\t", ";")
                    if tabCol < len(resTab[tabRow]) - 1:
                        ww.write(outCell + "\t")
                    else:
                        ww.write(outCell + "\n")
                # save the data
                if tabRow > 0:
                    if resTab[tabRow][rar_tar] not in verm_res[std_file]:
                        verm_res[std_file][resTab[tabRow][rar_tar]] = {}
                    if resTab[tabRow][rar_sample] not in verm_res[std_file][resTab[tabRow][rar_tar]]:
                        verm_res[std_file][resTab[tabRow][rar_tar]][resTab[tabRow][rar_sample ]] = {}
                        verm_res[std_file][resTab[tabRow][rar_tar]][resTab[tabRow][rar_sample ]]["mean Eff"] = []
                        verm_res[std_file][resTab[tabRow][rar_tar]][resTab[tabRow][rar_sample ]]["Ncopy"] = []
                    verm_res[std_file][resTab[tabRow][rar_tar]][resTab[tabRow][rar_sample ]]["mean Eff"].append(resTab[tabRow][rar_PCR_eff])
                    verm_res[std_file][resTab[tabRow][rar_tar]][resTab[tabRow][rar_sample ]]["Ncopy"].append(resTab[tabRow][rar_Ncopy])
    ww.close()
print("    Done.\n")

# Time the test
std_file = "vermeulen_all"
verm_res[std_file] = {}
startTime = time.time()
rdd = rdml.Rdml(vermeulen_file)

expList = rdd.experiments()
print("  Testing all data will run for up to 2 Minutes.")
if len(expList) < 1:
    print("No experiments found!")
    sys.exit(0)
ww = open(out_file, "w")
startLine = 0
reactionDataTrue = 0
reactionDataFalse = 0
for exp in expList:
    runList = exp.runs()
    if len(runList) < 1:
        print("No runs found!")
        sys.exit(0)
    for run in runList:
        if printExpRun:
            print("      Experiment: " + exp["id"] + " Run: " + run["id"])
        res = run.webAppLinRegPCR(pcrEfficiencyExl=0.05, updateTargetEfficiency=False, updateRDML=True, excludeNoPlateau=True, excludeEfficiency="outlier", excludeInstableBaseline=True)
        resTab = json.loads(res["LinRegPCR_Result_Table"])
        for tabRow in range(startLine, len(resTab)):
            if startLine == 0:
                ww.write("Experiment\tRun\t")
                startLine = 1
            else:
                ww.write(exp["id"] + "\t" + run["id"] + "\t")
            for tabCol in range(0, len(resTab[tabRow])):
                outCell = str(resTab[tabRow][tabCol]).replace("\t", ";")
                if tabCol < len(resTab[tabRow]) - 1:
                    ww.write(outCell + "\t")
                else:
                    ww.write(outCell + "\n")
            if startLine == 1:
                if resTab[tabRow][rar_sample_type] in ["unkn", "std"]:
                    if float(resTab[tabRow][rar_Ncopy]) > 5.0:
                        reactionDataTrue += 1
                    else:
                        reactionDataFalse += 1
            if tabRow > 0:
                if resTab[tabRow][rar_sample].startswith("STD_"):
                    if resTab[tabRow][rar_tar] not in verm_res[std_file]:
                        verm_res[std_file][resTab[tabRow][rar_tar]] = {}
                    if resTab[tabRow][rar_sample] not in verm_res[std_file][resTab[tabRow][rar_tar]]:
                        verm_res[std_file][resTab[tabRow][rar_tar]][resTab[tabRow][rar_sample ]] = {}
                        verm_res[std_file][resTab[tabRow][rar_tar]][resTab[tabRow][rar_sample ]]["mean Eff"] = []
                        verm_res[std_file][resTab[tabRow][rar_tar]][resTab[tabRow][rar_sample ]]["Ncopy"] = []
                    verm_res[std_file][resTab[tabRow][rar_tar]][resTab[tabRow][rar_sample ]]["mean Eff"].append(resTab[tabRow][rar_PCR_eff])
                    verm_res[std_file][resTab[tabRow][rar_tar]][resTab[tabRow][rar_sample ]]["Ncopy"].append(resTab[tabRow][rar_Ncopy])
print("    Done.\n")
ww.close()
rdd.save("temp_vermeulen_linregpcr.rdml")

curDa["reactionDataFalse"] = reactionDataFalse
curDa["reactionDataSum"] = reactionDataFalse + reactionDataTrue
printNice("Failing Reactions: ", curDa["reactionDataFalse"], laDa["reactionDataFalse"], 1, "-" , 0)
printNice("Sum Reactions: ", curDa["reactionDataSum"], laDa["reactionDataSum"], 1, "+", 0)

# Write the eff an Ncopy comparison
print("\n######################\n### Test PCR eff ###\n######################")
print("This test uses the Vermeulen data. The upper panel calculates the PCR efficiency and ")
print("the lower panel calculates the coressponding Ncopy of the highest concentration. ")
print("The left value uses all reactions of the 384 plate. The STD 5 only all standards ")
print("(15 reactions) downt to SDT 1 with only the 3 replicates. On the right side the ")
print("difference (upper) or the factor (lower) of each value to STD 5 is given.")
print("----------------------")

print("  Testing how the reduction of reactions affects PCR efficiency and Ncopy:")
vm_eff_diff = {}
vm_eff_val = {}
for cVer in verm_res:
    for tar in verm_res[cVer]:
        for sam in verm_res[cVer][tar]:
            verm_res[cVer][tar][sam]["calc mean Eff"] = np.mean(verm_res[cVer][tar][sam]["mean Eff"])
            verm_res[cVer][tar][sam]["calc Ncopy"] = np.mean(verm_res[cVer][tar][sam]["Ncopy"])
out_st = "\n"
out_st += "Target".ljust(23) + "all".ljust(8) + "STD 5".ljust(8) + "STD 4".ljust(8) + "STD 3".ljust(8) + "STD 2".ljust(8) + "STD 1".ljust(8)
out_st += " ".ljust(2) + "all".ljust(8) + "STD 5".ljust(8) + "STD 4".ljust(8) + "STD 3".ljust(8) + "STD 2".ljust(8) + "STD 1".ljust(8) + "\n"
out_st += "".ljust(23) + "PCR eff".ljust(8) + "PCR eff".ljust(8) + "PCR eff".ljust(8) + "PCR eff".ljust(8) + "PCR eff".ljust(8) + "PCR eff".ljust(8)
out_st += " ".ljust(2) + "diff".ljust(8) + "diff".ljust(8) + "diff".ljust(8) + "diff".ljust(8) + "diff".ljust(8) + "diff".ljust(8) + "\n"

for tar in verm_res[std_file]:
    if tar == "ALUsq(Eurogentec)_20080228":
        continue
    out_st += tar.ljust(22)
    for cVer in ["vermeulen_all", "vermeulen_std_5", "vermeulen_std_4", 
                 "vermeulen_std_3", "vermeulen_std_2", "vermeulen_std_1"]:
        out_st += "{:8.4f}".format(verm_res[cVer][tar]["STD_150000"]["calc mean Eff"])
    out_st += "  "
    for cVer in ["vermeulen_all", "vermeulen_std_5", "vermeulen_std_4", 
                 "vermeulen_std_3", "vermeulen_std_2", "vermeulen_std_1"]:
        if cVer not in vm_eff_diff:
            vm_eff_diff[cVer] = []
        if cVer not in vm_eff_val:
            vm_eff_val[cVer] = []
        vm_eff_val[cVer].append(verm_res[cVer][tar]["STD_150000"]["calc mean Eff"])
        vm_eff_diff[cVer].append(verm_res[cVer][tar]["STD_150000"]["calc mean Eff"] - verm_res["vermeulen_std_5"][tar]["STD_150000"]["calc mean Eff"])
        out_st += "{:8.4f}".format(verm_res[cVer][tar]["STD_150000"]["calc mean Eff"] - verm_res["vermeulen_std_5"][tar]["STD_150000"]["calc mean Eff"])
    out_st += "\n"
out_st += "Mean".ljust(22)
for cVer in ["vermeulen_all", "vermeulen_std_5", "vermeulen_std_4", 
                "vermeulen_std_3", "vermeulen_std_2", "vermeulen_std_1"]:
    curDa[cVer + "_eff_val_mean"] = np.nanmean(vm_eff_val[cVer])
    out_st += stringNice2(curDa, laDa, cVer + "_eff_val_mean", "{:8.4f}", 0.0001, "+", True)
out_st += " ".ljust(2)
for cVer in ["vermeulen_all", "vermeulen_std_5", "vermeulen_std_4", 
                "vermeulen_std_3", "vermeulen_std_2", "vermeulen_std_1"]:
    curDa[cVer + "_eff_dif_mean"] = np.nanmean(vm_eff_diff[cVer])
    out_st += stringNice2(curDa, laDa, cVer + "_eff_dif_mean", "{:8.4f}", 0.0001, "+", True)
out_st += "\n"
out_st += "Stored Mean".ljust(22)
for cVer in ["vermeulen_all", "vermeulen_std_5", "vermeulen_std_4", 
                "vermeulen_std_3", "vermeulen_std_2", "vermeulen_std_1"]:
    out_st += stringNice2(curDa, laDa, cVer + "_eff_val_mean", "{:8.4f}", 0.0001, "+", False)
out_st += " ".ljust(2)
for cVer in ["vermeulen_all", "vermeulen_std_5", "vermeulen_std_4", 
                "vermeulen_std_3", "vermeulen_std_2", "vermeulen_std_1"]:
    out_st += stringNice2(curDa, laDa, cVer + "_eff_dif_mean", "{:8.4f}", 0.0001, "+", False)
out_st += "\n"
out_st += "STD".ljust(22)
for cVer in ["vermeulen_all", "vermeulen_std_5", "vermeulen_std_4", 
                "vermeulen_std_3", "vermeulen_std_2", "vermeulen_std_1"]:
    curDa[cVer + "_eff_val_std"] = np.nanstd(vm_eff_val[cVer])
    out_st += stringNice2(curDa, laDa, cVer + "_eff_val_std", "{:8.4f}", 0.0001, "+", True)
out_st += " ".ljust(2)
for cVer in ["vermeulen_all", "vermeulen_std_5", "vermeulen_std_4", 
                "vermeulen_std_3", "vermeulen_std_2", "vermeulen_std_1"]:
    curDa[cVer + "_eff_dif_std"] = np.nanstd(vm_eff_diff[cVer])
    out_st += stringNice2(curDa, laDa, cVer + "_eff_dif_std", "{:8.4f}", 0.0001, "+", True)
out_st += "\n"
out_st += "Stored STD".ljust(22)
for cVer in ["vermeulen_all", "vermeulen_std_5", "vermeulen_std_4", 
                "vermeulen_std_3", "vermeulen_std_2", "vermeulen_std_1"]:
    out_st += stringNice2(curDa, laDa, cVer + "_eff_val_std", "{:8.4f}", 0.0001, "+", False)
out_st += " ".ljust(2)
for cVer in ["vermeulen_all", "vermeulen_std_5", "vermeulen_std_4", 
                "vermeulen_std_3", "vermeulen_std_2", "vermeulen_std_1"]:
    out_st += stringNice2(curDa, laDa, cVer + "_eff_dif_std", "{:8.4f}", 0.0001, "+", False)
out_st += "\n"

print(out_st)

vm_eff_diff = {}
vm_eff_val = {}
out_st = "\n"
out_st += "Target".ljust(23) + "all".ljust(9) + "STD 5".ljust(9) + "STD 4".ljust(9) + "STD 3".ljust(9) + "STD 2".ljust(9) + "STD 1".ljust(9)
out_st += " ".ljust(2) + "all".ljust(8) + "STD 5".ljust(8) + "STD 4".ljust(8) + "STD 3".ljust(8) + "STD 2".ljust(8) + "STD 1".ljust(8) + "\n"
out_st += "".ljust(23) + "Ncopy".ljust(9) + "Ncopy".ljust(9) + "Ncopy".ljust(9) + "Ncopy".ljust(9) + "Ncopy".ljust(9) + "Ncopy".ljust(9)
out_st += " ".ljust(2) + "fact".ljust(8) + "fact".ljust(8) + "fact".ljust(8) + "fact".ljust(8) + "fact".ljust(8) + "fact".ljust(8) + "\n"

for tar in verm_res[std_file]:
    if tar == "ALUsq(Eurogentec)_20080228":
        continue
   # if tar == "CLSTN1_20080227":
   #     continue
    out_st += tar.ljust(22)
    for cVer in ["vermeulen_all", "vermeulen_std_5", "vermeulen_std_4", 
                 "vermeulen_std_3", "vermeulen_std_2", "vermeulen_std_1"]:
        out_st += "{:9.0f}".format(verm_res[cVer][tar]["STD_150000"]["calc Ncopy"])
    out_st += "  "
    for cVer in ["vermeulen_all", "vermeulen_std_5", "vermeulen_std_4", 
                 "vermeulen_std_3", "vermeulen_std_2", "vermeulen_std_1"]:
        if cVer not in vm_eff_diff:
            vm_eff_diff[cVer] = []
        if cVer not in vm_eff_val:
            vm_eff_val[cVer] = []
        vm_eff_val[cVer].append(verm_res[cVer][tar]["STD_150000"]["calc Ncopy"])
        vm_eff_diff[cVer].append(verm_res[cVer][tar]["STD_150000"]["calc Ncopy"] / verm_res["vermeulen_std_5"][tar]["STD_150000"]["calc Ncopy"])
        out_st += "{:8.4f}".format(verm_res[cVer][tar]["STD_150000"]["calc Ncopy"] / verm_res["vermeulen_std_5"][tar]["STD_150000"]["calc Ncopy"])
    out_st += "\n"
out_st += "Mean".ljust(22)
for cVer in ["vermeulen_all", "vermeulen_std_5", "vermeulen_std_4", 
                "vermeulen_std_3", "vermeulen_std_2", "vermeulen_std_1"]:
    curDa[cVer + "_ncopy_val_mean"] = np.nanmean(vm_eff_val[cVer])
    out_st += stringNice2(curDa, laDa, cVer + "_ncopy_val_mean", "{:9.0f}", 1.0, "+", True)

out_st += " ".ljust(2)
for cVer in ["vermeulen_all", "vermeulen_std_5", "vermeulen_std_4", 
                "vermeulen_std_3", "vermeulen_std_2", "vermeulen_std_1"]:
    curDa[cVer + "_ncopy_dif_mean"] = np.nanmean(vm_eff_diff[cVer])
    out_st += stringNice2(curDa, laDa, cVer + "_ncopy_dif_mean", "{:8.4f}", 0.0001, "+", True)
out_st += "\n"
out_st += "Stored Mean".ljust(22)
for cVer in ["vermeulen_all", "vermeulen_std_5", "vermeulen_std_4", 
                "vermeulen_std_3", "vermeulen_std_2", "vermeulen_std_1"]:
    out_st += stringNice2(curDa, laDa, cVer + "_ncopy_val_mean", "{:9.0f}", 1.0, "+", False)
out_st += " ".ljust(2)
for cVer in ["vermeulen_all", "vermeulen_std_5", "vermeulen_std_4", 
                "vermeulen_std_3", "vermeulen_std_2", "vermeulen_std_1"]:
    out_st += stringNice2(curDa, laDa, cVer + "_ncopy_dif_mean", "{:8.4f}", 0.0001, "+", False)
out_st += "\n"
out_st += "STD".ljust(22)
for cVer in ["vermeulen_all", "vermeulen_std_5", "vermeulen_std_4", 
                "vermeulen_std_3", "vermeulen_std_2", "vermeulen_std_1"]:
    curDa[cVer + "_ncopy_val_std"] = np.nanstd(vm_eff_val[cVer])
    out_st += stringNice2(curDa, laDa, cVer + "_ncopy_val_std","{:9.0f}", 1.0, "+", True)
out_st += " ".ljust(2)
for cVer in ["vermeulen_all", "vermeulen_std_5", "vermeulen_std_4", 
                "vermeulen_std_3", "vermeulen_std_2", "vermeulen_std_1"]:
    curDa[cVer + "_ncopy_dif_std"] = np.nanstd(vm_eff_diff[cVer])
    out_st += stringNice2(curDa, laDa, cVer + "_ncopy_dif_std", "{:8.4f}", 0.0001, "+", True)
out_st += "\n"
out_st += "Stored STD".ljust(22)
for cVer in ["vermeulen_all", "vermeulen_std_5", "vermeulen_std_4", 
                "vermeulen_std_3", "vermeulen_std_2", "vermeulen_std_1"]:
    out_st += stringNice2(curDa, laDa, cVer + "_ncopy_val_std","{:9.0f}", 1.0, "+", False)
out_st += " ".ljust(2)
for cVer in ["vermeulen_all", "vermeulen_std_5", "vermeulen_std_4", 
                "vermeulen_std_3", "vermeulen_std_2", "vermeulen_std_1"]:
    out_st += stringNice2(curDa, laDa, cVer + "_ncopy_dif_std", "{:8.4f}", 0.0001, "+", False)
out_st += "\n"
print(out_st)

endTime = time.time()
runTime = endTime - startTime
runMin = math.floor(runTime / 60.0)
runSec = runTime - runMin* 60.0
print("\n######################\n### Test Vermeulen ###\n######################")
print("This test uses the Vermeulen data. The upper panel calculates the test like published (by order): ")
print("[Methods. 2013 Jan;59(1):32-46. doi: 10.1016/j.ymeth.2012.08.011.")
print("The lower panel calculates result as mean. ")
print("----------------------")
print("Runtime: " + str(runMin) + ":" + "{:.3f}".format(runSec))
json_object = json.dumps(curDa, indent=4)
with open(out_json, "w") as outfile:
    outfile.write(json_object)
print(" ")
os.system("python3 vermeulen_analyze_csv.py")
os.system("python3 vermeulen_analyze_compare.py")
sys.exit(0)
