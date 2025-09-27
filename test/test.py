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
rdml_mach_file = os.path.join(parent_dir, "experiments/untergasser/volume_machine.rdml")
rdml_amp_file = os.path.join(parent_dir, "experiments/untergasser/amplicon_primer_mix.rdml")
out_amp_file = "temp_amplicon_primer_resuls.csv"
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
                colSavSum[count][roNum] = sav[keyVal]
                direct = cur[keyVal] - sav[keyVal]
                diff = abs(cur[keyVal] - sav[keyVal])
                if diff > prec[count]:
                    if direct > 0.0:
                        res +=  '\033[44m'
                    else:
                        res +=  '\033[41m'
                if prin[count] == 0:
                    res += "{:8.0f}".format(cur[keyVal]) + " (" + "{:8.0f}".format(sav[keyVal]) + ") "
                elif prin[count] == 1:
                    res += "{:8.1f}".format(cur[keyVal]) + " (" + "{:8.1f}".format(sav[keyVal]) + ") "
                elif prin[count] == 2:
                    res += "{:8.2f}".format(cur[keyVal]) + " (" + "{:8.2f}".format(sav[keyVal]) + ") "
                elif prin[count] == 3:
                    res += "{:8.3f}".format(cur[keyVal]) + " (" + "{:8.3f}".format(sav[keyVal]) + ") "
                else:
                    res += "{:8.5f}".format(cur[keyVal]) + " (" + "{:8.5f}".format(sav[keyVal]) + ") "
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
        if keyVal in cur and keyVal in sav:
            direct = cur[keyVal] - sav[keyVal]
            diff = abs(cur[keyVal] - sav[keyVal])
            if diff > 0.01:
                if direct > 0.0:
                    res +=  '\033[44m'
                else:
                    res +=  '\033[41m'
            res += "{:8.5f}".format(cur[keyVal]) + " (" + "{:8.5f}".format(sav[keyVal]) + ") "
            if diff > 0.01:
                res +=  '\033[0m'
        else:
            res += "".ljust(20)
        print(res)


def printTable(rows, cols, prec, prin, cur, sav):
    res = "".ljust(20)
    for col in cols:
        res += col.ljust(20)
    print(res)
    colCurSum = np.zeros([len(cols), len(rows)], dtype=np.float64)
    colSavSum = np.zeros([len(cols), len(rows)], dtype=np.float64)
    colCurSum[:] = np.nan
    colSavSum[:] = np.nan
    roNum = 0
    for row in rows:
        res = row.ljust(20)
        count = 0
        for col in cols:
            keyVal = row.replace(" ", "_") + "_" + col
            colCurSum[count][roNum] = cur[keyVal]
            colSavSum[count][roNum] = sav[keyVal]
            direct = cur[keyVal] - sav[keyVal]
            diff = abs(cur[keyVal] - sav[keyVal])
            if diff > prec[count]:
                if direct > 0.0:
                    res +=  '\033[44m'
                else:
                    res +=  '\033[41m'
            if prin[count] == 0:
                res += "{:8.0f}".format(cur[keyVal]) + " (" + "{:8.0f}".format(sav[keyVal]) + ") "
            elif prin[count] == 1:
                res += "{:8.1f}".format(cur[keyVal]) + " (" + "{:8.1f}".format(sav[keyVal]) + ") "
            elif prin[count] == 2:
                res += "{:8.2f}".format(cur[keyVal]) + " (" + "{:8.2f}".format(sav[keyVal]) + ") "
            elif prin[count] == 3:
                res += "{:8.3f}".format(cur[keyVal]) + " (" + "{:8.3f}".format(sav[keyVal]) + ") "
            else:
                res += "{:8.5f}".format(cur[keyVal]) + " (" + "{:8.5f}".format(sav[keyVal]) + ") "
            if diff > prec[count]:
                res +=  '\033[0m'
            count += 1
        print(res)
        roNum += 1
    res = "Mean".ljust(20)
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
    res = "CV".ljust(20)
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
        res = "Diff cur dil".ljust(20)
        std_ar = ["curve_pcr_eff", "dilution_pcr_eff"]
        sum_eff_cur = 0.0
        sum_eff_sav = 0.0
        num_eff = 0.0
        for row in rows:
            keyCurv = row.replace(" ", "_") + "_curve_pcr_eff"
            keyDil = row.replace(" ", "_") + "_dilution_pcr_eff"
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

print("\n##############################\n### Remove Temporary Files ###\n###############################")
print('rm ' + os.path.join(parent_dir, "test/temp_*"))
os.system('rm ' + os.path.join(parent_dir, "test/temp_*"))

print("\n#######################\n### Technical Tests ###\n#######################")

print("Test: Fill Matrix Gaps:")

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
# Time the test
startTime = time.time()
rd = rdml.Rdml(rdml_mach_file)

expList = rd.experiments()
if len(expList) < 1:
    print("No experiments found!")
    sys.exit(0)

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
                        
                        linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]][resTab[tabRow][rar_sample]]["Ncopy"].append(resTab[tabRow][rar_Ncopy])  # Ncopy

        for tar in colTar:
            for sam in linRegRes[exp["id"]][run["id"]][tar]:
                linRegRes[exp["id"]][run["id"]][tar][sam]["Ncopy mean"] = np.mean(linRegRes[exp["id"]][run["id"]][tar][sam]["Ncopy"])

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
print("  5ul - 0.5ng/ul -  750 copies: " + "{:6.1f}".format(curDa["Test_Vol_AMC_384_05"]) + " (" + "{:6.1f}".format(saveNum(laDa, "Test_Vol_AMC_384_05")) + ")" )
print(" 10ul - 0.5ng/ul - 1500 copies: " + "{:6.1f}".format(curDa["Test_Vol_AMC_384_10"]) + " (" + "{:6.1f}".format(saveNum(laDa, "Test_Vol_AMC_384_10")) + ")" )
print(" 20ul - 0.5ng/ul - 3000 copies: " + "{:6.1f}".format(curDa["Test_Vol_AMC_384_20"]) + " (" + "{:6.1f}".format(saveNum(laDa, "Test_Vol_AMC_384_20")) + ")" )
print("\nEMBL Heidelberg DE - 384 Wells")
print("  5ul - 0.5ng/ul -  750 copies: " + "{:6.1f}".format(curDa["Test_Vol_EMBL_384_05"]) + " (" + "{:6.1f}".format(saveNum(laDa, "Test_Vol_EMBL_384_05")) + ")" )
print(" 10ul - 0.5ng/ul - 1500 copies: " + "{:6.1f}".format(curDa["Test_Vol_EMBL_384_10"]) + " (" + "{:6.1f}".format(saveNum(laDa, "Test_Vol_EMBL_384_10")) + ")" )
print(" 20ul - 0.5ng/ul - 3000 copies: " + "{:6.1f}".format(curDa["Test_Vol_EMBL_384_20"]) + " (" + "{:6.1f}".format(saveNum(laDa, "Test_Vol_EMBL_384_20")) + ")" )
print("\nEMBL Heidelberg DE - 96 Wells")
print(" 10ul - 0.5ng/ul - 1500 copies: " + "{:6.1f}".format(curDa["Test_Vol_EMBL_96_10"]) + " (" + "{:6.1f}".format(saveNum(laDa, "Test_Vol_EMBL_96_10")) + ")" )
print(" 20ul - 0.5ng/ul - 3000 copies: " + "{:6.1f}".format(curDa["Test_Vol_EMBL_96_20"]) + " (" + "{:6.1f}".format(saveNum(laDa, "Test_Vol_EMBL_96_20")) + ")" )
print(" 40ul - 0.5ng/ul - 6000 copies: " + "{:6.1f}".format(curDa["Test_Vol_EMBL_96_40"]) + " (" + "{:6.1f}".format(saveNum(laDa, "Test_Vol_EMBL_96_40")) + ")" )

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
print(" 10ul - 0.2ng/ul -   600 copies: " + "{:8.1f}".format(curDa["Test_Vol_AMC_384_c0.2"]) + " (" + "{:8.1f}".format(saveNum(laDa, "Test_Vol_AMC_384_c0.2")) + ")" )
print(" 10ul - 0.5ng/ul -  1500 copies: " + "{:8.1f}".format(curDa["Test_Vol_AMC_384_c0.5"]) + " (" + "{:8.1f}".format(saveNum(laDa, "Test_Vol_AMC_384_c0.5")) + ")" )
print(" 10ul - 1.0ng/ul -  3000 copies: " + "{:8.1f}".format(curDa["Test_Vol_AMC_384_c1"]) + " (" + "{:8.1f}".format(saveNum(laDa, "Test_Vol_AMC_384_c1")) + ")" )
print(" 10ul - 2.0ng/ul -  6000 copies: " + "{:8.1f}".format(curDa["Test_Vol_AMC_384_c2"]) + " (" + "{:8.1f}".format(saveNum(laDa, "Test_Vol_AMC_384_c2")) + ")" )
print(" 10ul - 4.0ng/ul - 12000 copies: " + "{:8.1f}".format(curDa["Test_Vol_AMC_384_c4"]) + " (" + "{:8.1f}".format(saveNum(laDa, "Test_Vol_AMC_384_c4")) + ")" )
print(" 10ul - 8.0ng/ul - 24000 copies: " + "{:8.1f}".format(curDa["Test_Vol_AMC_384_c8"]) + " (" + "{:8.1f}".format(saveNum(laDa, "Test_Vol_AMC_384_c8")) + ")" )
print("\nEMBL Heidelberg DE - 384 Wells")
print(" 10ul - 0.2ng/ul -   600 copies: " + "{:8.1f}".format(curDa["Test_Vol_EMBL_384_c0.2"]) + " (" + "{:8.1f}".format(saveNum(laDa, "Test_Vol_EMBL_384_c0.2")) + ")" )
print(" 10ul - 0.5ng/ul -  1500 copies: " + "{:8.1f}".format(curDa["Test_Vol_EMBL_384_c0.5"]) + " (" + "{:8.1f}".format(saveNum(laDa, "Test_Vol_EMBL_384_c0.5")) + ")" )
print(" 10ul - 1.0ng/ul -  3000 copies: " + "{:8.1f}".format(curDa["Test_Vol_EMBL_384_c1"]) + " (" + "{:8.1f}".format(saveNum(laDa, "Test_Vol_EMBL_384_c1")) + ")" )
print(" 10ul - 2.0ng/ul -  6000 copies: " + "{:8.1f}".format(curDa["Test_Vol_EMBL_384_c2"]) + " (" + "{:8.1f}".format(saveNum(laDa, "Test_Vol_EMBL_384_c2")) + ")" )
print(" 10ul - 4.0ng/ul - 12000 copies: " + "{:8.1f}".format(curDa["Test_Vol_EMBL_384_c4"]) + " (" + "{:8.1f}".format(saveNum(laDa, "Test_Vol_EMBL_384_c4")) + ")" )
print("\nEMBL Heidelberg DE - 96 Wells")
print(" 20ul - 0.2ng/ul -  1200 copies: " + "{:8.1f}".format(curDa["Test_Vol_EMBL_96_c0.2"]) + " (" + "{:8.1f}".format(saveNum(laDa, "Test_Vol_EMBL_96_c0.2")) + ")" )
print(" 20ul - 0.5ng/ul -  3000 copies: " + "{:8.1f}".format(curDa["Test_Vol_EMBL_96_c0.5"]) + " (" + "{:8.1f}".format(saveNum(laDa, "Test_Vol_EMBL_96_c0.5")) + ")" )
print(" 20ul - 1.0ng/ul -  6000 copies: " + "{:8.1f}".format(curDa["Test_Vol_EMBL_96_c1"]) + " (" + "{:8.1f}".format(saveNum(laDa, "Test_Vol_EMBL_96_c1")) + ")" )
print(" 20ul - 2.0ng/ul - 12000 copies: " + "{:8.1f}".format(curDa["Test_Vol_EMBL_96_c2"]) + " (" + "{:8.1f}".format(saveNum(laDa, "Test_Vol_EMBL_96_c2")) + ")" )
print(" 20ul - 4.0ng/ul - 24000 copies: " + "{:8.1f}".format(curDa["Test_Vol_EMBL_96_c4"]) + " (" + "{:8.1f}".format(saveNum(laDa, "Test_Vol_EMBL_96_c4")) + ")" )

print("\n#############################\n### Test Primer Dilutions ###\n#############################")
# Time the test
startTime = time.time()
rd = rdml.Rdml(rdml_amp_file)

expList = rd.experiments()
if len(expList) < 1:
    print("No experiments found!")
    sys.exit(0)
ww = open(out_amp_file, "w")
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
                            linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["PCR Eff indiv"] = []
                            linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["Ncopy indiv"] = []
                            linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["TD0"] = []
                        
                        linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["PCR Eff indiv"].append(resTab[tabRow][rar_indiv_PCR_eff])
                        linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["Ncopy indiv"].append(resTab[tabRow][rar_indiv_Ncopy])
                        linRegRes[exp["id"]][run["id"]][resTab[tabRow][rar_tar]]["TD0"].append(resTab[tabRow][rar_TD0])

                        reactionDataTrue += 1
                    else:
                        reactionDataFalse += 1

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

        print("Roche:")
        printPrimerTable(tars, "", curDa, laDa)
        print("\nSensi:")
        printPrimerTable(tars, "_SF", curDa, laDa)
        print("\nLC-Green:")
        printPrimerTable(tars, "_LC", curDa, laDa)
        
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
        print("Roche:")
        printCVTable(tar_roche, "", curDa, laDa)
        print("\nSensi:")
        printCVTable(tar_roche, "_SF", curDa, laDa)
        print("\nLC-Green:")
        printCVTable(tar_roche, "_LC", curDa, laDa)

        print("\n\nRoche:")
        printTable(tar_roche, head_a, prec_a, print_a, curDa, laDa)
        print("\nSensi:")
        printTable(tar_sensi, head_a, prec_a, print_a, curDa, laDa)
        print("\nLC-Green:")
        printTable(tar_lc, head_a, prec_a, print_a, curDa, laDa)
        print("\n\nRoche:")
        printTable(tar_roche, head_b, prec_b, print_b, curDa, laDa)
        print("\nSensi:")
        printTable(tar_sensi, head_b, prec_b, print_b, curDa, laDa)
        print("\nLC-Green:")
        printTable(tar_lc, head_b, prec_b, print_b, curDa, laDa)

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

print("\n######################\n### Test Vermeulen ###\n######################")
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
print("Runtime: " + str(runMin) + ":" + "{:.3f}".format(runSec))
json_object = json.dumps(curDa, indent=4)
with open(out_json, "w") as outfile:
    outfile.write(json_object)
print(" ")
os.system("python3 vermeulen_analyze_csv.py")
os.system("python3 vermeulen_analyze_compare.py")
sys.exit(0)
