#!/usr/bin/python

import csv
import sys

if len(sys.argv) != 3:
    print("compare_extra_columns.py requires two files to compare. Example use: ")
    print("   python3 compare_extra_columns.py temp_5_out_results.tsv test_5_out_results.tsv")
    sys.exit()

org_file = sys.argv[2]
tmp_file = sys.argv[1]
print("Comparing: " + org_file + " with " + tmp_file)
print("")

orgCSV = []
with open(org_file, newline='') as orgFH:  # add encoding='utf-8' ?
    orgCSV = list(csv.reader(orgFH, delimiter='\t'))
if len(orgCSV) == 0:
    print("File " + org_file + " is empty!")
    sys.exit()

tmpCSV = []
with open(tmp_file, newline='') as tmpFH:  # add encoding='utf-8' ?
    tmpCSV = list(csv.reader(tmpFH, delimiter='\t'))
if len(tmpCSV) == 0:
    print("File " + org_file + " is empty!")
    sys.exit()

commonHeaders = []
missingOrg = []
lookupOrg = {}
missingTmp = []
lookupTmp = {}
for col in range(0, len(orgCSV[0])):
    if orgCSV[0][col] in tmpCSV[0]:
        commonHeaders.append(orgCSV[0][col])
        lookupOrg[orgCSV[0][col]] = col
    else:
        missingOrg.append(orgCSV[0][col])
for col in range(0, len(tmpCSV[0])):
    if tmpCSV[0][col] in commonHeaders:
        lookupTmp[tmpCSV[0][col]] = col
    else:
        missingTmp.append(tmpCSV[0][col])

ll = 1
if len(missingOrg) > 0:
    print("Columns missing in " + org_file + ":")
    for cell in missingOrg:
        print("   " + str(ll) + ": " + cell)
        ll += 1
if ll > 1:
    print("")
ll = 1
if len(missingTmp) > 0:
    print("Columns missing in " + tmp_file + ":")
    for cell in missingTmp:
        print("   " + str(ll) + ": " + cell)
        ll += 1
if ll > 1:
    print("")

if len(tmpCSV) != len(orgCSV):
    print("Files have different length!")
    sys.exit()

found = False
for row in range(1, len(orgCSV)):
    for cellName in commonHeaders:
        if orgCSV[row][lookupOrg[cellName]] != tmpCSV[row][lookupTmp[cellName]]:
            wellName = ":"
            if "well" in commonHeaders:
                wellName = " in well: " + orgCSV[row][lookupOrg["well"]]
            print ("Difference in row " + str(row + 1) + " in column " + cellName + wellName)
            print ("   File " + org_file + ": " + orgCSV[row][lookupOrg[cellName]])
            print ("   File " + tmp_file + ": " + tmpCSV[row][lookupTmp[cellName]])
            found = True
if not found:
    print("Files are identical!")
