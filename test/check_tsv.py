#!/usr/bin/python

import csv
import sys

if len(sys.argv) != 2:
    print("checkt_sv.py requires a file to check. Example use: ")
    print("   python3 checkt_sv.py results.tsv")
    sys.exit()

org_file = sys.argv[1]
print("Checking: " + org_file)
print("")

orgCSV = []
with open(org_file, newline='') as orgFH:  # add encoding='utf-8' ?
    orgCSV = list(csv.reader(orgFH, delimiter='\t'))
if len(orgCSV) == 0:
    print("File " + org_file + " is empty!")
    sys.exit()

rowCount = len(orgCSV[0])

print("The header has " + str(rowCount) + " columns.\n")

for row in range(0, len(orgCSV)):
    if rowCount != len(orgCSV[row]):
        print("  Row " + str(row + 1) + " has " + str(len(orgCSV[row])) + " columns!")
