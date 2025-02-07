#!/usr/bin/python

from __future__ import division
from __future__ import print_function

import os
import json
import sys
import time
import math

parent_dir = os.path.abspath(os.path.join(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir), os.pardir))
sys.path.append(parent_dir)

import rdmlpython as rdml

rdml_file = "data_vermeulen_raw.rdml"
out_file = "temp_vermeulen_resuls.csv"

# Time the test
startTime = time.time()
rd = rdml.Rdml(rdml_file)

expList = rd.experiments()
print("Experiments in file \"" + rdml_file + "\":")
if len(expList) < 1:
    print("No experiments found!")
    sys.exit(0)
ww = open(out_file, "w")
startLine = 0
for exp in expList:
    print(exp["id"] + ":")
    runList = exp.runs()
    if len(runList) < 1:
        print("No runs found!")
        sys.exit(0)
    for run in runList:
        print("  - " + run["id"])
        res = run.webAppLinRegPCR(pcrEfficiencyExl=0.05, updateTargetEfficiency=False, updateRDML=True, excludeNoPlateau=True, excludeEfficiency="outlier", excludeInstableBaseline=True)
        resTab = json.loads(res["LinRegPCR_Result_Table"])
        for tabRow in range(startLine, len(resTab)):
            if startLine == 0:
                ww.write("Experiment\tRun\t")
            else:
                ww.write(exp["id"] + "\t" + run["id"] + "\t")
            for tabCol in range(0, len(resTab[tabRow])):
                outCell = str(resTab[tabRow][tabCol]).replace("\t", ";")
                if tabCol < len(resTab[tabRow]) - 1:
                    ww.write(outCell + "\t")
                else:
                    ww.write(outCell + "\n")
            startLine = 1
ww.close()
endTime = time.time()
runTime = endTime - startTime
runMin = math.floor(runTime / 60.0)
runSec = runTime - runMin* 60.0
print("Runtime: " + str(runMin) + ":" + "{:.3f}".format(runSec))
sys.exit(0)
