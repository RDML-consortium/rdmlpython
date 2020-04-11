#!/bin/bash


echo "Run Test 1:"
python3 rdml.py -lrp test/test_1_raw_data.rdml \
                -e "Experiment 1" \
                -r "Run 1" \
                --pcrEfficiencyExl 0.05 \
                --excludeNoPlateau \
                --excludeEfficiency "mean" \
                --ignoreExclusion \
                --saveRaw test/temp_1_out_raw_data.tsv \
                --saveBaslineCorr test/temp_1_out_baseline_corrected_data.tsv \
                --saveResults test/temp_1_out_results.tsv \
                --timeRun

python3 test/diff_table.py test/temp_1_out_raw_data.tsv test/test_1_out_raw_data.tsv "Test 1 - raw data" 20 N Y
python3 test/diff_table.py test/temp_1_out_baseline_corrected_data.tsv test/test_1_out_baseline_corrected_data.tsv "Test 1 - baseline corrected data" 20 N Y
python3 test/diff_table.py test/temp_1_out_results.tsv test/test_1_out_results.tsv "Test 1 - results" 20 N Y


echo "Run Test 2:"
python3 rdml.py -lrp test/test_2_raw_data.rdml \
                -e "Experiment 1" \
                -r "Run 1" \
                --pcrEfficiencyExl 0.05 \
                --excludeNoPlateau \
                --excludeEfficiency "mean" \
                --ignoreExclusion \
                --saveRaw test/temp_2_out_raw_data.tsv \
                --saveBaslineCorr test/temp_2_out_baseline_corrected_data.tsv \
                --saveResults test/temp_2_out_results.tsv \
                --timeRun

python3 test/diff_table.py test/temp_2_out_raw_data.tsv test/test_2_out_raw_data.tsv "Test 2 - raw data" 20 N Y
python3 test/diff_table.py test/temp_2_out_baseline_corrected_data.tsv test/test_2_out_baseline_corrected_data.tsv "Test 2 - baseline corrected data" 20 N Y
python3 test/diff_table.py test/temp_2_out_results.tsv test/test_2_out_results.tsv "Test 2 - results" 20 N Y


echo "Run Test 3:"
python3 rdml.py -lrp test/test_3_raw_data.rdml \
                -e "Experiment 1" \
                -r "Run 1" \
                --pcrEfficiencyExl 0.05 \
                --excludeNoPlateau \
                --excludeEfficiency "mean" \
                --ignoreExclusion \
                --saveRaw test/temp_3_out_raw_data.tsv \
                --saveBaslineCorr test/temp_3_out_baseline_corrected_data.tsv \
                --saveResults test/temp_3_out_results.tsv \
                --timeRun

python3 test/diff_table.py test/temp_3_out_raw_data.tsv test/test_3_out_raw_data.tsv "Test 3 - raw data" 20 N Y
python3 test/diff_table.py test/temp_3_out_baseline_corrected_data.tsv test/test_3_out_baseline_corrected_data.tsv "Test 3 - baseline corrected data" 20 N Y
python3 test/diff_table.py test/temp_3_out_results.tsv test/test_3_out_results.tsv "Test 3 - results" 20 N Y


echo "Run Test 4:"
python3 rdml.py -lrp test/test_4_raw_data.rdml \
                -e "Experiment 1" \
                -r "Run 1" \
                --pcrEfficiencyExl 0.05 \
                --excludeNoPlateau \
                --excludeEfficiency "mean" \
                --ignoreExclusion \
                --saveRaw test/temp_4_out_raw_data.tsv \
                --saveBaslineCorr test/temp_4_out_baseline_corrected_data.tsv \
                --saveResults test/temp_4_out_results.tsv \
                --timeRun

python3 test/diff_table.py test/temp_4_out_raw_data.tsv test/test_4_out_raw_data.tsv "Test 4 - raw data" 20 N Y
python3 test/diff_table.py test/temp_4_out_baseline_corrected_data.tsv test/test_4_out_baseline_corrected_data.tsv "Test 4 - baseline corrected data" 20 N Y
python3 test/diff_table.py test/temp_4_out_results.tsv test/test_4_out_results.tsv "Test 4 - results" 20 N Y

echo "Run Test 5:"
python3 rdml.py -lrp test/test_5_raw_data.rdml \
                -e "Experiment 1" \
                -r "Run 1" \
                --pcrEfficiencyExl 0.05 \
                --excludeNoPlateau \
                --excludeEfficiency "mean" \
                --ignoreExclusion \
                --saveRaw test/temp_5_out_raw_data.tsv \
                --saveBaslineCorr test/temp_5_out_baseline_corrected_data.tsv \
                --saveResults test/temp_5_out_results.tsv \
                --timeRun

python3 test/diff_table.py test/temp_5_out_raw_data.tsv test/test_5_out_raw_data.tsv "Test 5 - raw data" 20 N Y
python3 test/diff_table.py test/temp_5_out_baseline_corrected_data.tsv test/test_5_out_baseline_corrected_data.tsv "Test 5 - baseline corrected data" 20 N Y
python3 test/diff_table.py test/temp_5_out_results.tsv test/test_5_out_results.tsv "Test 5 - results" 20 N Y



