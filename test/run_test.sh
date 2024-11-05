#!/bin/bash


echo "Test: Fill Matrix Gaps:"
python3 rdml.py -lrp test/data_rdml_matrix_gaps.rdml \
                -e "Experiment 1" \
                -r "Run 1" \
                --pcrEfficiencyExl 0.05 \
                --excludeNoPlateau \
                --excludeEfficiency "outlier" \
                --ignoreExclusion \
                --saveRaw test/temp_out_raw_data_matrix_gaps.tsv \
                --saveBaslineCorr test/temp_out_baseline_corrected_data_matrix_gaps.tsv \
                --saveResults test/temp_out_results_matrix_gaps.tsv \
                --timeRun

python3 test/diff_table.py test/temp_out_raw_data_matrix_gaps.tsv test/out_raw_data_matrix_gaps.tsv "Test 1 - raw data" 20 N Y
python3 test/diff_table.py test/temp_out_baseline_corrected_data_matrix_gaps.tsv test/out_baseline_corrected_data_matrix_gaps.tsv "Test 1 - baseline corrected data" 20 N Y
python3 test/diff_table.py test/temp_out_results_matrix_gaps.tsv test/out_results_matrix_gaps.tsv "Test 1 - results" 20 N Y

