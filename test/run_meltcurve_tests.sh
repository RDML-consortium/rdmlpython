#!/bin/bash

echo "Run MCA Test 1:"
python3 rdml.py -mca test/test_mca_1_raw_data.rdml \
                -e "Experiment 1" \
                -r "Run 1" \
                --mcaNormMethod "exponential" \
                --mcaFluorSource "normalised" \
                --saveRaw test/temp_mca_1_out_raw_data.tsv \
                --saveDerivative test/temp_mca_1_out \
                --saveResults test/temp_mca_1_out_results.tsv \

python3 test/diff_table.py test/temp_mca_1_out_raw_data.tsv test/test_mca_1_out_raw_data.tsv "Test 1 - raw data" 20 N Y
python3 test/diff_table.py test/temp_mca_1_out_smoothed.tsv test/test_mca_1_out_smoothed.tsv "Test 1 - smoothed data" 20 N Y
python3 test/diff_table.py test/temp_mca_1_out_normalized.tsv test/test_mca_1_out_normalized.tsv "Test 1 - normalized data" 20 N Y
python3 test/diff_table.py test/temp_mca_1_out_firstDerivative.tsv test/test_mca_1_out_firstDerivative.tsv "Test 1 - first Derivative" 20 N Y
python3 test/diff_table.py test/temp_mca_1_out_secondDerivative.tsv test/test_mca_1_out_secondDerivative.tsv "Test 1 - second Derivative" 20 N Y
python3 test/diff_table.py test/temp_mca_1_out_thirdDerivative.tsv test/test_mca_1_out_thirdDerivative.tsv "Test 1 - third Derivative" 20 N Y
python3 test/diff_table.py test/temp_mca_1_out_results.tsv test/test_mca_1_out_results.tsv "Test 1 - results" 20 N Y

