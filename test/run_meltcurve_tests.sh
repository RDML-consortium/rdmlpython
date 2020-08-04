#!/bin/bash

echo "Run MCA Test 1: exponential"
python3 rdml.py -mca test/test_mca_1_raw_data.rdml \
                -e "Experiment 1" \
                -r "Run 1" \
                --mcaNormMethod "exponential" \
                --mcaFluorSource "normalised" \
                --mcaPeakCutoff 0.0 \
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

echo "Run MCA Test 2: bilinear"
python3 rdml.py -mca test/test_mca_1_raw_data.rdml \
                -e "Experiment 1" \
                -r "Run 1" \
                --mcaNormMethod "bilinear" \
                --mcaFluorSource "normalised" \
                --mcaPeakCutoff 0.0 \
                --saveRaw test/temp_mca_2_out_raw_data.tsv \
                --saveDerivative test/temp_mca_2_out \
                --saveResults test/temp_mca_2_out_results.tsv \

python3 test/diff_table.py test/temp_mca_2_out_raw_data.tsv test/test_mca_1_out_raw_data.tsv "Test 2 - raw data" 20 N Y
python3 test/diff_table.py test/temp_mca_2_out_smoothed.tsv test/test_mca_1_out_smoothed.tsv "Test 2 - smoothed data" 20 N Y
python3 test/diff_table.py test/temp_mca_2_out_normalized.tsv test/test_mca_2_out_normalized.tsv "Test 2 - normalized data" 20 N Y
python3 test/diff_table.py test/temp_mca_2_out_firstDerivative.tsv test/test_mca_2_out_firstDerivative.tsv "Test 2 - first Derivative" 20 N Y
python3 test/diff_table.py test/temp_mca_2_out_secondDerivative.tsv test/test_mca_2_out_secondDerivative.tsv "Test 2 - second Derivative" 20 N Y
python3 test/diff_table.py test/temp_mca_2_out_thirdDerivative.tsv test/test_mca_2_out_thirdDerivative.tsv "Test 2 - third Derivative" 20 N Y
python3 test/diff_table.py test/temp_mca_2_out_results.tsv test/test_mca_2_out_results.tsv "Test 2 - results" 20 N Y

echo "Run MCA Test 3: combined"
python3 rdml.py -mca test/test_mca_1_raw_data.rdml \
                -e "Experiment 1" \
                -r "Run 1" \
                --mcaNormMethod "combined" \
                --mcaFluorSource "normalised" \
                --mcaPeakCutoff 0.0 \
                --saveRaw test/temp_mca_3_out_raw_data.tsv \
                --saveDerivative test/temp_mca_3_out \
                --saveResults test/temp_mca_3_out_results.tsv \

python3 test/diff_table.py test/temp_mca_3_out_raw_data.tsv test/test_mca_1_out_raw_data.tsv "Test 2 - raw data" 20 N Y
python3 test/diff_table.py test/temp_mca_3_out_smoothed.tsv test/test_mca_1_out_smoothed.tsv "Test 2 - smoothed data" 20 N Y
python3 test/diff_table.py test/temp_mca_3_out_normalized.tsv test/test_mca_3_out_normalized.tsv "Test 2 - normalized data" 20 N Y
python3 test/diff_table.py test/temp_mca_3_out_firstDerivative.tsv test/test_mca_3_out_firstDerivative.tsv "Test 2 - first Derivative" 20 N Y
python3 test/diff_table.py test/temp_mca_3_out_secondDerivative.tsv test/test_mca_3_out_secondDerivative.tsv "Test 2 - second Derivative" 20 N Y
python3 test/diff_table.py test/temp_mca_3_out_thirdDerivative.tsv test/test_mca_3_out_thirdDerivative.tsv "Test 2 - third Derivative" 20 N Y
python3 test/diff_table.py test/temp_mca_3_out_results.tsv test/test_mca_3_out_results.tsv "Test 2 - results" 20 N Y

echo "Run MCA Test 5: exponential"
python3 rdml.py -mca test/test_mca_5_raw_data.rdml \
                -e "Experiment_1" \
                -r "Run_1" \
                --mcaNormMethod "exponential" \
                --mcaFluorSource "normalised" \
                --mcaPeakCutoff 0.0 \
                --saveRaw test/temp_mca_5_out_raw_data.tsv \
                --saveDerivative test/temp_mca_5_out \
                --saveResults test/temp_mca_5_out_results.tsv \

python3 test/diff_table.py test/temp_mca_5_out_raw_data.tsv test/test_mca_5_out_raw_data.tsv "Test 5 - raw data" 20 N Y
python3 test/diff_table.py test/temp_mca_5_out_smoothed.tsv test/test_mca_5_out_smoothed.tsv "Test 5 - smoothed data" 20 N Y
python3 test/diff_table.py test/temp_mca_5_out_normalized.tsv test/test_mca_5_out_normalized.tsv "Test 5 - normalized data" 20 N Y
python3 test/diff_table.py test/temp_mca_5_out_firstDerivative.tsv test/test_mca_5_out_firstDerivative.tsv "Test 5 - first Derivative" 20 N Y
python3 test/diff_table.py test/temp_mca_5_out_secondDerivative.tsv test/test_mca_5_out_secondDerivative.tsv "Test 5 - second Derivative" 20 N Y
python3 test/diff_table.py test/temp_mca_5_out_thirdDerivative.tsv test/test_mca_5_out_thirdDerivative.tsv "Test 5 - third Derivative" 20 N Y
python3 test/diff_table.py test/temp_mca_5_out_results.tsv test/test_mca_5_out_results.tsv "Test 5 - results" 20 N Y

