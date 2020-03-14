#!/bin/bash


echo "Run Test 1 chem hyd: "
python3 rdml.py -d test/test_1_chem_hyd.rdml

mv test/temp_out_raw_data.tsv test/temp_1_out_raw_data.tsv
mv test/temp_out_baseline_corrected_data.tsv test/temp_1_out_baseline_corrected_data.tsv
mv test/temp_out_results.tsv test/temp_1_chem_hyd_out_results.tsv

python3 test/diff_table.py test/temp_1_out_raw_data.tsv test/test_1_out_raw_data.tsv "Test 1 chem hyd - raw data" 20 N Y
python3 test/diff_table.py test/temp_1_out_baseline_corrected_data.tsv test/test_1_out_baseline_corrected_data.tsv "Test 1 chem hyd - baseline corrected data" 20 N Y
python3 test/diff_table.py test/temp_1_chem_hyd_out_results.tsv test/test_1_chem_hyd_out_results.tsv "Test 1 chem hyd - results" 20 N Y

echo "Run Test 1 chem dnazyme: "
python3 rdml.py -d test/test_1_chem_dnazyme.rdml

mv test/temp_out_raw_data.tsv test/temp_1_out_raw_data.tsv
mv test/temp_out_baseline_corrected_data.tsv test/temp_1_out_baseline_corrected_data.tsv
mv test/temp_out_results.tsv test/temp_1_chem_dnazyme_out_results.tsv

python3 test/diff_table.py test/temp_1_out_raw_data.tsv test/test_1_out_raw_data.tsv "Test 1 chem dnazyme - raw data" 20 N Y
python3 test/diff_table.py test/temp_1_out_baseline_corrected_data.tsv test/test_1_out_baseline_corrected_data.tsv "Test 1 chem dnazyme - baseline corrected data" 20 N Y
python3 test/diff_table.py test/temp_1_chem_dnazyme_out_results.tsv test/test_1_chem_dnazyme_out_results.tsv "Test 1 chem dnazyme - results" 20 N Y


echo "Run Test 2 chem labRP: "
python3 rdml.py -d test/test_2_chem_labRP.rdml

mv test/temp_out_raw_data.tsv test/temp_2_out_raw_data.tsv
mv test/temp_out_baseline_corrected_data.tsv test/temp_2_out_baseline_corrected_data.tsv
mv test/temp_out_results.tsv test/temp_2_chem_labRP_out_results.tsv

python3 test/diff_table.py test/temp_2_out_raw_data.tsv test/test_2_out_raw_data.tsv "Test 2 chem labRP - raw data" 20 N Y
python3 test/diff_table.py test/temp_2_out_baseline_corrected_data.tsv test/test_2_out_baseline_corrected_data.tsv "Test 2 chem labRP - baseline corrected data" 20 N Y
python3 test/diff_table.py test/temp_2_chem_labRP_out_results.tsv test/test_2_chem_labRP_out_results.tsv "Test 2 chem labRP - results" 20 N Y

echo "Run Test 2 chem dnazyme: "
python3 rdml.py -d test/test_2_chem_dnazyme.rdml

mv test/temp_out_raw_data.tsv test/temp_2_out_raw_data.tsv
mv test/temp_out_baseline_corrected_data.tsv test/temp_2_out_baseline_corrected_data.tsv
mv test/temp_out_results.tsv test/temp_2_chem_dnazyme_out_results.tsv

python3 test/diff_table.py test/temp_2_out_raw_data.tsv test/test_2_out_raw_data.tsv "Test 2 chem dnazyme - raw data" 20 N Y
python3 test/diff_table.py test/temp_2_out_baseline_corrected_data.tsv test/test_2_out_baseline_corrected_data.tsv "Test 2 chem dnazyme - baseline corrected data" 20 N Y
python3 test/diff_table.py test/temp_2_chem_dnazyme_out_results.tsv test/test_2_chem_dnazyme_out_results.tsv "Test 2 chem dnazyme - results" 20 N Y

echo "Run Test 2 chem dnazyme DS: "
python3 rdml.py -d test/test_2_chem_dnazymeDS.rdml

mv test/temp_out_raw_data.tsv test/temp_2_out_raw_data.tsv
mv test/temp_out_baseline_corrected_data.tsv test/temp_2_out_baseline_corrected_data.tsv
mv test/temp_out_results.tsv test/temp_2_chem_dnazymeDS_out_results.tsv

python3 test/diff_table.py test/temp_2_out_raw_data.tsv test/test_2_out_raw_data.tsv "Test 2 chem dnazyme DS - raw data" 20 N Y
python3 test/diff_table.py test/temp_2_out_baseline_corrected_data.tsv test/test_2_out_baseline_corrected_data.tsv "Test 2 chem dnazyme DS - baseline corrected data" 20 N Y
python3 test/diff_table.py test/temp_2_chem_dnazymeDS_out_results.tsv test/test_2_chem_dnazymeDS_out_results.tsv "Test 2 chem dnazyme DS - results" 20 N Y


echo "Run Test 3 chem hyd: "
python3 rdml.py -d test/test_3_chem_hyd.rdml

mv test/temp_out_raw_data.tsv test/temp_3_out_raw_data.tsv
mv test/temp_out_baseline_corrected_data.tsv test/temp_3_out_baseline_corrected_data.tsv
mv test/temp_out_results.tsv test/temp_3_chem_hyd_out_results.tsv

python3 test/diff_table.py test/temp_3_out_raw_data.tsv test/test_3_out_raw_data.tsv "Test 3 chem hyd - raw data" 20 N Y
python3 test/diff_table.py test/temp_3_out_baseline_corrected_data.tsv test/test_3_out_baseline_corrected_data.tsv "Test 3 chem hyd - baseline corrected data" 20 N Y
python3 test/diff_table.py test/temp_3_chem_hyd_out_results.tsv test/test_3_chem_hyd_out_results.tsv "Test 3 chem hyd - results" 20 N Y

echo "Run Test 3 chem labRP: "
python3 rdml.py -d test/test_3_chem_labRP.rdml

mv test/temp_out_raw_data.tsv test/temp_3_out_raw_data.tsv
mv test/temp_out_baseline_corrected_data.tsv test/temp_3_out_baseline_corrected_data.tsv
mv test/temp_out_results.tsv test/temp_3_chem_labRP_out_results.tsv

python3 test/diff_table.py test/temp_3_out_raw_data.tsv test/test_3_out_raw_data.tsv "Test 3 chem labRP - raw data" 20 N Y
python3 test/diff_table.py test/temp_3_out_baseline_corrected_data.tsv test/test_3_out_baseline_corrected_data.tsv "Test 3 chem labRP - baseline corrected data" 20 N Y
python3 test/diff_table.py test/temp_3_chem_labRP_out_results.tsv test/test_3_chem_labRP_out_results.tsv "Test 3 chem labRP - results" 20 N Y

echo "Run Test 3 chem dnazyme: "
python3 rdml.py -d test/test_3_chem_dnazyme.rdml

mv test/temp_out_raw_data.tsv test/temp_3_out_raw_data.tsv
mv test/temp_out_baseline_corrected_data.tsv test/temp_3_out_baseline_corrected_data.tsv
mv test/temp_out_results.tsv test/temp_3_chem_dnazyme_out_results.tsv

python3 test/diff_table.py test/temp_3_out_raw_data.tsv test/test_3_out_raw_data.tsv "Test 3 chem dnazyme - raw data" 20 N Y
python3 test/diff_table.py test/temp_3_out_baseline_corrected_data.tsv test/test_3_out_baseline_corrected_data.tsv "Test 3 chem dnazyme - baseline corrected data" 20 N Y
python3 test/diff_table.py test/temp_3_chem_dnazyme_out_results.tsv test/test_3_chem_dnazyme_out_results.tsv "Test 3 chem dnazyme - results" 20 N Y

echo "Run Test 3 chem dnazyme DS: "
python3 rdml.py -d test/test_3_chem_dnazymeDS.rdml

mv test/temp_out_raw_data.tsv test/temp_3_out_raw_data.tsv
mv test/temp_out_baseline_corrected_data.tsv test/temp_3_out_baseline_corrected_data.tsv
mv test/temp_out_results.tsv test/temp_3_chem_dnazymeDS_out_results.tsv

python3 test/diff_table.py test/temp_3_out_raw_data.tsv test/test_3_out_raw_data.tsv "Test 3 chem dnazyme DS - raw data" 20 N Y
python3 test/diff_table.py test/temp_3_out_baseline_corrected_data.tsv test/test_3_out_baseline_corrected_data.tsv "Test 3 chem dnazyme DS - baseline corrected data" 20 N Y
python3 test/diff_table.py test/temp_3_chem_dnazymeDS_out_results.tsv test/test_3_chem_dnazymeDS_out_results.tsv "Test 3 chem dnazyme DS - results" 20 N Y


echo "Run Test 4 chem labRP: "
python3 rdml.py -d test/test_4_chem_labRP.rdml

mv test/temp_out_raw_data.tsv test/temp_4_out_raw_data.tsv
mv test/temp_out_baseline_corrected_data.tsv test/temp_4_out_baseline_corrected_data.tsv
mv test/temp_out_results.tsv test/temp_4_chem_labRP_out_results.tsv

python3 test/diff_table.py test/temp_4_out_raw_data.tsv test/test_4_out_raw_data.tsv "Test 4 chem labRP - raw data" 20 N Y
python3 test/diff_table.py test/temp_4_out_baseline_corrected_data.tsv test/test_4_out_baseline_corrected_data.tsv "Test 4 chem labRP - baseline corrected data" 20 N Y
python3 test/diff_table.py test/temp_4_chem_labRP_out_results.tsv test/test_4_chem_labRP_out_results.tsv "Test 4 chem labRP - results" 20 N Y

echo "Run Test 4 chem dnazyme: "
python3 rdml.py -d test/test_4_chem_dnazyme.rdml

mv test/temp_out_raw_data.tsv test/temp_4_out_raw_data.tsv
mv test/temp_out_baseline_corrected_data.tsv test/temp_4_out_baseline_corrected_data.tsv
mv test/temp_out_results.tsv test/temp_4_chem_dnazyme_out_results.tsv

python3 test/diff_table.py test/temp_4_out_raw_data.tsv test/test_4_out_raw_data.tsv "Test 4 chem dnazyme - raw data" 20 N Y
python3 test/diff_table.py test/temp_4_out_baseline_corrected_data.tsv test/test_4_out_baseline_corrected_data.tsv "Test 4 chem dnazyme - baseline corrected data" 20 N Y
python3 test/diff_table.py test/temp_4_chem_dnazyme_out_results.tsv test/test_4_chem_dnazyme_out_results.tsv "Test 4 chem dnazyme - results" 20 N Y

