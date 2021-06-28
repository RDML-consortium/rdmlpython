# RDML-Python

A python library to handle RDML files.

Install with pip3
--------------------

`apt install python3 python3-pip`

`pip3 install rdmlpython`


Clone the repository
--------------------

`git clone --recursive https://github.com/RDML-consortium/rdmlpython.git rdmlpython`


Setup the environment
---------------------

Install the dependencies:

`sudo apt install python python-pip`

`pip install lxml numpy`


Command line functions
----------------------

The RDML-Python library offers a command line interface to use core functions. The 
command line interface is limited and should not replace the Python interface.

**Examples:**

Display a complete list of command line options:

`python3 rdml.py -h`

Validate RDML files:

`python3 rdml.py -v data.rdml`

List all experiments of a RDML file:

`python3 rdml.py -le data.rdml`

List all runs of the experiment 'exp_1' RDML file:

`python3 rdml.py -e exp_1 -lr data.rdml`

Run LinRegPCR:

`python3 rdml.py -lrp data.rdml --pcrEfficiencyExl 0.05 --excludeNoPlateau --excludeEfficiency "mean" -o out_data.rdml --saveResults temp_1_out_results.tsv --timeRun`

Run MeltCurveAnalysis:

`python3 rdml.py -mca data.rdml -e "Experiment 1" -r "Run 1" --mcaNormMethod "combined" --mcaFluorSource "normalised" --mcaPeakCutoff 0.0 --saveDerivative out_deriv --saveResults out_results.tsv`
