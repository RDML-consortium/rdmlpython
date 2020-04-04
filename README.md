# RDML-Python

A python library to handle RDML files.


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

Help on command line options:

`python3 rdml.py -h`

Validate RDML files:

`python3 rdml.py -v data.rdml`

Run LinRegPCR:

`python3 rdml.py -lrp data.rdml --pcrEfficiencyExl 0.05 --excludeNoPlateau --excludeEfficiency -o out_data.rdml --saveResults temp_1_out_results.tsv --timeRun`
