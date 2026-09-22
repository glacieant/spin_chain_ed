#!/bin/bash

##################################################
# random heisenberg spin chain                   #
# exact diagonalisation                          #
##################################################

# constants in constants.txt
# data in data
# figures in figs

set -e
cd "$(dirname "$0")"
mkdir -p data figs

source constants.txt

# virtual environment
python3 -m venv ENV
source ENV/bin/activate

pip install --upgrade pip
pip install -r requirements.txt

make

# working data directory
cd data
> randoutput.txt
KEY=$([ "$DISORDER" = field ] && echo FIELD_DELTA || echo DELTA)
for d in ${!KEY}; do
	sed "s/^$KEY=.*/$KEY=$d/" ../constants.txt > constants.txt
	../randspin
done
cp ../constants.txt .

# figures
python ../randplot.py "$DISORDER" "$YLABEL"
mv heisenberg_plot_power_distro.pdf ../figs/
