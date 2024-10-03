#!/bin/bash
rm nircam_*.txt
rm nirspec_*.txt
rm ns_layer_group.txt
rm nc_layer_group.txt
python3 parse_apt_pointing.py -i megasec.pointing --nirspec nirspec.txt --nircam nircam.txt --fits f200w_SCI.fits
