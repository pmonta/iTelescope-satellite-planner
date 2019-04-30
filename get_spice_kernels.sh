#!/bin/sh

mkdir spice-kernels
cd spice-kernels
wget https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls
wget https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de435.bsp
wget https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/earth_070425_370426_predict.bpc
