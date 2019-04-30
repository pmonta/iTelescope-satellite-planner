all:
	./plan.py --start 2019-04-30T04:20:00 --telescope T24 --objects plan-objects.txt --tle objects.tle >plan.txt

download:
	./get_spice_kernels.sh
