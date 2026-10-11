.PHONY: all test test-full bao cmb cmb-full sne sne-full data clean

all: bao cmb

test:        ## fast self-checks + regression vs canonical numbers
	python3 tests/run_checks.py

test-full:   ## full 8-angle CMB scan regression (minutes)
	python3 tests/run_checks.py --full

bao:         ## BAO: DESI DR2 vs bimetric baseline (full covariance)
	python3 -m probes.bao --cov full

bao-table4:  ## BAO with Table IV block-diagonal covariance (0.195)
	python3 -m probes.bao --cov table4

cmb:         ## CMB theta-scan, 8 angles (bimetric E^2)
	python3 -m probes.cmb_theta_scan

data:        ## download Pantheon+ data (needed for SNe)
	bash data/download_data.sh

sne:         ## SNe diagonal covariance (needs make data)
	python3 -m probes.sne

sne-full:    ## SNe full STAT+SYS covariance (needs make data)
	python3 -m probes.sne --full

clean:
	rm -f sne_results_v5.txt
