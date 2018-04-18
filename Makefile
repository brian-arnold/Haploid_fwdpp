.PHONY: all help clean

SHELL=/usr/bin/env bash -eo pipefail

F=fwdpp
P=haploid_ind

.SECONDARY:

.SUFFIXES:

all: $(P) #$(F)/Makefile

$(P): $(P).cc #$(F)/Makefile
	g++ -I. -I$(F) -I/usr/local/include/ -std=c++11 -lgsl -lgslcblas -o $@ $<


$(F)/Makefile:
	cd $(F) && ./configure


clean: $(F)/Makefile
	rm -f $(P)
	#$(MAKE) -C $(F) maintainer-clean

