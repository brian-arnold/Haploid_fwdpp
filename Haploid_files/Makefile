.PHONY: all help clean

SHELL=/usr/bin/env bash -eo pipefail

F=../fwdpp
P=../haploid_ind_seln ../haploid_ind_neutral
FLAGS=-O3 -I$(F) -I/usr/local/include/ -std=c++11 -lsequence -lgsl -lgslcblas 

.SECONDARY:

.SUFFIXES:

all: $(P) #$(F)/Makefile

$(F)/Makefile:
	cd $(F) && ./configure


clean: #  $(F)/Makefile
	rm -f $(P)
	#$(MAKE) -C $(F) maintainer-clean

../%: %.cc $(F)/Makefile
	g++ $(FLAGS) -o $@ $<


