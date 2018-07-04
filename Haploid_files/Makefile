.PHONY: all help clean

SHELL=/usr/bin/env bash -eo pipefail

F=../fwdpp
P=../haploid_ind_seln ../haploid_ind_neutral ../haploid_struct_neutral ../haploid_struct_signseln ../haploid_struct_seln_delayedRecMig ../haploid_struct_seln_disruptedRecMig ../haploid_struct_seln_disruptedIndMig ../haploid_ind_seln_epi ../haploid_struct_seln_disruptedIndMig_nfds ../haploid_struct_seln_epi_delayedRecMig
FLAGS=-O3 -I$(F) -I/usr/local/include/ -std=c++11 -lsequence -lgsl -lgslcblas 
## for odyssey cluster
#FLAGS=-O3 -I$(F) -I/n/home11/bjarnold/include/ -L/n/home11/bjarnold/lib -std=c++11 -lsequence -lgsl -lgslcblas

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


