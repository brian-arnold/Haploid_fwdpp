#include <iostream>
#include <gsl/gsl_randist.h>
#include <vector>
#include "fwdpp/fwdpp/sugar/GSLrng_t.hpp"

int main(){
	int x = 5 ;
	std::vector<double> fitnesses{1.,1.,1.,1.,1.} ;
	using GSLrng = fwdpp::GSLrng_t<fwdpp::GSL_RNG_MT19937>;

	gsl_ran_discrete_t *pdf ;
	pdf = gsl_ran_discrete_preproc(x, fitnesses.data());
	//gsl_ran_discrete_t * gsl_ran_discrete_preproc(x, fitnesses.data());

	size_t y ;
	GSLrng r(20845);

	x=gsl_ran_discrete(r ,pdf) ;
	std::cout << y << "\n" ;

	gsl_ran_discrete_free(pdf); 
	
	return 0 ;
}
