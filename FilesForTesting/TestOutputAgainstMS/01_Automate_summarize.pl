#!usr/bin/perl
use warnings ;
use strict ;

foreach my $rep (1..100){
	#system("perl ../../PerlScripts/Summarize_ms_results.pl ms_output.txt ${rep}") ;
	#system("perl ../../PerlScripts/Summarize_ms_results.pl haploid_ind_output.txt ${rep}") ;
	system("perl ../../PerlScripts/Summarize_ms_results.pl haploid_struct_output.txt ${rep}") ;
}
