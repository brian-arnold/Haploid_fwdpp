#!usr/bin/perl
use warnings ;
use strict ;

my $outputFile_name = $ARGV[0] ;
my $reps = $ARGV[1] ;

# JOB SUBMISSION PARAMETERS
my $queue = "general" ;
my $time = 1000 ;
my $mem = 1000 ;
my $n_cores = 1 ;

foreach my $rep ( 5..100 ){
	#system("cp ${simulator_dir}/${simulator} rep${rep}") ;

	my $script = "#!/bin/bash"."\n" ;
	$script = $script."#SBATCH -J fwdpp${rep}"."\n" ;
	$script = $script."#SBATCH -o out.${rep}"."\n" ; 
	$script = $script."#SBATCH -e err.${rep}"."\n" ;  
	$script = $script."#SBATCH -p ${queue}"."\n" ;
	$script = $script."#SBATCH -n ${n_cores}"."\n" ;
	$script = $script."#SBATCH -t ${time}"."\n" ;
	$script = $script."#SBATCH --mem ${mem}"."\n" ;
	$script = $script."#SBATCH -N 1"."\n\n" ;
	
	# execute script
	$script = $script."perl /n/regal/hanage_lab/ForBrian/Haploid_fwdpp/PerlScripts/Summarize_ms_results_LDslices.pl ${outputFile_name} ${rep}" ;
	$script = $script."\n" ;
	print $script ;

	open OUT, ">fwdpp.sh" ; #creating text file to house script
	print OUT $script, "\n" ;
	close OUT ;
	system ("sbatch < fwdpp.sh") ;  #actually run script.
	system ("rm fwdpp.sh") ;   

}


exit ;
