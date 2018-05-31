#!usr/bin/perl
use warnings ;
use strict ;

my $reps = $ARGV[0] ;

# JOB SUBMISSION PARAMETERS
my $queue = "serial_requeue" ;
my $time = 1400 ;
my $mem = 3700 ;
my $n_cores = 1 ;

foreach my $rep ( 1..$reps ){
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
	#$script = $script."perl /n/regal/hanage_lab/ForBrian/Haploid_fwdpp/PerlScripts/LD_neut_sel_perQuantile.pl results.txt selected_pos.txt ${rep}" ;
	$script = $script."perl /n/regal/hanage_lab/ForBrian/Haploid_fwdpp/PerlScripts/LD_between_neut_sel_bootstrap.pl results.txt selected_pos.txt 214 ${rep} 1000 5" ;
	$script = $script."\n" ;
	print $script ;

	open OUT, ">fwdpp.sh" ; #creating text file to house script
	print OUT $script, "\n" ;
	close OUT ;
	system ("sbatch < fwdpp.sh") ;  #actually run script.
	system ("rm fwdpp.sh") ;   

}


exit ;

