#!usr/bin/perl
use warnings ;
use strict ;

my $reps = $ARGV[0] ;
my $simulator = "haploid_ind_seln" ;
my $simulator_dir = "/n/regal/hanage_lab/ForBrian/Haploid_fwdpp" ;
unless($reps){
	die "Error: didn't specify # replicates\n" ;
}

# SIMULATOR PARAMETERS
my $N = 10000 ;
my $theta_neut = 1000 ;
my $theta_seln = 0 ;
my $rho = 1000 ;
my $fragsize = 100000 ;
my $trLen = 300 ;
my $selcoeff = 0.0 ;
my $gens = 100000 ;
my $samplesize = 100 ;
my $replicates = 1 ;

# JOB SUBMISSION PARAMETERS
my $queue = "general" ;
my $time = 200 ;
my $mem = 1000 ;
my $n_cores = 1 ;

foreach my $rep ( 1..$reps ){
	system("mkdir rep${rep}") ;
	#system("cp ${simulator_dir}/${simulator} rep${rep}") ;

	my $script = "#!/bin/bash"."\n" ;
	$script = $script."#SBATCH -J fwdpp${rep}"."\n" ;
	$script = $script."#SBATCH -o out.fwdpp${rep}"."\n" ; 
	$script = $script."#SBATCH -e err.fwdpp${rep}"."\n" ;  
	$script = $script."#SBATCH -p ${queue}"."\n" ;
	$script = $script."#SBATCH -n ${n_cores}"."\n" ;
	$script = $script."#SBATCH -t ${time}"."\n" ;
	$script = $script."#SBATCH --mem ${mem}"."\n" ;
	$script = $script."#SBATCH -N 1"."\n\n" ;
	# move into directory
	$script = $script."cd rep${rep}\n" ;
	# load modules
	$script = $script."source new-modules.sh\n" ;
	$script = $script."module load gcc/7.1.0-fasrc01\n" ;
	# specify path to linked library libsequence
	$script = $script."export LD_LIBRARY_PATH=/n/home11/bjarnold/lib:\$LD_LIBRARY_PATH\n" ;
	# execute script
	$script = $script."${simulator_dir}/${simulator} " ;
	$script = $script."${N} " ;
	$script = $script."${theta_neut} " ;
	$script = $script."${theta_seln} " ;
	$script = $script."${rho} " ;
	$script = $script."${fragsize} " ;
	$script = $script."${trLen} " ;
	$script = $script."${selcoeff} " ;
	$script = $script."${gens} " ;
	$script = $script."${samplesize} " ;
	$script = $script."${replicates} " ;
	$script = $script."\$RANDOM\n" ;
	#print $script ;

	open OUT, ">fwdpp.sh" ; #creating text file to house script
	print OUT $script, "\n" ;
	close OUT ;
	system ("sbatch < fwdpp.sh") ;  #actually run script.
	system ("rm fwdpp.sh") ;   

}


exit ;
