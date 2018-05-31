#!usr/bin/perl
use warnings ;
use strict ;

my $reps = $ARGV[0] ;

######
# SIMULATOR NAME AND LOCATION
#####
my $simulator = "haploid_struct_seln" ;
my $simulator_dir = "/n/regal/hanage_lab/ForBrian/Haploid_fwdpp" ;
my $Libseq_dir = "/n/home11/bjarnold/lib" ;

######
# SIMULATOR PARAMETERS
#####
my $N1 = 5000 ;		# size of population 1
my $N2 = 5000 ; 	# size of population 2
my $m12 = 0.5 ;		# prob recomb frag in pop1 comes from pop2
my $m21 = 0.5 ;		# prob recomb frag in pop2 comes from pop1
my $theta_neut = 1000 ; # population mutation rate of neutral mutations (same for both pops)
my $theta_seln = 0.1 ;	# population mutation rate of selected mutations (same for both pops)
my $rho = 3000 ;	# population recombination rate (same for both pops)
my $fragsize = 100000 ;	# size of fragment to simulate
my $trLen = 1000 ;	# mean size of recombination tract lengths (geometrically distributed)
my $selcoeff = 0.005 ;	# selection coefficient of selected mutations
my $gens = 100000 ;	# number of generations to simulate (~10N gens for equilibrium)
my $samplesize = 200 ;	# sample size taken from the population
my $replicates = 1 ;	# number of replicates

######
# JOB SUBMISSION PARAMETERS
######
my $queue = "general" ; # name of partition
my $time = 1000 ;	# amount of time for job
my $mem = 1000 ;	# amount of memory allocated for job
my $n_cores = 1 ;	# number of cores

######
# JOB SUBMISSION LOOP
######
unless($reps){
	die "Error: didn't specify # replicates\n" ;
}
foreach my $rep ( 1..$reps ){
	system("mkdir rep${rep}") ;

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
	$script = $script."export LD_LIBRARY_PATH=${Libseq_dir}:\$LD_LIBRARY_PATH\n" ;
	# execute script
	$script = $script."${simulator_dir}/${simulator} " ;
	$script = $script."${N1} " ;
	$script = $script."${N2} " ;
	$script = $script."${m12} " ;
	$script = $script."${m21} " ;
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
	print $script ;

	open OUT, ">fwdpp.sh" ; #creating text file to house script
	print OUT $script, "\n" ;
	close OUT ;
	system ("sbatch < fwdpp.sh") ;  #actually run script.
	system ("rm fwdpp.sh") ;   

}


exit ;

