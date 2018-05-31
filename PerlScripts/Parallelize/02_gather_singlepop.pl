#!usr/bin/perl
use warnings ;
use strict ;

my $reps = $ARGV[0] ;
unless($reps){
	die "Error: didn't specify # replicates\n" ;
}

my $CpResults = "cat " ;
my $CpResults2 = "cat " ;
my $CpResultsBoth = "cat " ;
my $CpSelpos = "cat " ;
my $CpTime = "cat " ;

foreach my $rep ( 1..$reps ){
	$CpResults = $CpResults."rep${rep}/results.txt " ; 
	$CpResults2 = $CpResults2."rep${rep}/results2.txt " ; 
	$CpResultsBoth = $CpResultsBoth."rep${rep}/resultsboth.txt " ; 
	$CpSelpos = $CpSelpos."rep${rep}/selected_pos.txt " ; 
	$CpTime = $CpTime."rep${rep}/time.txt " ; 
}

$CpResults = $CpResults." > results.txt" ;
$CpResults2 = $CpResults2." > results2.txt" ;
$CpResultsBoth = $CpResultsBoth." > resultsboth.txt" ;
$CpSelpos = $CpSelpos." > selected_pos.txt" ;
$CpTime = $CpTime." > time.txt" ;

system("${CpResults}") ;
system("${CpResults2}") ;
system("${CpResultsBoth}") ;
system("${CpSelpos}") ;
system("${CpTime}") ;

system("rm -r rep*") ;
system("rm -r out.fwdpp*") ;
system("rm -r err.fwdpp*") ;


exit ;

