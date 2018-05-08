#!usr/bin/perl
use warnings ;
use strict ;

my $reps = $ARGV[0] ;
unless($reps){
	die "Error: didn't specify # replicates\n" ;
}

my $CpResults = "cat " ;
my $CpSelpos = "cat " ;
my $CpTime = "cat " ;

foreach my $rep ( 1..$reps ){
	$CpResults = $CpResults."rep${rep}/results.txt " ; 
	$CpSelpos = $CpSelpos."rep${rep}/selected_pos.txt " ; 
	$CpTime = $CpTime."rep${rep}/time.txt " ; 
}

$CpResults = $CpResults." > results.txt" ;
$CpSelpos = $CpSelpos." > selected_pos.txt" ;
$CpTime = $CpTime." > time.txt" ;

system("${CpResults}") ;
system("${CpSelpos}") ;
system("${CpTime}") ;

system("rm -r rep*") ;
system("rm -r out.fwdpp*") ;
system("rm -r err.fwdpp*") ;


exit ;
