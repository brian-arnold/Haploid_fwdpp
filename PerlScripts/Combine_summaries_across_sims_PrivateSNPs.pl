#!usr/bin/perl
use warnings ;
use strict ;

my $reps = $ARGV[0] ;
my @PrivatePolyN1 ;
my @PrivatePolyN2 ;

#####
## COLLECT INFORMATION, KEEPING TRACK OF REPLICATE SINCE SOME MAY NOT HAVE DATA
#####
foreach my $rep ( 1..$reps ){
	open IN, "<./rep${rep}_summaries/PrivateSharedPoly.txt" ;
	while(<IN>){
		chomp $_ ;
		my @line = split(/\t/, $_) ;
		if($_ =~ m/PrivatePolyN1_Density/){
			push @PrivatePolyN1, $line[1] ;
		}
		if($_ =~ m/PrivatePolyN2_Density/){
			push @PrivatePolyN2, $line[1] ;
		}	
	}
	close IN ;
}

print "PrivateSNP_N1", "\t", "PrivateSNP_N2", "\n" ; 
foreach my $index (0..$#PrivatePolyN1){
	print $PrivatePolyN1[$index], "\t", $PrivatePolyN2[$index], "\n" ;
}

exit ;


