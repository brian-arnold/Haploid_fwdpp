#!usr/bin/perl
use warnings ;
use strict ;

my $reps = $ARGV[0] ;

my %Results ; # @{$Results{$rep}{$quant}} = array ;
my %LowerQuants ;
my %UpperQuants ;

foreach my $rep (1 .. $reps){
	#@{$Results{$rep}} = () ;
	open IN, "<./rep${rep}_summaries/WithinGeneD_Bootreps.txt" ;
	while(<IN>){
		chomp $_ ;
		if($_ !~ m/Bootrep/){
			my @line = split(/\t/, $_) ;
			my $quant = $line[0] ;
			my $meanSynD = $line[2] ;
			my $meanNonsynD = $line[3] ;
			my $diff = $meanNonsynD - $meanSynD ;
			push @{$Results{$rep}{$quant}}, $diff ;
		}
	}
	close IN ;
}
foreach my $rep (1 .. $reps){
	foreach my $quant (sort {$a <=> $b} keys %{$Results{$rep}} ){
		my @tmp = @{$Results{$rep}{$quant}} ;
		@{$Results{$rep}{$quant}} = sort{$a <=> $b} @tmp ;
		$LowerQuants{$rep}{$quant} = ((scalar @{$Results{$rep}{$quant}})*0.05) - 1 ;
		$UpperQuants{$rep}{$quant} = ((scalar @{$Results{$rep}{$quant}})*0.95) - 1 ;
	}
}

open OUT, ">./WithinGeneD_Bootreps_Quantiles.txt" ;
print OUT "Rep\tQuant\tNumbObs\t0.05_Quant\t0.95_Quant\n" ;
my %NonsynGreaterThanSyn  ;
my %NonsynLessThanSyn  ;
my %NonsynEqualToSyn ;
my %Totals  ;
foreach my $rep (1 .. $reps){
	foreach my $quant (sort {$a <=> $b} keys %{$Results{$rep}} ){
		my $lowerQuantile = $LowerQuants{$rep}{$quant} ;
		my $upperQuantile = $UpperQuants{$rep}{$quant} ;
		print OUT $rep, "\t", $quant, "\t", scalar @{$Results{$rep}{$quant}}, "\t", ${$Results{$rep}{$quant}}[ $lowerQuantile ], "\t", ${$Results{$rep}{$quant}}[ $upperQuantile ], "\n" ;
		
		$Totals{$quant}++ ;
		if( ${$Results{$rep}{$quant}}[ $upperQuantile ] < 0 ){
			$NonsynLessThanSyn{$quant}++ ;
		}elsif( ${$Results{$rep}{$quant}}[ $lowerQuantile ] > 0){
			$NonsynGreaterThanSyn{$quant}++ ;
		}elsif( ${$Results{$rep}{$quant}}[ $upperQuantile ] >= 0 && ${$Results{$rep}{$quant}}[ $lowerQuantile ] <= 0 ){
			$NonsynEqualToSyn{$quant}++ ;
		}
	}
}
close OUT ;

open OUT, ">./WithinGeneD_SummarizeBootreps.txt" ;
print OUT "Quant", "\t", "NonsynLessThanSyn", "\t", "NonsynGreaterThanSyn", "\t", "NonsynEqualToSyn", "\n" ;
foreach my $quant (sort{$a <=> $b} keys %Totals){
	print OUT $quant, "\t", $NonsynLessThanSyn{$quant}/$Totals{$quant}, "\t", $NonsynGreaterThanSyn{$quant}/$Totals{$quant}, "\t", $NonsynEqualToSyn{$quant}/$Totals{$quant}, "\n" ;
}

close OUT ;

exit ;


