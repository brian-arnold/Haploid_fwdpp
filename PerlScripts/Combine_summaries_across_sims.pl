#!usr/bin/perl
use warnings ;
use strict ;

my $reps = $ARGV[0] ;

my @AvgD_withinGene ;
my @AvgD_BetweenGene ;

my %AvgD_withinGene_byQuantile ;
my %AvgD_betweenGene_byQuantile ;

my %AvgDprime_withinGene_byQuantile ;
my %AvgDprime_betweenGene_byQuantile ;

foreach my $rep ( 1..$reps ){
	open IN, "<./rep${rep}_summaries/QC.txt" ;
	while(<IN>){
		chomp $_ ;
		if($_ =~ m/AVG LD SYN WITHIN GENE/){
			my @line = split(/:/, $_) ;
			$line[1] =~ s/\s//g ;
			push @AvgD_withinGene, $line[1] ;
		}
		if($_ =~ m/AVG LD SYN BETWEEN GENE/){
			my @line = split(/:/, $_) ;
			$line[1] =~ s/\s//g ;
			push @AvgD_BetweenGene, $line[1] ;
		}	
	}
	close IN ;
	
	open IN, "<./rep${rep}_summaries/D.Dprime.vs.Quant.SYN_WithinGene" ;
	while(<IN>){
		chomp $_ ;
		my @line = split(/\t/, $_) ;
		my $quant = $line[0] ;
		my $D = $line[2] ;
		my $Dprime = $line[3] ;
		push @{$AvgD_withinGene_byQuantile{$quant}}, $D ;
		push @{$AvgDprime_withinGene_byQuantile{$quant}}, $Dprime ;
	}
	close IN ;
	
	open IN, "<./rep${rep}_summaries/D.Dprime.vs.Quant.SYN_BetweenGene" ;
	while(<IN>){
		chomp $_ ;
		my @line = split(/\t/, $_) ;
		my $quant = $line[0] ;
		my $D = $line[2] ;
		my $Dprime = $line[3] ;
		push @{$AvgD_betweenGene_byQuantile{$quant}}, $D ;
		push @{$AvgDprime_betweenGene_byQuantile{$quant}}, $Dprime ;
	}
	close IN ;
}


# SORT ARRAYS
@AvgD_withinGene = sort{$a <=> $b} @AvgD_withinGene ;
@AvgD_BetweenGene = sort{$a <=> $b} @AvgD_BetweenGene ;

foreach my $quant (sort{$a<=>$b} keys %AvgD_withinGene_byQuantile){
	@{$AvgD_withinGene_byQuantile{$quant}} = sort{$a <=> $b} @{$AvgD_withinGene_byQuantile{$quant}} ;
	@{$AvgDprime_withinGene_byQuantile{$quant}} = sort{$a <=> $b} @{$AvgDprime_withinGene_byQuantile{$quant}} ;
	@{$AvgD_betweenGene_byQuantile{$quant}} = sort{$a <=> $b} @{$AvgD_betweenGene_byQuantile{$quant}} ;
	@{$AvgDprime_betweenGene_byQuantile{$quant}} = sort{$a <=> $b} @{$AvgDprime_betweenGene_byQuantile{$quant}} ;
}

my $quant1 = int(0.25*$reps - 1) ;
my $quant2 = int(0.5*$reps - 1) ;
my $quant3 = int(0.75*$reps - 1) ;

print "AvgLD_withinGene quartile1: ", $AvgD_withinGene[$quant1], "\n" ;
print "AvgLD_withinGene median: ", $AvgD_withinGene[$quant2], "\n" ;
print "AvgLD_withinGene quartile3: ", $AvgD_withinGene[$quant3], "\n" ;
print "\n" ;
print "AvgLD_BetweenGene quartile1: ", $AvgD_BetweenGene[$quant1], "\n" ;
print "AvgLD_BetweenGene median: ", $AvgD_BetweenGene[$quant2], "\n" ;
print "AvgLD_BetweenGene quartile3: ", $AvgD_BetweenGene[$quant3], "\n" ;
print "\n" ;

foreach my $quant (sort{$a<=>$b} keys %AvgD_withinGene_byQuantile){
	print "AvgD_withinGene_byQuantile QUANT${quant} quartile1: ", ${$AvgD_withinGene_byQuantile{$quant}}[$quant1], "\n" ;
	print "AvgD_withinGene_byQuantile QUANT${quant} median: ", ${$AvgD_withinGene_byQuantile{$quant}}[$quant2], "\n" ;
	print "AvgD_withinGene_byQuantile QUANT${quant} quartile3: ", ${$AvgD_withinGene_byQuantile{$quant}}[$quant3], "\n" ;
	print "#######\n" ;
}
print "\n\n" ;

foreach my $quant (sort{$a<=>$b} keys %AvgD_betweenGene_byQuantile){
	print "AvgD_betweenGene_byQuantile QUANT${quant} quartile1: ", ${$AvgD_betweenGene_byQuantile{$quant}}[$quant1], "\n" ;
	print "AvgD_betweenGene_byQuantile QUANT${quant} median: ", ${$AvgD_betweenGene_byQuantile{$quant}}[$quant2], "\n" ;
	print "AvgD_betweenGene_byQuantile QUANT${quant} quartile3: ", ${$AvgD_betweenGene_byQuantile{$quant}}[$quant3], "\n" ;
	print "#######\n" ;
}
print "\n\n" ;

foreach my $quant (sort{$a<=>$b} keys %AvgDprime_withinGene_byQuantile){
	print "AvgDprime_withinGene_byQuantile QUANT${quant} quartile1: ", ${$AvgDprime_withinGene_byQuantile{$quant}}[$quant1], "\n" ;
	print "AvgDprime_withinGene_byQuantile QUANT${quant} median: ", ${$AvgDprime_withinGene_byQuantile{$quant}}[$quant2], "\n" ;
	print "AvgDprime_withinGene_byQuantile QUANT${quant} quartile3: ", ${$AvgDprime_withinGene_byQuantile{$quant}}[$quant3], "\n" ;
	print "#######\n" ;
}
print "\n\n" ;

foreach my $quant (sort{$a<=>$b} keys %AvgDprime_betweenGene_byQuantile){
	print "AvgDprime_betweenGene_byQuantile QUANT${quant} quartile1: ", ${$AvgDprime_betweenGene_byQuantile{$quant}}[$quant1], "\n" ;
	print "AvgDprime_betweenGene_byQuantile QUANT${quant} median: ", ${$AvgDprime_betweenGene_byQuantile{$quant}}[$quant2], "\n" ;
	print "AvgDprime_betweenGene_byQuantile QUANT${quant} quartile3: ", ${$AvgDprime_betweenGene_byQuantile{$quant}}[$quant3], "\n" ;
	print "#######\n" ;
}


exit ;


