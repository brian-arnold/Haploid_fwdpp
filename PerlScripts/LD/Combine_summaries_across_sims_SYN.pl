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
	if(-e "./rep${rep}_summaries/QC.txt"){
		open IN, "<./rep${rep}_summaries/QC.txt" ;
	}elsif(-e "./rep${rep}_summaries/QC_SharedSNPs.txt"){
		open IN, "<./rep${rep}_summaries/QC_SharedSNPs.txt" ;
	}else{
		print "Couldn't open QC file\n" ;
		exit ;
	}
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


###############
# PRINT OVERALL NUMBER OF OBSERVATIONS
###############
open OUT, ">./LinkageSummaries_NumObservations_SYN.txt" ;
print OUT "NumbObservations:\n" ;
print OUT "AvgD_withinGene\t", scalar @AvgD_withinGene, "\n" ;
print OUT "AvgD_BetweenGene\t", scalar @AvgD_BetweenGene, "\n" ;
foreach my $quant (sort{$a<=>$b} keys %AvgD_withinGene_byQuantile){
	@{$AvgD_withinGene_byQuantile{$quant}} = sort{$a <=> $b} @{$AvgD_withinGene_byQuantile{$quant}} ;
	@{$AvgDprime_withinGene_byQuantile{$quant}} = sort{$a <=> $b} @{$AvgDprime_withinGene_byQuantile{$quant}} ;
	@{$AvgD_betweenGene_byQuantile{$quant}} = sort{$a <=> $b} @{$AvgD_betweenGene_byQuantile{$quant}} ;
	@{$AvgDprime_betweenGene_byQuantile{$quant}} = sort{$a <=> $b} @{$AvgDprime_betweenGene_byQuantile{$quant}} ;
	
	print OUT "AvgD_withinGene_byQuantile_QUANT${quant}\t", scalar @{$AvgD_withinGene_byQuantile{$quant}}, "\n" ;
	print OUT "AvgDprime_withinGene_byQuantile${quant}\t", scalar @{$AvgDprime_withinGene_byQuantile{$quant}}, "\n" ;
	print OUT "AvgD_betweenGene_byQuantile${quant}\t", scalar @{$AvgD_betweenGene_byQuantile{$quant}}, "\n" ;
	print OUT "AvgDprime_betweenGene_byQuantile${quant}\t", scalar @{$AvgDprime_betweenGene_byQuantile{$quant}}, "\n" ;

}
print OUT "\n" ;
close OUT ;


###############
# PRINT MEAN LD WITHIN AND BETWEEN GENES
###############
open OUT, ">./LinkageSummaries_MeanD_SYN.txt" ;
print OUT "SimsQuantile\tWithinGene\tBetweenGene\n" ;
print OUT ".05", "\t" ;
my $lower90CI_tmp = int ((scalar @AvgD_withinGene)*0.05) - 1 ;
print OUT $AvgD_withinGene[$lower90CI_tmp], "\t" ;
$lower90CI_tmp = int ((scalar @AvgD_BetweenGene)*0.05) - 1 ;
print OUT $AvgD_BetweenGene[$lower90CI_tmp], "\n" ;

print OUT ".5", "\t" ;
my $median_tmp = int ((scalar @AvgD_withinGene)*0.5) - 1 ;
print OUT $AvgD_withinGene[$median_tmp], "\t" ;
$median_tmp = int ((scalar @AvgD_BetweenGene)*0.5) - 1 ;
print OUT $AvgD_BetweenGene[$median_tmp], "\n" ;

print OUT ".95", "\t" ;
my $upper90CI_tmp = int ((scalar @AvgD_withinGene)*0.95) - 1 ;
print OUT $AvgD_withinGene[$upper90CI_tmp], "\t" ;
$upper90CI_tmp = int ((scalar @AvgD_BetweenGene)*0.95) - 1 ;
print OUT $AvgD_BetweenGene[$upper90CI_tmp], "\n" ;
close OUT ;


open OUT, ">./LinkageSummaries_DbyQuantile_SYN.txt" ;
print OUT "NONSYN_density_quant\tSims_quant\tWithinGene\tBetweenGene\n" ; 
foreach my $quant (sort{$a<=>$b} keys %AvgD_withinGene_byQuantile){
	print OUT $quant, "\t" ;
	print OUT ".05", "\t" ;
	my $lower90CI = int ((scalar @{$AvgD_withinGene_byQuantile{$quant}})*0.05) - 1 ;
	print OUT ${$AvgD_withinGene_byQuantile{$quant}}[$lower90CI], "\t" ;
	$lower90CI = int ((scalar @{$AvgD_betweenGene_byQuantile{$quant}})*0.05) - 1 ;
	print OUT ${$AvgD_betweenGene_byQuantile{$quant}}[$lower90CI], "\n" ;
	
	print OUT $quant, "\t" ;
	print OUT ".5", "\t" ;
	my $median = int ((scalar @{$AvgD_withinGene_byQuantile{$quant}})*0.5) - 1 ;
	print OUT ${$AvgD_withinGene_byQuantile{$quant}}[$median], "\t" ;
	$median = int ((scalar @{$AvgD_betweenGene_byQuantile{$quant}})*0.5) - 1 ;
	print OUT ${$AvgD_betweenGene_byQuantile{$quant}}[$median], "\n" ;
	
	print OUT $quant, "\t" ;
	print OUT ".95", "\t" ;	
	my $upper90CI = int ((scalar @{$AvgD_withinGene_byQuantile{$quant}})*0.95) - 1 ;
	print OUT ${$AvgD_withinGene_byQuantile{$quant}}[$upper90CI], "\t" ;
	$upper90CI = int ((scalar @{$AvgD_betweenGene_byQuantile{$quant}})*0.95) - 1 ;
	print OUT ${$AvgD_betweenGene_byQuantile{$quant}}[$upper90CI], "\n" ;

}
close OUT ;

open OUT, ">./LinkageSummaries_DprimebyQuantile_SYN.txt" ;
print OUT "NONSYN_density_quant\tSims_quant\tWithinGene\tBetweenGene\n" ; 
foreach my $quant (sort{$a<=>$b} keys %AvgDprime_withinGene_byQuantile){
	print OUT $quant, "\t" ;
	print OUT ".05", "\t" ;
	my $lower90CI = int ((scalar @{$AvgDprime_withinGene_byQuantile{$quant}})*0.05) - 1 ;
	print OUT ${$AvgDprime_withinGene_byQuantile{$quant}}[$lower90CI], "\t" ;
	$lower90CI = int ((scalar @{$AvgDprime_betweenGene_byQuantile{$quant}})*0.05) - 1 ;
	print OUT ${$AvgDprime_betweenGene_byQuantile{$quant}}[$lower90CI], "\n" ;
	
	print OUT $quant, "\t" ;
	print OUT ".5", "\t" ;
	my $median = int ((scalar @{$AvgDprime_withinGene_byQuantile{$quant}})*0.5) - 1 ;
	print OUT ${$AvgDprime_withinGene_byQuantile{$quant}}[$median], "\t" ;
	$median = int ((scalar @{$AvgDprime_betweenGene_byQuantile{$quant}})*0.5) - 1 ;
	print OUT ${$AvgDprime_betweenGene_byQuantile{$quant}}[$median], "\n" ;
	
	print OUT $quant, "\t" ;
	print OUT ".95", "\t" ;	
	my $upper90CI = int ((scalar @{$AvgDprime_withinGene_byQuantile{$quant}})*0.95) - 1 ;
	print OUT ${$AvgDprime_withinGene_byQuantile{$quant}}[$upper90CI], "\t" ;
	$upper90CI = int ((scalar @{$AvgDprime_betweenGene_byQuantile{$quant}})*0.95) - 1 ;
	print OUT ${$AvgDprime_betweenGene_byQuantile{$quant}}[$upper90CI], "\n" ;

}
close OUT ;


exit ;
