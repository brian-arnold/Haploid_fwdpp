#!usr/bin/perl
use warnings ;
use strict ;

my $reps = $ARGV[0] ;



#####
## MAIN DATA STRUCTURES, HAVE EQUAL NUMBERS OF REPLICATES FOR SYN AND NONSYN
#####
my @AvgD_withinGene_SYN ;
my @AvgD_betweenGene_SYN ;
my @AvgD_withinGene_NONSYN ;
my @AvgD_betweenGene_NONSYN ;

my %AvgD_withinGene_byQuantile_SYN ;
my %AvgD_betweenGene_byQuantile_SYN ;
my %AvgD_withinGene_byQuantile_NONSYN ;
my %AvgD_betweenGene_byQuantile_NONSYN ;

my %AvgDprime_withinGene_byQuantile_SYN ;
my %AvgDprime_betweenGene_byQuantile_SYN ;
my %AvgDprime_withinGene_byQuantile_NONSYN ;
my %AvgDprime_betweenGene_byQuantile_NONSYN ;

#####
## TEMPORARY DATA STRUCTURES
#####
my %AvgD_withinGene_SYN_tmp ;
my %AvgD_betweenGene_SYN_tmp ;
my %AvgD_withinGene_NONSYN_tmp ;
my %AvgD_betweenGene_NONSYN_tmp ;

my %AvgD_withinGene_byQuantile_SYN_tmp ;
my %AvgD_betweenGene_byQuantile_SYN_tmp ;
my %AvgD_withinGene_byQuantile_NONSYN_tmp ;
my %AvgD_betweenGene_byQuantile_NONSYN_tmp ;

my %AvgDprime_withinGene_byQuantile_SYN_tmp ;
my %AvgDprime_betweenGene_byQuantile_SYN_tmp ;
my %AvgDprime_withinGene_byQuantile_NONSYN_tmp ;
my %AvgDprime_betweenGene_byQuantile_NONSYN_tmp ;

#####
## COLLECT INFORMATION, KEEPING TRACK OF REPLICATE SINCE SOME MAY NOT HAVE DATA
#####
foreach my $rep ( 1..$reps ){
	open IN, "<./rep${rep}_summaries/QC.txt" ;
	while(<IN>){
		chomp $_ ;
		if($_ =~ m/AVG LD SYN WITHIN GENE/){
			my @line = split(/:/, $_) ;
			$line[1] =~ s/\s//g ;
			$AvgD_withinGene_SYN_tmp{$rep} = $line[1] ;
		}
		if($_ =~ m/AVG LD SYN BETWEEN GENE/){
			my @line = split(/:/, $_) ;
			$line[1] =~ s/\s//g ;
			$AvgD_betweenGene_SYN_tmp{$rep} = $line[1] ;
		}	
		if($_ =~ m/AVG LD NONSYN WITHIN GENE/){
			my @line = split(/:/, $_) ;
			$line[1] =~ s/\s//g ;
			$AvgD_withinGene_NONSYN_tmp{$rep} = $line[1] ;
		}
		if($_ =~ m/AVG LD NONSYN BETWEEN GENE/){
			my @line = split(/:/, $_) ;
			$line[1] =~ s/\s//g ;
			$AvgD_betweenGene_NONSYN_tmp{$rep} = $line[1] ;
		}
	}
	close IN ;
	
	#### SYN
	open IN, "<./rep${rep}_summaries/D.Dprime.vs.Quant.SYN_WithinGene" ;
	while(<IN>){
		chomp $_ ;
		my @line = split(/\t/, $_) ;
		my $quant = $line[0] ;
		my $D = $line[2] ;
		my $Dprime = $line[3] ;
		$AvgD_withinGene_byQuantile_SYN_tmp{$quant}{$rep} = $D ;
		$AvgDprime_withinGene_byQuantile_SYN_tmp{$quant}{$rep} = $Dprime ;
	}
	close IN ;
	open IN, "<./rep${rep}_summaries/D.Dprime.vs.Quant.SYN_BetweenGene" ;
	while(<IN>){
		chomp $_ ;
		my @line = split(/\t/, $_) ;
		my $quant = $line[0] ;
		my $D = $line[2] ;
		my $Dprime = $line[3] ;
		$AvgD_betweenGene_byQuantile_SYN_tmp{$quant}{$rep} = $D ;
		$AvgDprime_betweenGene_byQuantile_SYN_tmp{$quant}{$rep} = $Dprime ;
	}
	close IN ;
	
	#### NONSYN
	open IN, "<./rep${rep}_summaries/D.Dprime.vs.Quant.NONSYN_WithinGene" ;
	while(<IN>){
		chomp $_ ;
		my @line = split(/\t/, $_) ;
		my $quant = $line[0] ;
		my $D = $line[2] ;
		my $Dprime = $line[3] ;
		$AvgD_withinGene_byQuantile_NONSYN_tmp{$quant}{$rep} = $D ;
		$AvgDprime_withinGene_byQuantile_NONSYN_tmp{$quant}{$rep} = $Dprime ;
	}
	close IN ;
	open IN, "<./rep${rep}_summaries/D.Dprime.vs.Quant.NONSYN_BetweenGene" ;
	while(<IN>){
		chomp $_ ;
		my @line = split(/\t/, $_) ;
		my $quant = $line[0] ;
		my $D = $line[2] ;
		my $Dprime = $line[3] ;
		$AvgD_betweenGene_byQuantile_NONSYN_tmp{$quant}{$rep} = $D ;
		$AvgDprime_betweenGene_byQuantile_NONSYN_tmp{$quant}{$rep} = $Dprime ;
	}
	close IN ;
}

###############
# NORMALIZE REPLICATES ACROSS SYN AND NONSYN CATEGORIES
###############
foreach my $rep ( sort{$a<=>$b} keys %AvgD_withinGene_SYN_tmp ){
	if(exists $AvgD_withinGene_NONSYN_tmp{$rep} ){
		if($AvgD_withinGene_SYN_tmp{$rep} ne "NA" && $AvgD_withinGene_NONSYN_tmp{$rep} ne "NA"){
			push @AvgD_withinGene_SYN, $AvgD_withinGene_SYN_tmp{$rep} ;
			push @AvgD_withinGene_NONSYN, $AvgD_withinGene_NONSYN_tmp{$rep} ;
		}
	}
}
foreach my $rep ( sort{$a<=>$b} keys %AvgD_betweenGene_SYN_tmp ){
	if(exists $AvgD_betweenGene_NONSYN_tmp{$rep} ){
		if($AvgD_betweenGene_SYN_tmp{$rep} ne "NA" && $AvgD_betweenGene_NONSYN_tmp{$rep} ne "NA"){
			push @AvgD_betweenGene_SYN, $AvgD_betweenGene_SYN_tmp{$rep} ;
			push @AvgD_betweenGene_NONSYN, $AvgD_betweenGene_NONSYN_tmp{$rep} ;
		}
	}
}

#D by quantile
foreach my $quant (sort{$a<=>$b} keys %AvgD_withinGene_byQuantile_SYN_tmp){
	foreach my $rep (sort{$a<=>$b} keys %{$AvgD_withinGene_byQuantile_SYN_tmp{$quant}}){
		if( exists $AvgD_withinGene_byQuantile_NONSYN_tmp{$quant}{$rep} ){
			push @{$AvgD_withinGene_byQuantile_SYN{$quant}}, $AvgD_withinGene_byQuantile_SYN_tmp{$quant}{$rep} ;
			push @{$AvgD_withinGene_byQuantile_NONSYN{$quant}}, $AvgD_withinGene_byQuantile_NONSYN_tmp{$quant}{$rep} ;
		}
	}
}
foreach my $quant (sort{$a<=>$b} keys %AvgD_betweenGene_byQuantile_SYN_tmp){
	foreach my $rep (sort{$a<=>$b} keys %{$AvgD_betweenGene_byQuantile_SYN_tmp{$quant}}){
		if( exists $AvgD_betweenGene_byQuantile_NONSYN_tmp{$quant}{$rep} ){
			push @{$AvgD_betweenGene_byQuantile_SYN{$quant}}, $AvgD_betweenGene_byQuantile_SYN_tmp{$quant}{$rep} ;
			push @{$AvgD_betweenGene_byQuantile_NONSYN{$quant}}, $AvgD_betweenGene_byQuantile_NONSYN_tmp{$quant}{$rep} ;
		}
	}
}

#Dprime by quantile
foreach my $quant (sort{$a<=>$b} keys %AvgDprime_withinGene_byQuantile_SYN_tmp){
	foreach my $rep (sort{$a<=>$b} keys %{$AvgDprime_withinGene_byQuantile_SYN_tmp{$quant}}){
		if( exists $AvgDprime_withinGene_byQuantile_NONSYN_tmp{$quant}{$rep} ){
			push @{$AvgDprime_withinGene_byQuantile_SYN{$quant}}, $AvgDprime_withinGene_byQuantile_SYN_tmp{$quant}{$rep} ;
			push @{$AvgDprime_withinGene_byQuantile_NONSYN{$quant}}, $AvgDprime_withinGene_byQuantile_NONSYN_tmp{$quant}{$rep} ;
		}
	}
}
foreach my $quant (sort{$a<=>$b} keys %AvgDprime_betweenGene_byQuantile_SYN_tmp){
	foreach my $rep (sort{$a<=>$b} keys %{$AvgDprime_betweenGene_byQuantile_SYN_tmp{$quant}}){
		if( exists $AvgDprime_betweenGene_byQuantile_NONSYN_tmp{$quant}{$rep} ){
			push @{$AvgDprime_betweenGene_byQuantile_SYN{$quant}}, $AvgDprime_betweenGene_byQuantile_SYN_tmp{$quant}{$rep} ;
			push @{$AvgDprime_betweenGene_byQuantile_NONSYN{$quant}}, $AvgDprime_betweenGene_byQuantile_NONSYN_tmp{$quant}{$rep} ;
		}
	}
}


###############
# PRINT OVERALL NUMBER OF OBSERVATIONS
###############
open OUT, ">./LinkageSummaries_NumObservations_SYN.txt" ;
print OUT "NumbObservations:\n" ;
print OUT "AvgD_withinGene_SYN\t", scalar @AvgD_withinGene_SYN, "\n" ;
print OUT "AvgD_betweenGene_SYN\t", scalar @AvgD_betweenGene_SYN, "\n" ;
foreach my $quant (sort{$a<=>$b} keys %AvgD_withinGene_byQuantile_SYN){
	print OUT "AvgD_withinGene_byQuantile_SYN_QUANT${quant}\t", scalar @{$AvgD_withinGene_byQuantile_SYN{$quant}}, "\n" ;
	print OUT "AvgDprime_withinGene_byQuantile_SYN${quant}\t", scalar @{$AvgDprime_withinGene_byQuantile_SYN{$quant}}, "\n" ;
	print OUT "AvgD_betweenGene_byQuantile_SYN${quant}\t", scalar @{$AvgD_betweenGene_byQuantile_SYN{$quant}}, "\n" ;
	print OUT "AvgDprime_betweenGene_byQuantile_SYN${quant}\t", scalar @{$AvgDprime_betweenGene_byQuantile_SYN{$quant}}, "\n" ;
}
print OUT "\n" ;
close OUT ;

open OUT, ">./LinkageSummaries_NumObservations_NONSYN.txt" ;
print OUT "NumbObservations:\n" ;
print OUT "AvgD_withinGene_NONSYN\t", scalar @AvgD_withinGene_NONSYN, "\n" ;
print OUT "AvgD_betweenGene_NONSYN\t", scalar @AvgD_betweenGene_NONSYN, "\n" ;
foreach my $quant (sort{$a<=>$b} keys %AvgD_withinGene_byQuantile_NONSYN){
	print OUT "AvgD_withinGene_byQuantile_NONSYN_QUANT${quant}\t", scalar @{$AvgD_withinGene_byQuantile_NONSYN{$quant}}, "\n" ;
	print OUT "AvgDprime_withinGene_byQuantile_NONSYN${quant}\t", scalar @{$AvgDprime_withinGene_byQuantile_NONSYN{$quant}}, "\n" ;
	print OUT "AvgD_betweenGene_byQuantile_NONSYN${quant}\t", scalar @{$AvgD_betweenGene_byQuantile_NONSYN{$quant}}, "\n" ;
	print OUT "AvgDprime_betweenGene_byQuantile_NONSYN${quant}\t", scalar @{$AvgDprime_betweenGene_byQuantile_NONSYN{$quant}}, "\n" ;
}
print OUT "\n" ;
close OUT ;
## CALCULATE DIFFS BETWEEN SYN AND NONSYN WITHIN EACH SIMULATION, THEN SORT DATA STRUCTS 
## TO GET QUANTILES

###############
# CALCULATE DIFFS BETWEEN SYN AND NONSYN
###############
my @AvgD_diff_withinGene ;
my @AvgD_diff_betweenGene ;
my %AvgD_diff_withinGene_byQuantile ;
my %AvgD_diff_betweenGene_byQuantile ;

foreach my $index (0..$#AvgD_withinGene_SYN){
	push @AvgD_diff_withinGene, $AvgD_withinGene_SYN[$index]-$AvgD_withinGene_NONSYN[$index] ;
}
foreach my $index (0..$#AvgD_betweenGene_SYN){
	push @AvgD_diff_betweenGene, $AvgD_betweenGene_SYN[$index]-$AvgD_betweenGene_NONSYN[$index] ;
}

foreach my $quant (sort{$a<=>$b} keys %AvgD_withinGene_byQuantile_SYN){
	foreach my $index (0..$#{$AvgD_withinGene_byQuantile_SYN{$quant}}){
		my $diff = ${$AvgD_withinGene_byQuantile_SYN{$quant}}[$index] - ${$AvgD_withinGene_byQuantile_NONSYN{$quant}}[$index] ;
		push @{$AvgD_diff_withinGene_byQuantile{$quant}}, $diff ;
	}
	foreach my $index (0..$#{$AvgD_betweenGene_byQuantile_SYN{$quant}}){
		my $diff = ${$AvgD_betweenGene_byQuantile_SYN{$quant}}[$index] - ${$AvgD_betweenGene_byQuantile_NONSYN{$quant}}[$index] ;
		push @{$AvgD_diff_betweenGene_byQuantile{$quant}}, $diff ;
	}	
}

###############
# SORT ARRAYS
###############
@AvgD_withinGene_SYN = sort{$a <=> $b} @AvgD_withinGene_SYN ;
@AvgD_betweenGene_SYN = sort{$a <=> $b} @AvgD_betweenGene_SYN ;
@AvgD_withinGene_NONSYN = sort{$a <=> $b} @AvgD_withinGene_NONSYN ;
@AvgD_betweenGene_NONSYN = sort{$a <=> $b} @AvgD_betweenGene_NONSYN ;

@AvgD_diff_withinGene = sort{$a <=> $b} @AvgD_diff_withinGene ;
@AvgD_diff_betweenGene = sort{$a <=> $b} @AvgD_diff_betweenGene ;

foreach my $quant (sort{$a<=>$b} keys %AvgD_withinGene_byQuantile_SYN){
	@{$AvgD_withinGene_byQuantile_SYN{$quant}} = sort{$a <=> $b} @{$AvgD_withinGene_byQuantile_SYN{$quant}} ;
	@{$AvgDprime_withinGene_byQuantile_SYN{$quant}} = sort{$a <=> $b} @{$AvgDprime_withinGene_byQuantile_SYN{$quant}} ;
	@{$AvgD_betweenGene_byQuantile_SYN{$quant}} = sort{$a <=> $b} @{$AvgD_betweenGene_byQuantile_SYN{$quant}} ;
	@{$AvgDprime_betweenGene_byQuantile_SYN{$quant}} = sort{$a <=> $b} @{$AvgDprime_betweenGene_byQuantile_SYN{$quant}} ;

	@{$AvgD_withinGene_byQuantile_NONSYN{$quant}} = sort{$a <=> $b} @{$AvgD_withinGene_byQuantile_NONSYN{$quant}} ;
	@{$AvgDprime_withinGene_byQuantile_NONSYN{$quant}} = sort{$a <=> $b} @{$AvgDprime_withinGene_byQuantile_NONSYN{$quant}} ;
	@{$AvgD_betweenGene_byQuantile_NONSYN{$quant}} = sort{$a <=> $b} @{$AvgD_betweenGene_byQuantile_NONSYN{$quant}} ;
	@{$AvgDprime_betweenGene_byQuantile_NONSYN{$quant}} = sort{$a <=> $b} @{$AvgDprime_betweenGene_byQuantile_NONSYN{$quant}} ;
	
	@{$AvgD_diff_withinGene_byQuantile{$quant}} = sort{$a <=> $b} @{$AvgD_diff_withinGene_byQuantile{$quant}} ;
	@{$AvgD_diff_betweenGene_byQuantile{$quant}} = sort{$a <=> $b} @{$AvgD_diff_betweenGene_byQuantile{$quant}} ;
}


###############
# PRINT MEAN LD WITHIN AND BETWEEN GENES
###############

my $lower90CI_tmp ;
my $median_tmp ;
my $upper90CI_tmp ;

### SYN
open OUT, ">./LinkageSummaries_MeanD_SYN.txt" ;
print OUT "SimsQuantile\tWithinGene\tbetweenGene\n" ;
print OUT ".05", "\t" ;
$lower90CI_tmp = int ((scalar @AvgD_withinGene_SYN)*0.05) - 1 ;
print OUT $AvgD_withinGene_SYN[$lower90CI_tmp], "\t" ;
$lower90CI_tmp = int ((scalar @AvgD_betweenGene_SYN)*0.05) - 1 ;
print OUT $AvgD_betweenGene_SYN[$lower90CI_tmp], "\n" ;

print OUT ".5", "\t" ;
$median_tmp = int ((scalar @AvgD_withinGene_SYN)*0.5) - 1 ;
print OUT $AvgD_withinGene_SYN[$median_tmp], "\t" ;
$median_tmp = int ((scalar @AvgD_betweenGene_SYN)*0.5) - 1 ;
print OUT $AvgD_betweenGene_SYN[$median_tmp], "\n" ;

print OUT ".95", "\t" ;
$upper90CI_tmp = int ((scalar @AvgD_withinGene_SYN)*0.95) - 1 ;
print OUT $AvgD_withinGene_SYN[$upper90CI_tmp], "\t" ;
$upper90CI_tmp = int ((scalar @AvgD_betweenGene_SYN)*0.95) - 1 ;
print OUT $AvgD_betweenGene_SYN[$upper90CI_tmp], "\n" ;
close OUT ;

### NONSYN
open OUT, ">./LinkageSummaries_MeanD_NONSYN.txt" ;
print OUT "SimsQuantile\tWithinGene\tbetweenGene\n" ;
print OUT ".05", "\t" ;
$lower90CI_tmp = int ((scalar @AvgD_withinGene_NONSYN)*0.05) - 1 ;
print OUT $AvgD_withinGene_NONSYN[$lower90CI_tmp], "\t" ;
$lower90CI_tmp = int ((scalar @AvgD_betweenGene_NONSYN)*0.05) - 1 ;
print OUT $AvgD_betweenGene_NONSYN[$lower90CI_tmp], "\n" ;

print OUT ".5", "\t" ;
$median_tmp = int ((scalar @AvgD_withinGene_NONSYN)*0.5) - 1 ;
print OUT $AvgD_withinGene_NONSYN[$median_tmp], "\t" ;
$median_tmp = int ((scalar @AvgD_betweenGene_NONSYN)*0.5) - 1 ;
print OUT $AvgD_betweenGene_NONSYN[$median_tmp], "\n" ;

print OUT ".95", "\t" ;
$upper90CI_tmp = int ((scalar @AvgD_withinGene_NONSYN)*0.95) - 1 ;
print OUT $AvgD_withinGene_NONSYN[$upper90CI_tmp], "\t" ;
$upper90CI_tmp = int ((scalar @AvgD_betweenGene_NONSYN)*0.95) - 1 ;
print OUT $AvgD_betweenGene_NONSYN[$upper90CI_tmp], "\n" ;
close OUT ;

### SYN-NONSYN
open OUT, ">./LinkageSummaries_MeanD_SYNminusNONSYN.txt" ;
print OUT "SimsQuantile\tWithinGene\tbetweenGene\n" ;
print OUT ".05", "\t" ;
$lower90CI_tmp = int ((scalar @AvgD_diff_withinGene)*0.05) - 1 ;
print OUT $AvgD_diff_withinGene[$lower90CI_tmp], "\t" ;
$lower90CI_tmp = int ((scalar @AvgD_diff_betweenGene)*0.05) - 1 ;
print OUT $AvgD_diff_betweenGene[$lower90CI_tmp], "\n" ;

print OUT ".5", "\t" ;
$median_tmp = int ((scalar @AvgD_diff_withinGene)*0.5) - 1 ;
print OUT $AvgD_diff_withinGene[$median_tmp], "\t" ;
$median_tmp = int ((scalar @AvgD_diff_betweenGene)*0.5) - 1 ;
print OUT $AvgD_diff_betweenGene[$median_tmp], "\n" ;

print OUT ".95", "\t" ;
$upper90CI_tmp = int ((scalar @AvgD_diff_withinGene)*0.95) - 1 ;
print OUT $AvgD_diff_withinGene[$upper90CI_tmp], "\t" ;
$upper90CI_tmp = int ((scalar @AvgD_diff_betweenGene)*0.95) - 1 ;
print OUT $AvgD_diff_betweenGene[$upper90CI_tmp], "\n" ;
close OUT ;




###############
# PRINT MEAN D BY QUANTILE
###############

### SYN
open OUT, ">./LinkageSummaries_DbyQuantile_SYN.txt" ;
print OUT "NONSYN_density_quant\tSims_quant\tWithinGene\tbetweenGene\n" ; 
foreach my $quant (sort{$a<=>$b} keys %AvgD_withinGene_byQuantile_SYN){
	print OUT $quant, "\t" ;
	print OUT ".05", "\t" ;
	my $lower90CI = int ((scalar @{$AvgD_withinGene_byQuantile_SYN{$quant}})*0.05) - 1 ;
	print OUT ${$AvgD_withinGene_byQuantile_SYN{$quant}}[$lower90CI], "\t" ;
	$lower90CI = int ((scalar @{$AvgD_betweenGene_byQuantile_SYN{$quant}})*0.05) - 1 ;
	print OUT ${$AvgD_betweenGene_byQuantile_SYN{$quant}}[$lower90CI], "\n" ;
	
	print OUT $quant, "\t" ;
	print OUT ".5", "\t" ;
	my $median = int ((scalar @{$AvgD_withinGene_byQuantile_SYN{$quant}})*0.5) - 1 ;
	print OUT ${$AvgD_withinGene_byQuantile_SYN{$quant}}[$median], "\t" ;
	$median = int ((scalar @{$AvgD_betweenGene_byQuantile_SYN{$quant}})*0.5) - 1 ;
	print OUT ${$AvgD_betweenGene_byQuantile_SYN{$quant}}[$median], "\n" ;
	
	print OUT $quant, "\t" ;
	print OUT ".95", "\t" ;	
	my $upper90CI = int ((scalar @{$AvgD_withinGene_byQuantile_SYN{$quant}})*0.95) - 1 ;
	print OUT ${$AvgD_withinGene_byQuantile_SYN{$quant}}[$upper90CI], "\t" ;
	$upper90CI = int ((scalar @{$AvgD_betweenGene_byQuantile_SYN{$quant}})*0.95) - 1 ;
	print OUT ${$AvgD_betweenGene_byQuantile_SYN{$quant}}[$upper90CI], "\n" ;

}
close OUT ;

### NONSYN
open OUT, ">./LinkageSummaries_DbyQuantile_NONSYN.txt" ;
print OUT "NONSYN_density_quant\tSims_quant\tWithinGene\tbetweenGene\n" ; 
foreach my $quant (sort{$a<=>$b} keys %AvgD_withinGene_byQuantile_NONSYN){
	print OUT $quant, "\t" ;
	print OUT ".05", "\t" ;
	my $lower90CI = int ((scalar @{$AvgD_withinGene_byQuantile_NONSYN{$quant}})*0.05) - 1 ;
	print OUT ${$AvgD_withinGene_byQuantile_NONSYN{$quant}}[$lower90CI], "\t" ;
	$lower90CI = int ((scalar @{$AvgD_betweenGene_byQuantile_NONSYN{$quant}})*0.05) - 1 ;
	print OUT ${$AvgD_betweenGene_byQuantile_NONSYN{$quant}}[$lower90CI], "\n" ;
	
	print OUT $quant, "\t" ;
	print OUT ".5", "\t" ;
	my $median = int ((scalar @{$AvgD_withinGene_byQuantile_NONSYN{$quant}})*0.5) - 1 ;
	print OUT ${$AvgD_withinGene_byQuantile_NONSYN{$quant}}[$median], "\t" ;
	$median = int ((scalar @{$AvgD_betweenGene_byQuantile_NONSYN{$quant}})*0.5) - 1 ;
	print OUT ${$AvgD_betweenGene_byQuantile_NONSYN{$quant}}[$median], "\n" ;
	
	print OUT $quant, "\t" ;
	print OUT ".95", "\t" ;	
	my $upper90CI = int ((scalar @{$AvgD_withinGene_byQuantile_NONSYN{$quant}})*0.95) - 1 ;
	print OUT ${$AvgD_withinGene_byQuantile_NONSYN{$quant}}[$upper90CI], "\t" ;
	$upper90CI = int ((scalar @{$AvgD_betweenGene_byQuantile_NONSYN{$quant}})*0.95) - 1 ;
	print OUT ${$AvgD_betweenGene_byQuantile_NONSYN{$quant}}[$upper90CI], "\n" ;

}
close OUT ;

### SYN - NONSYN
open OUT, ">./LinkageSummaries_DbyQuantile_SYNminusNONSYN.txt" ;
print OUT "NONSYN_density_quant\tSims_quant\tWithinGene\tbetweenGene\n" ; 
foreach my $quant (sort{$a<=>$b} keys %AvgD_diff_withinGene_byQuantile){
	print OUT $quant, "\t" ;
	print OUT ".05", "\t" ;
	my $lower90CI = int ((scalar @{$AvgD_diff_withinGene_byQuantile{$quant}})*0.05) - 1 ;
	print OUT ${$AvgD_diff_withinGene_byQuantile{$quant}}[$lower90CI], "\t" ;
	$lower90CI = int ((scalar @{$AvgD_diff_betweenGene_byQuantile{$quant}})*0.05) - 1 ;
	print OUT ${$AvgD_diff_betweenGene_byQuantile{$quant}}[$lower90CI], "\n" ;
	
	print OUT $quant, "\t" ;
	print OUT ".5", "\t" ;
	my $median = int ((scalar @{$AvgD_diff_withinGene_byQuantile{$quant}})*0.5) - 1 ;
	print OUT ${$AvgD_diff_withinGene_byQuantile{$quant}}[$median], "\t" ;
	$median = int ((scalar @{$AvgD_diff_betweenGene_byQuantile{$quant}})*0.5) - 1 ;
	print OUT ${$AvgD_diff_betweenGene_byQuantile{$quant}}[$median], "\n" ;
	
	print OUT $quant, "\t" ;
	print OUT ".95", "\t" ;	
	my $upper90CI = int ((scalar @{$AvgD_diff_withinGene_byQuantile{$quant}})*0.95) - 1 ;
	print OUT ${$AvgD_diff_withinGene_byQuantile{$quant}}[$upper90CI], "\t" ;
	$upper90CI = int ((scalar @{$AvgD_diff_betweenGene_byQuantile{$quant}})*0.95) - 1 ;
	print OUT ${$AvgD_diff_betweenGene_byQuantile{$quant}}[$upper90CI], "\n" ;

}
close OUT ;









###############
# PRINT MEAN Dprime BY QUANTILE
###############

# SYN
open OUT, ">./LinkageSummaries_DprimebyQuantile_SYN.txt" ;
print OUT "NONSYN_density_quant\tSims_quant\tWithinGene\tbetweenGene\n" ; 
foreach my $quant (sort{$a<=>$b} keys %AvgDprime_withinGene_byQuantile_SYN){
	print OUT $quant, "\t" ;
	print OUT ".05", "\t" ;
	my $lower90CI = int ((scalar @{$AvgDprime_withinGene_byQuantile_SYN{$quant}})*0.05) - 1 ;
	print OUT ${$AvgDprime_withinGene_byQuantile_SYN{$quant}}[$lower90CI], "\t" ;
	$lower90CI = int ((scalar @{$AvgDprime_betweenGene_byQuantile_SYN{$quant}})*0.05) - 1 ;
	print OUT ${$AvgDprime_betweenGene_byQuantile_SYN{$quant}}[$lower90CI], "\n" ;
	
	print OUT $quant, "\t" ;
	print OUT ".5", "\t" ;
	my $median = int ((scalar @{$AvgDprime_withinGene_byQuantile_SYN{$quant}})*0.5) - 1 ;
	print OUT ${$AvgDprime_withinGene_byQuantile_SYN{$quant}}[$median], "\t" ;
	$median = int ((scalar @{$AvgDprime_betweenGene_byQuantile_SYN{$quant}})*0.5) - 1 ;
	print OUT ${$AvgDprime_betweenGene_byQuantile_SYN{$quant}}[$median], "\n" ;
	
	print OUT $quant, "\t" ;
	print OUT ".95", "\t" ;	
	my $upper90CI = int ((scalar @{$AvgDprime_withinGene_byQuantile_SYN{$quant}})*0.95) - 1 ;
	print OUT ${$AvgDprime_withinGene_byQuantile_SYN{$quant}}[$upper90CI], "\t" ;
	$upper90CI = int ((scalar @{$AvgDprime_betweenGene_byQuantile_SYN{$quant}})*0.95) - 1 ;
	print OUT ${$AvgDprime_betweenGene_byQuantile_SYN{$quant}}[$upper90CI], "\n" ;

}
close OUT ;

# NONSYN
open OUT, ">./LinkageSummaries_DprimebyQuantile_NONSYN.txt" ;
print OUT "NONSYN_density_quant\tSims_quant\tWithinGene\tbetweenGene\n" ; 
foreach my $quant (sort{$a<=>$b} keys %AvgDprime_withinGene_byQuantile_NONSYN){
	print OUT $quant, "\t" ;
	print OUT ".05", "\t" ;
	my $lower90CI = int ((scalar @{$AvgDprime_withinGene_byQuantile_NONSYN{$quant}})*0.05) - 1 ;
	print OUT ${$AvgDprime_withinGene_byQuantile_NONSYN{$quant}}[$lower90CI], "\t" ;
	$lower90CI = int ((scalar @{$AvgDprime_betweenGene_byQuantile_NONSYN{$quant}})*0.05) - 1 ;
	print OUT ${$AvgDprime_betweenGene_byQuantile_NONSYN{$quant}}[$lower90CI], "\n" ;
	
	print OUT $quant, "\t" ;
	print OUT ".5", "\t" ;
	my $median = int ((scalar @{$AvgDprime_withinGene_byQuantile_NONSYN{$quant}})*0.5) - 1 ;
	print OUT ${$AvgDprime_withinGene_byQuantile_NONSYN{$quant}}[$median], "\t" ;
	$median = int ((scalar @{$AvgDprime_betweenGene_byQuantile_NONSYN{$quant}})*0.5) - 1 ;
	print OUT ${$AvgDprime_betweenGene_byQuantile_NONSYN{$quant}}[$median], "\n" ;
	
	print OUT $quant, "\t" ;
	print OUT ".95", "\t" ;	
	my $upper90CI = int ((scalar @{$AvgDprime_withinGene_byQuantile_NONSYN{$quant}})*0.95) - 1 ;
	print OUT ${$AvgDprime_withinGene_byQuantile_NONSYN{$quant}}[$upper90CI], "\t" ;
	$upper90CI = int ((scalar @{$AvgDprime_betweenGene_byQuantile_NONSYN{$quant}})*0.95) - 1 ;
	print OUT ${$AvgDprime_betweenGene_byQuantile_NONSYN{$quant}}[$upper90CI], "\n" ;

}
close OUT ;



exit ;


