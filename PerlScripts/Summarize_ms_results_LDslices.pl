#!usr/bin/perl
use warnings ;
use strict ;
use lib "/n/holyscratch/bomblies_lab/bjarnold/List-MoreUtils-0.33/lib" ;
use List::Util qw(first max maxstr min minstr reduce shuffle sum) ;

#############################################
# THIS SCRIPT LOOKS THROUGH A FILE OF MS-LIKE
# OUTPUT, ANALYZES A SINGLE REPLICATE
#############################################

my $file = $ARGV[0] ;
my $rep = $ARGV[1] ;

open QC, ">./QC.txt" ;

my $pop_size ;
my $Sample_size ;
my $num_reps ;
my $num_sites ;
my $rec_rate ;
my $theta ;
my $gene_conv_rate ;
my $TrLen ;
my $min_AF = 0.05 ;
my $max_AF = 0.95 ;

my $lower ;
my $upper ;
my $middle ;


my %DistSlices ;
@{$DistSlices{1}} = (1,50) ;
@{$DistSlices{2}} = (100,150) ;
@{$DistSlices{3}} = (300,350) ;
@{$DistSlices{4}} = (500,550) ;
@{$DistSlices{5}} = (5000,5050) ;


#############################################
# DATA STRUCTURES
#############################################
my $sum ;
my $sq_sum ;

my %multi_allele_sites ; #  $multi_allele_sites{num hits} = count ;
my %Seq_data ; #  $Seq_data{index} = scalar of 0's and 1's
my %Segsites ; # $Segsites{index} = site ;
my %Segsites_Bi_Unfiltered; ; # $Segsites{index} = site ;
my @Segsites_Bi_Unfiltered_indices ; # for mapping back to positin in string of 0's and 1's
my %Segsites_Bi_FilterAF ; # $Segsites_Bi_FilterAF{index} = site ;
my @Segsites_Bi_FilterAF_indices ; # for mapping back to positin in string of 0's and 1's
my %AFS ; # $AFS{1}{site} = AF ; 
my %AFS_forPrinting ;
my %Pairwise_Dprime  ; # $Pairwise_Dprime{site1}{site2} = D' ;

my %SumStat_Results ;
my %PairWise_LD_vs_dist ;
my %PairWise_PC_vs_dist ;
my %PairWise_LD_dist_slice ;
my %PairWise_PC_dist_slice ;

#############################################
# ASSEMBLE SEQUENCES PER GENE PER IND
#############################################
my $replicate = 0 ;
open IN, "<./${file}" ;
my $num_segsites_raw = 0 ;
while(<IN>){
	chomp $_ ;
	if($_ =~ m/haploid/ ){ # command used with arguments
		my @line = split(" ", $_) ;
		$pop_size = $line[1] ;
		$Sample_size = $line[7] ;
		$num_reps = $line[8] ;
		$theta = $line[2] ;
		$rec_rate = $line[3] ;
		$num_sites = $line[4] ;
		$TrLen = $line[5] ;
	}
	if($_ =~ m/haploid_ind_seln/ ){ # command used with arguments
		my @line = split(" ", $_) ;
		$Sample_size = $line[9] ;
		$num_reps = $line[10] ;
		$num_sites = $line[5] ;
	}
	if($_ =~ m/ms/ ){ # command used with arguments
		my @line = split(" ", $_) ;
		$Sample_size = $line[1] ;
		$num_reps = $line[2] ;
		$theta = $line[4] ;
		$rec_rate = $line[9] ;
		$num_sites = $line[7] ;
		$TrLen = $line[10] ;
	}
	if($_ =~ m/haploid_struct_neutral/ ){ # command used with arguments
		my @line = split(" ", $_) ;
		$Sample_size = $line[10] ;
		$num_reps = $line[11] ;
		$theta = $line[5] ;
		$rec_rate = $line[6] ;
		$num_sites = $line[7] ;
		$TrLen = $line[8] ;
	}
	if($_ =~ m/haploid_struct_seln/ ){ # command used with arguments
		my @line = split(" ", $_) ;
		$Sample_size = $line[12] ;
		$num_reps = $line[13] ;
		$theta = $line[5]+$line[6] ;
		$rec_rate = $line[7] ;
		$num_sites = $line[8] ;
		$TrLen = $line[9] ;
	}
	if($_ =~ m/\/\//){ # beginning of replicate
		$replicate++ ;
		if($replicate == $rep){
			my $segsitesline = <IN> ; chomp $segsitesline ;
			my @line = split(" ", $segsitesline) ;
			$num_segsites_raw = $line[1] ;
		
			my $positionsline = <IN> ; chomp $positionsline ;
			@line = split(" ", $positionsline) ;
			foreach my $x ( 1..$num_segsites_raw ){
				$Segsites{($x-1)} = (int($line[$x]*$num_sites)) ;
			}
		
			foreach my $ind (0..($Sample_size-1)){
				my $x = <IN> ; chomp $x ;
				$Seq_data{ $ind } = $x ;
			}
		}
	}

}
close IN ;
$SumStat_Results{"S_star"} = $num_segsites_raw/$num_sites ;

print QC "######## number of individuals with seqs\n" ;
print QC scalar keys %Seq_data, "\n" ;
print QC "num segsites ", scalar keys %Segsites, "\n" ;
foreach my $key (sort{$a<=>$b} keys %Segsites){
	print QC $key, "\t", $Segsites{$key}, "\n" ;
}
print QC "\n" ;


#############################################
### delete any seg site index which has same site, multiallelic
#############################################
foreach my $i ( 0..$num_segsites_raw-2 ){
	if(exists $Segsites{$i}){
		my $j = $i+1 ;
		if($Segsites{$j} == $Segsites{$i}){
			while($Segsites{$j} == $Segsites{$i}){
				$j++ ; #to capture potential runs of multiple hits at same site
				if(!exists $Segsites{$j}){ #reached end
					last ;
				}
			}
			## delete $i through $j-1 (-1 b.c. while loop added extra)
			foreach my $x ($i..($j-1)){
				delete $Segsites{$x} ;
			}
			$multi_allele_sites{($j-1-$i)+1} ++ ;
		}else{
			next ;
		}
	}
}
my $count = 0 ;
foreach my $key (sort{$a<=>$b} keys %Segsites){
    push @Segsites_Bi_Unfiltered_indices, $key ; # key is position of segsite in 0/1 string
    $Segsites_Bi_Unfiltered{$count} = $Segsites{$key} ;
    $count++ ;
}
undef %Segsites ; # ensure it's not used again
print QC "num segsites unfiltered ", scalar keys %Segsites_Bi_Unfiltered, "\n" ;
#foreach my $key (sort{$a<=>$b} keys %Segsites_Bi_Unfiltered){
#	print QC $key, "\t", $Segsites_Bi_Unfiltered_indices[$key], "\t", $Segsites_Bi_Unfiltered{$key}, "\n" ;
#}
#############################################
## CONSTRUCT SINGLE SITE %AFS FOR SUMSTATS AND LD
#############################################
foreach my $seg_site_index ( keys %Segsites_Bi_Unfiltered ){
    my $site = $Segsites_Bi_Unfiltered_indices[$seg_site_index] ;
	my $af = 0 ;
	foreach my $ind (keys %Seq_data){
		if( substr($Seq_data{$ind}, $site, 1) == 1 ){
			$af++ ;
		}
	}
	$AFS{1}{$Segsites_Bi_Unfiltered{$seg_site_index}} = $af ;
	$AFS_forPrinting{$af/$Sample_size}++ ;
}

#############################################
### Construct Segsites_Bi_FilterAF hash
### Continuously indexed biallelic sites
### Rare variants removed
#############################################
$count = 0 ;
foreach my $key (sort{$a<=>$b} keys %Segsites_Bi_Unfiltered){
	if( ($AFS{1}{$Segsites_Bi_Unfiltered{$key}}/$Sample_size) >= $min_AF && ($AFS{1}{$Segsites_Bi_Unfiltered{$key}}/$Sample_size) <= $max_AF ){
		$Segsites_Bi_FilterAF{$count} = $Segsites_Bi_Unfiltered{$key} ;
		$count++ ;
		push @Segsites_Bi_FilterAF_indices, $Segsites_Bi_Unfiltered_indices[$key] ;
	}
}
print QC "num segsites filtered: ", scalar keys %Segsites_Bi_FilterAF, "\n" ;
#foreach my $key (sort{$a<=>$b} keys %Segsites_Bi_FilterAF){
#	print QC $key, " ", $Segsites_Bi_FilterAF_indices[$key], " ", $Segsites_Bi_FilterAF{$key}, "\n" ;
#}


#############################################
## CONSTRUCT %Pairwise_LD_LongRange, %Pairwise_Compatibility_LongRange values used for all subsequent calculations
#############################################

OUTER: for ($lower = 0; $lower<((scalar keys %Segsites_Bi_FilterAF)-1); $lower++){
	my $site_lower_AF = $AFS{1}{$Segsites_Bi_FilterAF{$lower}}/$Sample_size ;
	INNER: for ($upper = ($lower+1); $upper<(scalar keys %Segsites_Bi_FilterAF); $upper++){
		my $dist = ($Segsites_Bi_FilterAF{$upper} - $Segsites_Bi_FilterAF{$lower}) ;
		my $site_upper_AF = $AFS{1}{$Segsites_Bi_FilterAF{$upper}}/$Sample_size ;
		# calc pairwise LD
		my $LD = 0 ;
		my $derived_haplo = 0 ;
		my %num_haplotypes ;
		my $compatibility ;
		#print QC "indices test filtered", "\t", $Segsites_Bi_FilterAF_indices[$lower], "\t", $Segsites_Bi_FilterAF_indices[$upper], "\n"  ;
		foreach my $ind (keys %Seq_data){
			if( (substr($Seq_data{$ind}, ($Segsites_Bi_FilterAF_indices[$upper]), 1)) == 0  &&  (substr($Seq_data{$ind}, ($Segsites_Bi_FilterAF_indices[$lower]), 1)) == 0 ){
				$num_haplotypes{1} ++ ;
			
			}
			if( (substr($Seq_data{$ind}, ($Segsites_Bi_FilterAF_indices[$upper]), 1)) == 0  &&  (substr($Seq_data{$ind}, ($Segsites_Bi_FilterAF_indices[$lower]), 1)) == 1 ){
				$num_haplotypes{2} ++ ;
			
			}
			if( (substr($Seq_data{$ind}, ($Segsites_Bi_FilterAF_indices[$upper]), 1)) == 1  &&  (substr($Seq_data{$ind}, ($Segsites_Bi_FilterAF_indices[$lower]), 1)) == 0 ){
				$num_haplotypes{3} ++ ;
			
			}
			if( (substr($Seq_data{$ind}, ($Segsites_Bi_FilterAF_indices[$upper]), 1)) == 1  &&  (substr($Seq_data{$ind}, ($Segsites_Bi_FilterAF_indices[$lower]), 1)) == 1 ){
				$derived_haplo++ ;
				$num_haplotypes{4} ++ ;
			}
		}
		if( (scalar keys %num_haplotypes) == 4 ){
			$compatibility = 0 ;
		}else{
			$compatibility = 1 ;
		}
		$LD = ($derived_haplo/$Sample_size) - ($site_lower_AF*$site_upper_AF) ;
		my $Dmax = 0 ;
		if ($LD>=0){
			if( $site_upper_AF*(1-$site_lower_AF) < $site_lower_AF*(1-$site_upper_AF) ){
				$Dmax = $site_upper_AF*(1-$site_lower_AF) ;
			}else{
				$Dmax = $site_lower_AF*(1-$site_upper_AF) ;
			}
		}elsif($LD<0){
			if( ($site_upper_AF*$site_lower_AF) < (1-$site_upper_AF)*(1-$site_lower_AF) ){
				$Dmax = $site_upper_AF*$site_lower_AF ;
			}else{
				$Dmax = (1-$site_upper_AF)*(1-$site_lower_AF) ;
			}
		}
		push @{$PairWise_PC_vs_dist{$dist}}, $compatibility ;
		push @{$PairWise_LD_vs_dist{$dist}}, abs($LD/$Dmax) ;
	}
}
my %MeanLDPerSlice ;
my %Totals ;
foreach my $dist (sort{$a<=>$b} keys %PairWise_LD_vs_dist){
	my $distSliceCat = 0 ;
	foreach my $distSlice (sort{$a<=>$b} keys %DistSlices){
		if($dist >= ${$DistSlices{$distSlice}}[0] && $dist <= ${$DistSlices{$distSlice}}[1]){
			$distSliceCat = $distSlice ;
		}
	}
	if($distSliceCat){
		foreach ( @{$PairWise_LD_vs_dist{$dist}} ){
			$MeanLDPerSlice{$distSliceCat} += $_ ;
			$Totals{$distSliceCat}++ ;
		}
	}
}


my $FiniteSites_Theta = 0 ;

# constants for wattersons estimators, see Tajima 1996, genetics, pp 1457-1465
my $a1 = 0 ;
foreach ( 1 .. ($Sample_size-1) ){ # sample size may vary depending on $max_miss_data
	$a1 += (1/$_) ;
}
my $a1_square = 0 ;
foreach ( 1 .. ($Sample_size-1) ){
	$a1_square += (1/$_**2) ;
}
my $a2 = (($a1**2) - $a1_square)/2 ;
my $c2 = (4*$a1/3) - (7*$a2/(3*$a1)) ;
$FiniteSites_Theta = $SumStat_Results{"S_star"}/($a1 - ($c2*$SumStat_Results{"S_star"})) ;


open OUT, ">./Summaries_rep${rep}.txt" ;
print OUT $rep, "\t", "FiniteSites_Theta", "\t", $FiniteSites_Theta, "\n" ;
foreach my $distSliceCat (sort{$a<=>$b} keys %MeanLDPerSlice){
	 print OUT "distSlice_", $distSliceCat, "\t", $MeanLDPerSlice{$distSliceCat}/$Totals{$distSliceCat}, "\n" ;
}
close OUT ;

exit ;

