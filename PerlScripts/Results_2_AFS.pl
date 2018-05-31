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

my $pop_size ;
my $Sample_size ;
my $num_reps ;
my $num_sites ;
my $rec_rate ;
my $theta ;
my $gene_conv_rate ;
my $TrLen ;
my $min_AF = 0.0 ;
my $max_AF = 1.0 ;

my $lower ;
my $upper ;
my $middle ;

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
	$AFS{$af/$Sample_size}++ ;
}

foreach my $af (sort {$a <=> $b} keys %AFS){
	print $af, "\t", $AFS{$af}, "\n" ;
}


exit ;

