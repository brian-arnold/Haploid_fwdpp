#!usr/bin/perl
use warnings ;
use strict ;
use lib "/n/holyscratch/bomblies_lab/bjarnold/List-MoreUtils-0.33/lib" ;
use List::Util qw(first max maxstr min minstr reduce shuffle sum) ;

#############################################
# THIS SCRIPT LOOKS THROUGH A FILE OF MS-LIKE
# OUTPUT, ANALYZES A SINGLE REPLICATE
# CATEGORIZES DATA BY NONSYNONYMOUS SNP DENSITY
#############################################

my $results_file = $ARGV[0] ;
my $results2_file = $ARGV[1] ;
my $selpos_file = $ARGV[2] ;
my $rep = $ARGV[3] ;
my $quantiles = 5 ;
unless(-e "rep${rep}_summaries"){
	system("mkdir rep${rep}_summaries") ;
}
open QC, ">./rep${rep}_summaries/QC_PrivateSNPs.txt" ;

my $pop_size ;
my $Sample_size ;
my $num_reps ;
my $num_sites ;
my $rec_rate ;
my $theta ;
my $gene_conv_rate ;
my $TrLen ;

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
my %Seq_data2 ; #  $Seq_data{index} = scalar of 0's and 1's
my %Segsites_hash ; # $Segsites{site} ++
my %Segsites2_hash ; # $Segsites{site} ++
my %Segsites ; # $Segsites{site} ++
my %Segsites2 ; # $Segsites{site} ++
my %Segsites_Total_hash ; # $Segsites{site} ++
my %Segsites_Total ; # $Segsites{site} ++

my %SumStat_Results ;

my %SelPos_hash ;
my %SelPos ;

#############################################
# ASSEMBLE SEQUENCES PER GENE PER IND
#############################################
my $replicate = 0 ;
open IN, "<./${results_file}" ;
my $num_segsites_raw = 0 ;
while(<IN>){
	chomp $_ ;
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
				$Segsites_hash{$line[$x]}++ ; # to check if sel positions exist
			}
		
			foreach my $ind (0..($Sample_size-1)){
				my $x = <IN> ; chomp $x ;
				$Seq_data{ $ind } = $x ;
			}
		}
	}
}
close IN ;

$replicate = 0 ;
open IN, "<./${results2_file}" ;
$num_segsites_raw = 0 ;
while(<IN>){
	chomp $_ ;
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
				$Segsites2_hash{$line[$x]}++ ; # to check if sel positions exist
			}
		
			foreach my $ind (0..($Sample_size-1)){
				my $x = <IN> ; chomp $x ;
				$Seq_data2{ $ind } = $x ;
			}
		}
	}
}
unless($Sample_size){
	print "didn't load simulation parameters properly; no sample size\n" ;
	exit ;
}
close IN ;


print QC "SampleSize_resultsFile: ", $Sample_size, "\n" ;
print QC "######## number of individuals with seqs\n" ;
print QC "N1\t", scalar keys %Seq_data, "\n" ;
print QC "N2\t", scalar keys %Seq_data2, "\n" ;
print QC "Num segsites N1\t", scalar keys %Segsites_hash, "\n" ;
print QC "Num segsites N2\t", scalar keys %Segsites2_hash, "\n" ;

$replicate = 0 ;
open IN, "<./${selpos_file}" ;
while(<IN>){
	chomp $_ ;
	if($_ =~ m/\/\//){ # beginning of replicate
		$replicate++ ;
		if($replicate == $rep){
			while(my $line = <IN>){
				chomp $line ;
				if($line !~ m/POSITIONS/ && $line !~ m/\/\//){
					my @info = split(/:/, $line) ;
					my $pos = $info[0] ;
					my $freq = $info[2] ;
					if($freq > 0){
						$SelPos_hash{ $pos }++ ;
					}
				}
			}
		}
	}
}
close IN ;
print QC "Num SelPos\t", scalar keys %SelPos_hash, "\n" ;

#print "SELPOS\t", scalar keys %SelPos_tmphash, "\n" ;
#############################################
# DELETE SELECTED POSITIONS [0,1] NOT IN SAMPLE
# this helps prevent extra selected positions from 
# getting created when making them into integers later
#############################################
foreach my $pos ( keys %SelPos_hash ){
	if( !exists $Segsites_hash{$pos} && !exists $Segsites2_hash{$pos}){
		delete $SelPos_hash{$pos} ;
	}
}
#print "SELPOS\t", scalar keys %SelPos_tmphash, "\n" ;

#############################################
# MAKE TOTAL SEGSITES HASH TO REMOVE MULTIALLELIC SITES
#############################################
#positions need to be ordered, so make master hash first
#positions can be shared among populations, but proximal ones
#on chromosome will be interpreted as multiallelic
foreach my $pos ( sort{$a <=> $b} keys %Segsites_hash ){
	$Segsites_Total_hash{$pos}++ ;
}
foreach my $pos ( sort{$a <=> $b} keys %Segsites2_hash ){
	$Segsites_Total_hash{$pos}++ ;
}
my $cnt2=0 ;
foreach my $pos ( sort{$a <=> $b} keys %Segsites_Total_hash ){
	$Segsites_Total{$cnt2} = int($pos*$num_sites) ;
	$cnt2++ ;
}


print "Segsites_Total_hash b4 filter: ", scalar keys %Segsites_Total_hash, "\n" ;
print "Segsites_Total b4 filter: ", scalar keys %Segsites_Total, "\n" ;
#############################################
### delete any seg site index which has same site, multiallelic
#############################################
my %SitesToExclude ;
foreach my $i ( 0.. (scalar keys %Segsites_Total)-2 ){
	if(exists $Segsites_Total{$i}){
		my $j = $i+1 ;
		if($Segsites_Total{$j} == $Segsites_Total{$i}){
			while($Segsites_Total{$j} == $Segsites_Total{$i}){
				$j++ ; #to capture potential runs of multiple hits at same site
				if(!exists $Segsites_Total{$j}){ #reached end
					last ;
				}
			}
			## delete $i through $j-1 (-1 b.c. while loop added extra)
			foreach my $x ($i..($j-1)){
				$SitesToExclude{$Segsites_Total{$x}}++ ;
				delete $Segsites_Total{$x} ;
			}
			$multi_allele_sites{($j-1-$i)+1} ++ ;
		}else{
			next ;
		}
	}
}
print "Segsites_Total aft filter: ", scalar keys %Segsites_Total, "\n" ;
foreach my $pos ( sort{$a <=> $b} keys %Segsites_hash ){
	my $chromopos = int($pos*$num_sites) ;
	if(!exists $SitesToExclude{$chromopos} ){
		$Segsites{$chromopos}++ ;
	}
}
foreach my $pos ( sort{$a <=> $b} keys %Segsites2_hash ){
	my $chromopos = int($pos*$num_sites) ;
	if(!exists $SitesToExclude{$chromopos} ){
		$Segsites2{$chromopos}++ ;
	}
}
foreach my $pos ( keys %SelPos_hash ){
	my $chromopos = int($pos*$num_sites) ;
	if(!exists $SitesToExclude{$chromopos} ){
		$SelPos{ $chromopos }++ ;
	}
}

print "segsites b4 and after filtering \n" ;
print scalar keys %Segsites_hash, "\t", scalar keys %Segsites, "\n" ;
print scalar keys %Segsites2_hash, "\t", scalar keys %Segsites2, "\n" ;

my $SharedPoly = 0 ;
my $PrivatePolyN1 = 0 ;
my $PrivatePolyN2 = 0 ;
foreach my $chromopos (sort{$a <=> $b} keys %Segsites){
	if(exists $Segsites2{$chromopos} ){
		$SharedPoly++
	}else{
		$PrivatePolyN1++ ;
	}
}
foreach my $chromopos (sort{$a <=> $b} keys %Segsites2){
	if(!exists $Segsites{$chromopos} ){
		$PrivatePolyN2++
	}
}

open OUT, ">rep${rep}_summaries/PrivateSharedPoly.txt" ;
print OUT "SharedPoly_Density\t", $SharedPoly/$num_sites, "\n" ;
print OUT "PrivatePolyN1_Density\t", $PrivatePolyN1/$num_sites, "\n" ;
print OUT "PrivatePolyN2_Density\t", $PrivatePolyN2/$num_sites, "\n" ;

close QC ;


exit ;


