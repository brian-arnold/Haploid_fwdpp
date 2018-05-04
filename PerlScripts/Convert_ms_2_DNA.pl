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

#############################################
# DATA STRUCTURES
#############################################
my $sum ;
my $sq_sum ;

my %multi_allele_sites ; #  $multi_allele_sites{num hits} = count ;
my %Seq_data ; #  $Seq_data{index} = scalar of 0's and 1's
my %Segsites ; # $Segsites{index} = site ;
my %Segsites_filtered_positionsInDNA ;


#############################################
# ASSEMBLE SEQUENCES PER GENE PER IND
#############################################
my $replicate = 0 ;
open IN, "<./${file}" ;
my $num_segsites_raw = 0 ;
while(<IN>){
	chomp $_ ;
	if($_ =~ m/^\.\/haploid_ind_neut/ ){ # command used with arguments
		my @line = split(" ", $_) ;
		$pop_size = $line[1] ;
		$Sample_size = $line[7] ;
		$num_reps = $line[8] ;
		$theta = $line[2] ;
		$rec_rate = $line[3] ;
		$num_sites = $line[4] ;
		$TrLen = $line[5] ;
	}
	if($_ =~ m/^\.\/haploid_ind_seln/ ){ # command used with arguments
		my @line = split(" ", $_) ;
		$Sample_size = $line[9] ;
		$num_reps = $line[10] ;
		$num_sites = $line[5] ;
	}
	if($_ =~ m/^\.\/ms/ ){ # command used with arguments
		my @line = split(" ", $_) ;
		$Sample_size = $line[1] ;
		$num_reps = $line[2] ;
		$theta = $line[4] ;
		$rec_rate = $line[9] ;
		$num_sites = $line[7] ;
		$TrLen = $line[10] ;
	}
	if($_ =~ m/^\.\/haploid_struct_neutral/ ){ # command used with arguments
		my @line = split(" ", $_) ;
		$Sample_size = $line[10] ;
		$num_reps = $line[11] ;
		$theta = $line[5] ;
		$rec_rate = $line[6] ;
		$num_sites = $line[7] ;
		$TrLen = $line[8] ;
	}
	if($_ =~ m/^\.\/haploid_struct_seln/ ){ # command used with arguments
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
			# load positions
			my $positionsline = <IN> ; chomp $positionsline ;
			@line = split(" ", $positionsline) ;
			foreach my $x ( 1..$num_segsites_raw ){
				$Segsites{($x-1)} = (int($line[$x]*$num_sites)) ;
			}
			# load sequences
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
#open OUT, ">./positions${rep}.fasta" ;
my $count = 0 ;
foreach my $key (sort{$a<=>$b} keys %Segsites){    
	$Segsites_filtered_positionsInDNA{$Segsites{$key}} = $key ;  # key is position of segsite in 0/1 string
	#print OUT $Segsites{$key}, "\n" ;
	$count++ ;
}
undef %Segsites ; # ensure it's not used again
#close OUT ;


open OUT, ">./DNA${rep}.fasta" ;
foreach my $ind (sort {$a <=> $b} keys %Seq_data ){
	print OUT ">Individual_", $ind, "\n" ;
	foreach my $DNApos (0 .. $num_sites-1){ # start at 0 since int() cuts off decimal when making position
		if(exists $Segsites_filtered_positionsInDNA{$DNApos}){
			my $posInSegsites = $Segsites_filtered_positionsInDNA{$DNApos} ;
			# the value of  $Positions{$pos} is the index of the array corresponding to that position
			if( substr($Seq_data{$ind}, $posInSegsites, 1) == 1 ){
				print OUT "G" ;
			}else{
				print OUT "C" ;
			}					
		}else{
			print OUT "A" ;
		}
	}
	print OUT "\n" ;
}
close OUT ;




exit ;
