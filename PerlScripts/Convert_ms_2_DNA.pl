#!usr/bin/perl
use warnings ;
use strict ;

#############################################
# THIS SCRIPT LOOKS THROUGH A FILE OF MS-LIKE
# OUTPUT, CONVERTS TO DNA SEQUENCE DATA
#############################################

my $results_file = $ARGV[0] ;
my $selpos_file = $ARGV[1] ;
my $rep = $ARGV[2] ;

#open QC, ">./QC.txt" ;

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
my %Segsites_tmphash ; # $Segsites{site} ++
my %Segsites_Bi_PosInChromo; ; # $Segsites{index} = site ;
my %Segsites_Bi_PosInGenostring ; # for mapping back to positin in string of 0's and 1's
my %AFS ; # $AFS{1}{site} = AF ; 
my %AFS_forPrinting ;
my %Pairwise_Dprime  ; # $Pairwise_Dprime{site1}{site2} = D' ;

my %Segsites_filtered_positionsInDNA ;

my %SelPos ;
my %SelPos_tmphash ;

#############################################
# ASSEMBLE SEQUENCES PER GENE PER IND
# WHICH PROGRAM WAS USED?
#############################################
my $replicate = 0 ;
open IN, "<./${results_file}" ;
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
				$Segsites{($x-1)} = ($line[$x]) ;
				$Segsites_tmphash{$line[$x]}++ ; # to check if sel positions exist
			}
		
			foreach my $ind (0..($Sample_size-1)){
				my $x = <IN> ; chomp $x ;
				$Seq_data{ $ind } = $x ;
			}
		}
	}
}
close IN ;

#print QC "SampleSize_resultsFile: ", $Sample_size, "\n" ;
#print QC "######## number of individuals with seqs\n" ;
#print QC scalar keys %Seq_data, "\n" ;
#print QC "num segsites ", scalar keys %Segsites, "\n" ;
#foreach my $key (sort{$a<=>$b} keys %Segsites){
#	print QC $key, "\t", $Segsites{$key}, "\n" ;
#}
#print QC "\n" ;

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
						$SelPos_tmphash{ $pos }++ ;
					}
				}
			}
		}
	}
}
close IN ;

=begin
foreach my $site (sort {$a <=> $b} keys %Segsites){
	print $site, "\t", $Segsites{$site}, "\t" ;
	if(exists $SelPos{$Segsites{$site}}){
		print "SELECTED\n" ;
	}else{
		print "\n" ;
	}
}
=cut

#print "SELPOS\t", scalar keys %SelPos_tmphash, "\n" ;
#############################################
# DELETE SELECTED POSITIONS [0,1] NOT IN SAMPLE
# this helps prevent extra selected positions from 
# getting created when making them into integers later
#############################################
foreach my $pos ( keys %SelPos_tmphash ){
	if( !exists $Segsites_tmphash{$pos} ){
		delete $SelPos_tmphash{$pos} ;
	}
}
#print "SELPOS\t", scalar keys %SelPos_tmphash, "\n" ;

#############################################
# CONVERT POSITIONS TO INTEGERS
#############################################
foreach my $cnt ( keys %Segsites ){
	$Segsites{$cnt} = int($Segsites{$cnt}*$num_sites) ;
}
foreach my $pos ( keys %SelPos_tmphash ){
	$SelPos{ int($pos*$num_sites) }++ ;
}


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
foreach my $key (sort{$a<=>$b} keys %Segsites){
    $Segsites_filtered_positionsInDNA{$Segsites{$key}} = $key ;
}
undef %Segsites ; # ensure it's not used again
#print QC "num segsites unfiltered ", scalar keys %Segsites_filtered_positionsInDNA, "\n" ;
#foreach my $key (sort{$a<=>$b} keys %Segsites_Bi_Unfiltered){
#	print QC $key, "\t", $Segsites_Bi_Unfiltered_indices[$key], "\t", $Segsites_Bi_Unfiltered{$key}, "\n" ;
#}

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

open OUT, ">./VariablePositions${rep}.txt" ;
foreach my $cnt ( sort{$a <=> $b} keys %Segsites ){
	print OUT $Segsites{$cnt}, "\t" ;
	if(exists $SelPos{$Segsites{$cnt}}){
		print OUT "SELECTED\n" ;
	}else{
		print OUT "NEUTRAL\n" ;
	}
}

close OUT ;


#close QC ;




exit ;
