#!usr/bin/perl
use warnings ;
use strict ;
use lib "/n/holyscratch/bomblies_lab/bjarnold/List-MoreUtils-0.33/lib" ;
use List::Util qw(first max maxstr min minstr reduce shuffle sum) ;

#############################################
# THIS SCRIPT LOOKS THROUGH A FILE OF MS-LIKE
# OUTPUT, ANALYZES A SINGLE REPLICATE
#############################################

my $results_file = $ARGV[0] ;
my $selpos_file = $ARGV[1] ;
my $max_mutationFreq = $ARGV[2] ;
my $rep = $ARGV[3] ;
my $gene_size = $ARGV[4] ;
my $quantiles = $ARGV[5] ;
my $Bootreps = 100 ;
unless(scalar @ARGV == 6){
	print "not enough args\n" ; exit ;
}

unless(-e "rep${rep}_summaries"){
	system("mkdir rep${rep}_summaries") ;
}
open QC, ">./rep${rep}_summaries/QC_bootstrapLD.txt" ;

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

my %SumStat_Results ;
my %PairWise_LD_vs_dist ;
my %PairWise_PC_vs_dist ;

my %SelPos ;
my %SelPos_tmphash ;

#############################################
# ASSEMBLE SEQUENCES PER GENE PER IND
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
$SumStat_Results{"S_star"} = $num_segsites_raw/$num_sites ;

if($max_mutationFreq > $Sample_size){
	$max_mutationFreq = $Sample_size ;
}

print QC "SampleSize_resultsFile: ", $Sample_size, "\n" ;
print QC "######## number of individuals with seqs\n" ;
print QC scalar keys %Seq_data, "\n" ;
print QC "num segsites ", scalar keys %Segsites, "\n" ;
#foreach my $key (sort{$a<=>$b} keys %Segsites){
#	print QC $key, "\t", $Segsites{$key}, "\n" ;
#}
print QC "\n" ;

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
=begin
foreach my $site (sort {$Segsites{$a} <=> $Segsites{$b}} keys %Segsites){
	print $site, "\t", $Segsites{$site}, "\t", $Segsites{$site}*$num_sites, "\t" ;
	if(exists $SelPos{$Segsites{$site}}){
		print "SELECTED\n" ;
	}else{
		print "\n" ;
	}
}
=cut

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
    $Segsites_Bi_PosInGenostring{$count} = $key ; # key is position of segsite in 0/1 string
    $Segsites_Bi_PosInChromo{$count} = $Segsites{$key} ; # position in chromosome
    $count++ ;
}
undef %Segsites ; # ensure it's not used again
print QC "num segsites unfiltered ", scalar keys %Segsites_Bi_PosInChromo, "\n" ;
#foreach my $key (sort{$a<=>$b} keys %Segsites_Bi_Unfiltered){
#	print QC $key, "\t", $Segsites_Bi_Unfiltered_indices[$key], "\t", $Segsites_Bi_Unfiltered{$key}, "\n" ;
#}

=begin
foreach my $site (sort {$Segsites_Bi_PosInChromo{$a} <=> $Segsites_Bi_PosInChromo{$b}} keys %Segsites_Bi_PosInChromo){
#foreach my $site (sort {$Segsites{$a} <=> $Segsites{$b}} keys %Segsites){
	print "FILTER\t", $site, "\t", $Segsites_Bi_PosInChromo{$site}, "\t", $Segsites_Bi_PosInChromo{$site}*$num_sites, "\t" ;
	if(exists $SelPos{$Segsites_Bi_PosInChromo{$site}}){
		print "SELECTED\n" ;
	}else{
		print "\n" ;
	}
}
=cut 

#############################################
## CONSTRUCT SINGLE SITE %AFS FOR SUMSTATS AND LD
#############################################
foreach my $seg_site_index ( keys %Segsites_Bi_PosInChromo ){
    my $site = $Segsites_Bi_PosInGenostring{$seg_site_index} ;
	my $af = 0 ;
	foreach my $ind (keys %Seq_data){
		if( substr($Seq_data{$ind}, $site, 1) == 1 ){
			$af++ ;
		}
	}
	$AFS{$Segsites_Bi_PosInGenostring{$seg_site_index}} = $af ;
	$AFS_forPrinting{$af/$Sample_size}++ ;
}

#############################################
##	MAKE FRAGMENT INTO GENES
#############################################
my %Functional_Effect_BiAllelic ;
my $numSynSNPs = 0 ;
my $numNonSynSNPs = 0 ;

foreach my $count ( sort{$a<=>$b} keys %Segsites_Bi_PosInChromo ){
	my $gene = int($Segsites_Bi_PosInChromo{$count}/$gene_size) ; # genes indexed from 0 to ...
	my $chromo_site = $Segsites_Bi_PosInChromo{$count} ;
	my $genostring_site = $Segsites_Bi_PosInGenostring{$count} ;

	if( exists $SelPos{$chromo_site} ){
		$Functional_Effect_BiAllelic{$gene}{$genostring_site} = "NONSYNONYMOUS" ;
		$numNonSynSNPs++ ;
	}else{
		$Functional_Effect_BiAllelic{$gene}{$genostring_site} = "SYNONYMOUS" ;
		$numSynSNPs++ ;
	}
}
print QC "SYN_SNPS: ", $numSynSNPs, "\n" ;
print QC "NONSYN_SNPS: ", $numNonSynSNPs, "\n" ;






open OUTWITHIN, ">./rep${rep}_summaries/WithinGeneD_Bootreps.txt" ;	
print OUTWITHIN "Quant", "\t", "Bootrep", "\t", "MeanSynD", "\t", "MeanNonsynD", "\n" ;

open OUTBETWEEN, ">./rep${rep}_summaries/BetweenGeneD_Bootreps.txt" ;	
print OUTBETWEEN "Quant", "\t", "Bootrep", "\t", "MeanSynD", "\t", "MeanNonsynD", "\n" ;

foreach my $bootrep ( 1 .. $Bootreps ){

	## due to nonsyn SNP density ties between genes
	## quantiling genes needs to be done within this loop
	## to guarantee
	#############################################
	##	QUANTILE GENES BY SNP DENSITY
	#############################################
	my %SYN_SNP_densities ; # $SNP_densities{$gene} = density
	my %NONSYN_SNP_densities ; # $SNP_densities{$gene} = density
	my %Genes_Quantiled ; # $Genes_Quantiled{gene} = quantile
	my %Quant_means ;

	foreach my $gene( sort {$a <=> $b} keys %Functional_Effect_BiAllelic ){
		$SYN_SNP_densities{$gene} = 0 ;
		$NONSYN_SNP_densities{$gene} = 0 ;
		foreach my $genostring_site (sort {$a <=> $b} keys %{$Functional_Effect_BiAllelic{$gene}} ){
			if( $Functional_Effect_BiAllelic{$gene}{$genostring_site} eq "NONSYNONYMOUS"){
				$NONSYN_SNP_densities{$gene} += 1/$gene_size ;
			}else{
				$SYN_SNP_densities{$gene} += 1/$gene_size ;
			}
		}
	}

	my $counter = 0 ;
	my $numGenesPerQuantile = (scalar keys %NONSYN_SNP_densities)/$quantiles ;
	foreach my $gene ( sort{ $NONSYN_SNP_densities{$a} <=> $NONSYN_SNP_densities{$b} } keys %NONSYN_SNP_densities){
		my $quant = int($counter/$numGenesPerQuantile) ;
		$Genes_Quantiled{$gene} = $quant ;
		$Quant_means{$quant} += $NONSYN_SNP_densities{$gene}/$numGenesPerQuantile ;
		$counter++ ;
	}
	my $tmp_sum = 0 ;
	foreach my $gene (sort{$a <=> $b} keys %NONSYN_SNP_densities){
		if( $Genes_Quantiled{$gene} == 0 ){
			$tmp_sum += $gene ;
		}
	}	
	#print $tmp_sum, "\n" ;
=begin	
	print QC "GENE\tQUANT\tNONSYNSNPS\tSYNSNPS\n" ;
	foreach my $quant ( 0 .. $quantiles-1 ){
		foreach my $gene (sort{$a <=> $b} keys %NONSYN_SNP_densities){
			if( $Genes_Quantiled{$gene} == $quant ){
				print QC $gene, "\t", $quant, "\t", $NONSYN_SNP_densities{$gene}*$gene_size, "\t" ;
				if(exists $SYN_SNP_densities{$gene} ){
					print QC $SYN_SNP_densities{$gene}*$gene_size, "\n" ;
				}else{
					print QC "0\n" ;
				}
			}
		}
	}
=cut

=begin
	foreach my $gene (sort {$a <=> $b} keys %Functional_Effect_BiAllelic){
		print "######GENE: ", $gene, "\t", $Genes_Quantiled{$gene}, "\n" ;
		print $SYN_SNP_densities{$gene}, "\n" ;
		#foreach my $genostring_site (sort {$a <=> $b} keys %{$Functional_Effect_BiAllelic{$gene}}){
		#	print $genostring_site, "\t", $Functional_Effect_BiAllelic{$gene}{$genostring_site}, "\n" ;
		#}
	}
=cut
	#print QC "Quantile\tMean\n" ;
	#foreach my $quant (keys %Quant_means){
	#	print QC $quant, "\t", $Quant_means{$quant}, "\n" ;
	#}

	#############################################
	##	CALC RESAMPLE SIZE PER ALLELE FREQ
	#############################################
	my %WithinGeneD ;
	my %BetweenGeneD ;
	my %TotalNumberSNPs ; # TotalNumberSNPs{$quant} = total num nonsyn for this quantile
	my %Nonsyn_SNP_indexer ; # @{$Nonsyn_SNP_indexer{dNdS_quantile}{$freq}{counter}} = ($gene, $site, $freq, $base) ;
	my %Syn_SNP_indexer ; # @{$Syn_SNP_indexer{dNdS_quantile}{$freq}{counter}} = ($gene, $site, $freq, $base) ;
	############
	my %nonsyn_counter ; # $nonsyn_counter{$quant}{$freq}
	my %syn_counter  ; # $syn_counter{$freq}
	############
	my %ResampleSizePerAF ; # $ResampleSizePerAF{$quant}{$freq}


	foreach my $freq (1 .. $max_mutationFreq){
		foreach my $quant (0 .. $quantiles-1){
			$nonsyn_counter{$quant}{$freq} = 0 ;
			$syn_counter{$quant}{$freq} = 0 ;
		}
	}
	foreach my $gene (keys %Functional_Effect_BiAllelic){
		my $quant = $Genes_Quantiled{$gene} ;
		foreach my $genostring_site (keys %{$Functional_Effect_BiAllelic{$gene}}){
			if($AFS{$genostring_site} <= $max_mutationFreq){
				my $freq = $AFS{$genostring_site} ;
				if( $Functional_Effect_BiAllelic{$gene}{$genostring_site} eq "NONSYNONYMOUS" ){
					@{$Nonsyn_SNP_indexer{ $quant }{$freq}{ $nonsyn_counter{$quant}{$freq} } } = ($gene, $genostring_site, $freq) ;
					$nonsyn_counter{$quant}{$freq}++ ;
					#$TotalNonSynSFS{$quant}{$freq}++ ;
				}elsif( $Functional_Effect_BiAllelic{$gene}{$genostring_site} eq "SYNONYMOUS" ){
					@{$Syn_SNP_indexer{$quant}{$freq}{ $syn_counter{$quant}{$freq} } } = ($gene, $genostring_site, $freq) ;
					$syn_counter{$quant}{$freq}++ ;
					#$TotalSynSFS{$freq}++ ;
				}					
			}
		}
	}
	foreach my $freq (1 .. $max_mutationFreq){
		foreach my $quant (0 .. $quantiles-1){
			if(exists $syn_counter{$quant}{$freq} && exists $nonsyn_counter{$quant}{$freq}){
				if( $syn_counter{$quant}{$freq} > $nonsyn_counter{$quant}{$freq}){
					$ResampleSizePerAF{$quant}{$freq} = $nonsyn_counter{$quant}{$freq} ;
				}else{
					$ResampleSizePerAF{$quant}{$freq} = $syn_counter{$quant}{$freq} ;
				}
			}else{
				$ResampleSizePerAF{$quant}{$freq} = 0 ;
			}
		}
	}

=begin
	print QC "Quant\tfreq\tNumSyn\tNumNonSyn\tResampleSize\n" ;
	foreach my $quant (sort{$a <=> $b} keys %syn_counter){
		foreach my $freq (sort{$a <=> $b} keys %{$syn_counter{$quant}}){
			print QC $quant, "\t", $freq, "\t", $syn_counter{$quant}{$freq}, "\t" ;
			if(exists $nonsyn_counter{$quant}{$freq}){
				print QC $nonsyn_counter{$quant}{$freq}, "\t" ;
			}else{
				print QC "0\t" ;
			}
			print QC $ResampleSizePerAF{$quant}{$freq}, "\n" ;
		}
	}
=cut




	my %nonsyn_resamp ;
	my %syn_resamp ;
	foreach my $quant ( 0 .. $quantiles-1 ){
		# SUBSAMPLE NONSYN AND SYN SNPs
		
		## bootstrapping nonsyn muts, so resample these and equal number of syn
		
		## resamp size as lower of NONSYN or SYN count		
		foreach my $freq (1 .. $max_mutationFreq){
			if( $ResampleSizePerAF{$quant}{$freq} > 0 ){
				my @indices_left_nonsyn = ( 0 .. ((scalar keys %{$Nonsyn_SNP_indexer{$quant}{$freq}}) - 1) ) ;	
				my @indices_left_syn = ( 0 .. ((scalar keys %{$Syn_SNP_indexer{$quant}{$freq}}) - 1) ) ;	
				foreach ( 1 .. $ResampleSizePerAF{$quant}{$freq} ){ ## sample size could be determined by syn or nonsyn, depending on which is lower
					my $tmp = int(rand( scalar @indices_left_nonsyn )) ;
					my $rand_index = $indices_left_nonsyn[$tmp] ;
					push @{$nonsyn_resamp{$quant}{$freq}}, $rand_index ;
					## for quantile determining $resamp_size, this will delete all indices
					splice @indices_left_nonsyn, $tmp, 1 ;
				
					$tmp = int(rand( scalar @indices_left_syn )) ;
					$rand_index = $indices_left_syn[$tmp] ;
					push @{$syn_resamp{$quant}{$freq}}, $rand_index ;
					## for quantile determining $resamp_size, this will delete all indices
					splice @indices_left_syn, $tmp, 1 ;
				}		
			}
		}

		# Now that you have random indices, calculate mutation burden
		foreach my $freq (1 .. $max_mutationFreq){
			if( $ResampleSizePerAF{$quant}{$freq} > 0 ){
				my $lower ;
				my $upper ;
				for($lower=0; $lower < ((scalar @{$nonsyn_resamp{$quant}{$freq}})-1); $lower++ ){
					my $index1 = ${$nonsyn_resamp{$quant}{$freq}}[$lower] ;
					my $gene1_rand = ${$Nonsyn_SNP_indexer{ $quant }{$freq}{$index1}}[0] ;
					my $site1_rand = ${$Nonsyn_SNP_indexer{ $quant }{$freq}{$index1}}[1] ;
					my $freq1_rand = ${$Nonsyn_SNP_indexer{ $quant }{$freq}{$index1}}[2] ;			
				
					for ($upper = ($lower+1); $upper<(scalar @{$nonsyn_resamp{$quant}{$freq}}); $upper++){
						my $index2 = ${$nonsyn_resamp{$quant}{$freq}}[$upper] ;
						my $gene2_rand = ${$Nonsyn_SNP_indexer{ $quant }{$freq}{$index2}}[0] ;
						my $site2_rand = ${$Nonsyn_SNP_indexer{ $quant }{$freq}{$index2}}[1] ;
						my $freq2_rand = ${$Nonsyn_SNP_indexer{ $quant }{$freq}{$index2}}[2] ;			
				
						## CALCULATE LD
						my $site1_site2_dblDAF_haplo = 0 ;
						my $site1_DAF = 0 ;
						my $site2_DAF = 0 ;
						foreach my $ind (keys %Seq_data){
							my $b1 = substr($Seq_data{$ind}, $site1_rand, 1) ;
							my $b2 = substr($Seq_data{$ind}, $site2_rand, 1) ;
							if( $b1 == 1 ){
								$site1_DAF++ ;
								if( $b2 == 1 ){
									$site1_site2_dblDAF_haplo++ ;
								}
							}
							if( $b2 == 1 ){
								$site2_DAF++ ;
							}
						}
						$site1_DAF = $site1_DAF/$Sample_size ;
						$site2_DAF = $site2_DAF/$Sample_size ;
						$site1_site2_dblDAF_haplo = $site1_site2_dblDAF_haplo/$Sample_size ;
						my $LD = $site1_site2_dblDAF_haplo - ($site1_DAF*$site2_DAF) ;
						if($gene1_rand eq $gene2_rand){
							push @{$WithinGeneD{$quant}{$bootrep}{"NONSYNONYMOUS"}}, $LD ;
						}else{
							push @{$BetweenGeneD{$quant}{$bootrep}{"NONSYNONYMOUS"}}, $LD ;
						}
					}
				}
				
				
				for($lower=0; $lower < ((scalar @{$syn_resamp{$quant}{$freq}})-1); $lower++ ){
					my $index1 = ${$syn_resamp{$quant}{$freq}}[$lower] ;
					my $gene1_rand = ${$Syn_SNP_indexer{ $quant }{$freq}{$index1}}[0] ;
					my $site1_rand = ${$Syn_SNP_indexer{ $quant }{$freq}{$index1}}[1] ;
					my $freq1_rand = ${$Syn_SNP_indexer{ $quant }{$freq}{$index1}}[2] ;			
				
					for ($upper = ($lower+1); $upper<(scalar @{$syn_resamp{$quant}{$freq}}); $upper++){
						my $index2 = ${$syn_resamp{$quant}{$freq}}[$upper] ;
						my $gene2_rand = ${$Syn_SNP_indexer{ $quant }{$freq}{$index2}}[0] ;
						my $site2_rand = ${$Syn_SNP_indexer{ $quant }{$freq}{$index2}}[1] ;
						my $freq2_rand = ${$Syn_SNP_indexer{ $quant }{$freq}{$index2}}[2] ;			
				
						## CALCULATE LD
						my $site1_site2_dblDAF_haplo = 0 ;
						my $site1_DAF = 0 ;
						my $site2_DAF = 0 ;
						foreach my $ind (keys %Seq_data){
							my $b1 = substr($Seq_data{$ind}, $site1_rand, 1) ;
							my $b2 = substr($Seq_data{$ind}, $site2_rand, 1) ;
							if( $b1 == 1 ){
								$site1_DAF++ ;
								if( $b2 == 1 ){
									$site1_site2_dblDAF_haplo++ ;
								}
							}
							if( $b2 == 1 ){
								$site2_DAF++ ;
							}
								
							
						}
						$site1_DAF = $site1_DAF/$Sample_size ;
						$site2_DAF = $site2_DAF/$Sample_size ;
						$site1_site2_dblDAF_haplo = $site1_site2_dblDAF_haplo/$Sample_size ;
						my $LD = $site1_site2_dblDAF_haplo - ($site1_DAF*$site2_DAF) ;
						if($gene1_rand eq $gene2_rand){
							push @{$WithinGeneD{$quant}{$bootrep}{"SYNONYMOUS"}}, $LD ;
						}else{
							push @{$BetweenGeneD{$quant}{$bootrep}{"SYNONYMOUS"}}, $LD ;
						}
					}
				}
			}
		}
	}

	foreach my $quant ( sort{$a <=> $b} keys %WithinGeneD ){
		foreach my $bootrep ( sort{$a <=> $b} keys %{$WithinGeneD{$quant}} ){
			print OUTWITHIN $quant, "\t", $bootrep, "\t" ;
			my $sum=0 ;
			foreach ( @{$WithinGeneD{$quant}{$bootrep}{"SYNONYMOUS"}} ){
				$sum += $_ ;
			}
			print OUTWITHIN $sum/(scalar @{$WithinGeneD{$quant}{$bootrep}{"SYNONYMOUS"}} ), "\t" ;
			$sum=0 ;
			foreach ( @{$WithinGeneD{$quant}{$bootrep}{"NONSYNONYMOUS"}} ){
				$sum += $_ ;
			}		
			print OUTWITHIN $sum/(scalar @{$WithinGeneD{$quant}{$bootrep}{"NONSYNONYMOUS"}} ), "\n"  ;		
		}
	}

	foreach my $quant ( sort{$a <=> $b} keys %BetweenGeneD ){
		foreach my $bootrep ( sort{$a <=> $b} keys %{$BetweenGeneD{$quant}} ){
			print OUTBETWEEN $quant, "\t", $bootrep, "\t" ;
			my $sum=0 ;
			foreach ( @{$BetweenGeneD{$quant}{$bootrep}{"SYNONYMOUS"}} ){
				$sum += $_ ;
			}
			print OUTBETWEEN $sum/(scalar @{$BetweenGeneD{$quant}{$bootrep}{"SYNONYMOUS"}} ), "\t" ;
			$sum=0 ;
			foreach ( @{$BetweenGeneD{$quant}{$bootrep}{"NONSYNONYMOUS"}} ){
				$sum += $_ ;
			}
			print OUTBETWEEN $sum/(scalar @{$BetweenGeneD{$quant}{$bootrep}{"NONSYNONYMOUS"}} ), "\n"  ;		
		}
	}	
}
	
close OUTWITHIN ;
close OUTBETWEEN ;
	

exit ;








