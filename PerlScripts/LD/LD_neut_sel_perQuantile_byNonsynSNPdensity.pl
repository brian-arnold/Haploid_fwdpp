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
my $selpos_file = $ARGV[1] ;
my $rep = $ARGV[2] ;
my $gene_size =  $ARGV[3] ;
my $quantiles = $ARGV[4] ;
unless(-e "rep${rep}_summaries"){
	system("mkdir rep${rep}_summaries") ;
}
open QC, ">./rep${rep}_summaries/QC.txt" ;
print QC "Gene size:\t", $gene_size, "\n" ;
print QC "Quantiles:\t", $quantiles, "\n" ;

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
	if($_ =~ m/haploid_ind_seln_epi/ ){ # command used with arguments
		my @line = split(" ", $_) ;
		$Sample_size = $line[10] ;
		$num_reps = $line[11] ;
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
my %GenostringSite_2_Gene ;

foreach my $count ( sort{$a<=>$b} keys %Segsites_Bi_PosInChromo ){
	my $gene = int($Segsites_Bi_PosInChromo{$count}/$gene_size) ; # genes indexed from 0 to ...
	my $chromo_site = $Segsites_Bi_PosInChromo{$count} ;
	my $genostring_site = $Segsites_Bi_PosInGenostring{$count} ;

	$GenostringSite_2_Gene{$genostring_site} = $gene ;
	if( exists $SelPos{$chromo_site} ){
		$Functional_Effect_BiAllelic{$gene}{$genostring_site} = "NONSYNONYMOUS" ;
	}else{
		$Functional_Effect_BiAllelic{$gene}{$genostring_site} = "SYNONYMOUS" ;
	}
}

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

=begin
foreach my $gene (sort {$a <=> $b} keys %Functional_Effect_BiAllelic){
	print "######GENE: ", $gene, "\t", $Genes_Quantiled{$gene}, "\n" ;
	print $SYN_SNP_densities{$gene}, "\n" ;
	#foreach my $genostring_site (sort {$a <=> $b} keys %{$Functional_Effect_BiAllelic{$gene}}){
	#	print $genostring_site, "\t", $Functional_Effect_BiAllelic{$gene}{$genostring_site}, "\n" ;
	#}
}
=cut
print QC "Quantile\tMean\n" ;
foreach my $quant (keys %Quant_means){
	print QC $quant, "\t", $Quant_means{$quant}, "\n" ;
}

#############################################
##	MEAN LD PER QUANTILE
#############################################



my %PairWise_D_SYN_WithinGene ; # @{$PairWise_D_SYN{$dist}}, ($LD) ;
my %PairWise_Dprime_SYN_WithinGene ; # @{$PairWise_Dprime_SYN{$dist}}, $compatibility ;
my %PairWise_D_NONSYN_WithinGene ; # @{$PairWise_D_vs_dist_SYN{$dist}}, ($LD) ;
my %PairWise_Dprime_NONSYN_WithinGene ; # @{$PairWise_Dprime_SYN{$dist}}, $compatibility ;

my %PairWise_D_SYN_BetweenGene ; # @{$PairWise_D_SYN{$dist}}, ($LD) ;
my %PairWise_Dprime_SYN_BetweenGene ; # @{$PairWise_Dprime_SYN{$dist}}, $compatibility ;
my %PairWise_D_NONSYN_BetweenGene ; # @{$PairWise_D_vs_dist_SYN{$dist}}, ($LD) ;
my %PairWise_Dprime_NONSYN_BetweenGene ; # @{$PairWise_Dprime_SYN{$dist}}, $compatibility ;
	
my %PairWise_D_SYN_ByDist ; # PairWise_D_SYN_ByDist{$quant}{$dist} ;
my %PairWise_Dprime_SYN_ByDist ;
my %PairWise_D_NONSYN_ByDist ;
my %PairWise_Dprime_NONSYN_ByDist ;

my @count_Syn ;
my @count_NonSyn ;
foreach my $cnt ( sort{$a <=> $b} keys %Segsites_Bi_PosInGenostring ){
	my $chromo_site = $Segsites_Bi_PosInChromo{$cnt} ;
	my $genostring_site = $Segsites_Bi_PosInGenostring{$cnt} ;
	if( exists $SelPos{$chromo_site} ){
		push @count_NonSyn, $genostring_site ;
	}else{
		push @count_Syn, $genostring_site ;
	}
}
print QC "\@count_Syn elements: ", scalar @count_Syn, "\n" ;
print QC "\@count_NonSyn elements: ", scalar @count_NonSyn, "\n";  

### SYNONYMOUS
for ($lower = 0; $lower<((scalar @count_Syn)-1); $lower++){
	my $gene1 = $GenostringSite_2_Gene{$count_Syn[$lower]} ;
	my $site1 = $count_Syn[$lower] ;

	for ($upper = ($lower+1); $upper<(scalar @count_Syn); $upper++){
		my $gene2 = $GenostringSite_2_Gene{$count_Syn[$upper]} ;
		my $site2 = $count_Syn[$upper] ;
		if( $Genes_Quantiled{$gene1} == $Genes_Quantiled{$gene2} ){
			my $quant = $Genes_Quantiled{$gene1} ;
	
			my $site1_site2_dblDAF_haplo = 0 ;
			my $site1_DAF = 0 ;
			my $site2_DAF = 0 ;
			# go through individuals of lower samp size gene, if they differ
			my $samp_size_for_pair = 0 ;
			foreach my $ind (keys %Seq_data){
				$samp_size_for_pair++ ;
				my $b1 = substr($Seq_data{$ind}, $site1, 1) ;
				my $b2 = substr($Seq_data{$ind}, $site2, 1) ;
				if( $b1 == 1){
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
			# if Dab positive, Dmax = min( Pa*(1-Pb), Pb*(1-Pa) )
			# if Dab negative, Dmax = min( Pa*(Pb), (1-Pa)*(1-Pb) )
			my $Dmax = 0 ;
			if ($LD>=0){
				if( $site1_DAF*(1-$site2_DAF) < $site2_DAF*(1-$site1_DAF) ){
					$Dmax = $site1_DAF*(1-$site2_DAF) ;
				}else{
					$Dmax = $site2_DAF*(1-$site1_DAF) ;
				}
			}elsif($LD<0){
				if( ($site1_DAF*$site2_DAF) < (1-$site1_DAF)*(1-$site2_DAF) ){
					$Dmax = $site1_DAF*$site2_DAF ;
				}else{
					$Dmax = (1-$site1_DAF)*(1-$site2_DAF) ;
				}
			}
			if($gene1 eq $gene2){
				push @{$PairWise_D_SYN_WithinGene{$quant}}, ($LD) ;
				push @{$PairWise_Dprime_SYN_WithinGene{$quant}}, ($LD)/$Dmax ;
				#push @{$PairWise_Rsq_vs_dist_SYN{$dist}}, ($LD**2)/($site1_DAF*(1-$site1_DAF)*$site2_DAF*(1-$site2_DAF))  ;
			}else{
				push @{$PairWise_D_SYN_BetweenGene{$quant}}, ($LD) ;
				push @{$PairWise_Dprime_SYN_BetweenGene{$quant}}, ($LD)/$Dmax ;
		
			}
		}
		
	}
}


### NONSYNONYMOUS
for ($lower = 0; $lower<((scalar @count_NonSyn)-1); $lower++){
	
	my $gene1 = $GenostringSite_2_Gene{$count_NonSyn[$lower]} ;
	my $site1 = $count_NonSyn[$lower] ;
		
	for ($upper = ($lower+1); $upper<(scalar @count_NonSyn); $upper++){
		my $gene2 = $GenostringSite_2_Gene{$count_NonSyn[$upper]} ;
		my $site2 = $count_NonSyn[$upper] ;

		if( $Genes_Quantiled{$gene1} == $Genes_Quantiled{$gene2} ){
			my $quant = $Genes_Quantiled{$gene1} ;
						
			my $site1_site2_dblDAF_haplo = 0 ;
			my $site1_DAF = 0 ;
			my $site2_DAF = 0 ;
			# go through individuals of lower samp size gene, if they differ
			foreach my $ind (keys %Seq_data){
				my $b1 = substr($Seq_data{$ind}, $site1, 1) ;
				my $b2 = substr($Seq_data{$ind}, $site2, 1) ;
				if( $b1 == 1){
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
			# if Dab positive, Dmax = min( Pa*(1-Pb), Pb*(1-Pa) )
			# if Dab negative, Dmax = min( Pa*(Pb), (1-Pa)*(1-Pb) )
			my $Dmax = 0 ;
			if ($LD>=0){
				if( $site1_DAF*(1-$site2_DAF) < $site2_DAF*(1-$site1_DAF) ){
					$Dmax = $site1_DAF*(1-$site2_DAF) ;
				}else{
					$Dmax = $site2_DAF*(1-$site1_DAF) ;
				}
			}elsif($LD<0){
				if( ($site1_DAF*$site2_DAF) < (1-$site1_DAF)*(1-$site2_DAF) ){
					$Dmax = $site1_DAF*$site2_DAF ;
				}else{
					$Dmax = (1-$site1_DAF)*(1-$site2_DAF) ;
				}
			}
			if($gene1 eq $gene2){
				push @{$PairWise_D_NONSYN_WithinGene{$quant}}, ($LD) ;
				push @{$PairWise_Dprime_NONSYN_WithinGene{$quant}}, ($LD)/$Dmax ;
			}else{
				push @{$PairWise_D_NONSYN_BetweenGene{$quant}}, ($LD) ;
				push @{$PairWise_Dprime_NONSYN_BetweenGene{$quant}}, ($LD)/$Dmax ;
			}
	
		}
	}
}

#################
## PRINT WITHIN GENE RESULTS
#################


my $total_D_SYN = 0 ;
my $total_D_SYN_tally = 0 ;
open OUT1, ">./rep${rep}_summaries/D.Dprime.vs.Quant.SYN_WithinGene" ;
foreach my $quant (sort{$a <=> $b} keys %PairWise_D_SYN_WithinGene){
	my $sum_D = 0 ;
	my $sum_Dpr = 0 ;
	my $sum_Rsq = 0 ;			
	foreach my $x ( @{$PairWise_D_SYN_WithinGene{$quant}} ){
		$sum_D += $x ;
		$total_D_SYN += $x ;
		$total_D_SYN_tally++ ;
	}
	foreach my $x ( @{$PairWise_Dprime_SYN_WithinGene{$quant}} ){
		$sum_Dpr += $x ;
	}
	#foreach my $x ( @{$PairWise_Rsq_vs_dist_SYN{$dist}} ){
	#	$sum_Rsq += $x ;
	#}			
	print OUT1 $quant, "\t", scalar @{$PairWise_D_SYN_WithinGene{$quant}}, "\t" ;
	print OUT1 $sum_D/(scalar @{$PairWise_D_SYN_WithinGene{$quant}}), "\t", $sum_Dpr/(scalar @{$PairWise_Dprime_SYN_WithinGene{$quant}}), "\n" ;
}
close OUT1 ;
print QC "AVG LD SYN WITHIN GENE: " ;
if($total_D_SYN_tally){
	print QC $total_D_SYN/$total_D_SYN_tally, "\n"
}else{
	print QC "NA\n"
}


#######
my $total_D_NONSYN = 0 ;
my $total_D_NONSYN_tally = 0 ;
open OUT1, ">./rep${rep}_summaries/D.Dprime.vs.Quant.NONSYN_WithinGene" ;
foreach my $quant (sort{$a <=> $b} keys %PairWise_D_NONSYN_WithinGene){
	my $sum_D = 0 ;
	my $sum_Dpr = 0 ;
	my $sum_Rsq = 0 ;			
	foreach my $x ( @{$PairWise_D_NONSYN_WithinGene{$quant}} ){
		$sum_D += $x ;
		$total_D_NONSYN += $x ;
		$total_D_NONSYN_tally++ ;
	}
	foreach my $x ( @{$PairWise_Dprime_NONSYN_WithinGene{$quant}} ){
		$sum_Dpr += $x ;
	}
	#foreach my $x ( @{$PairWise_Rsq_vs_dist_NONSYN_WithinGene{$quant}} ){
	#	$sum_Rsq += $x ;
	#}			
	print OUT1 $quant, "\t", scalar @{$PairWise_D_NONSYN_WithinGene{$quant}}, "\t" ;
	print OUT1 $sum_D/(scalar @{$PairWise_D_NONSYN_WithinGene{$quant}}), "\t", $sum_Dpr/(scalar @{$PairWise_Dprime_NONSYN_WithinGene{$quant}}), "\n" ;
}
close OUT1 ;

print QC "AVG LD NONSYN WITHIN GENE: " ;
if($total_D_NONSYN_tally){
	print QC $total_D_NONSYN/$total_D_NONSYN_tally, "\n"
}else{
	print QC "NA\n"
}

#################
## PRINT BETWEEN GENE RESULTS
#################
$total_D_SYN = 0 ;
$total_D_SYN_tally = 0 ;
open OUT1, ">./rep${rep}_summaries/D.Dprime.vs.Quant.SYN_BetweenGene" ;
foreach my $quant (sort{$a <=> $b} keys %PairWise_D_SYN_BetweenGene){
	my $sum_D = 0 ;
	my $sum_Dpr = 0 ;
	my $sum_Rsq = 0 ;			
	foreach my $x ( @{$PairWise_D_SYN_BetweenGene{$quant}} ){
		$sum_D += $x ;
		$total_D_SYN += $x ;
		$total_D_SYN_tally++ ;
	}
	foreach my $x ( @{$PairWise_Dprime_SYN_BetweenGene{$quant}} ){
		$sum_Dpr += $x ;
	}
	#foreach my $x ( @{$PairWise_Rsq_vs_dist_SYN{$dist}} ){
	#	$sum_Rsq += $x ;
	#}			
	print OUT1 $quant, "\t", scalar @{$PairWise_D_SYN_BetweenGene{$quant}}, "\t" ;
	print OUT1 $sum_D/(scalar @{$PairWise_D_SYN_BetweenGene{$quant}}), "\t", $sum_Dpr/(scalar @{$PairWise_Dprime_SYN_BetweenGene{$quant}}), "\n" ;
}
close OUT1 ;
print QC "AVG LD SYN BETWEEN GENE:\t" ;
if($total_D_SYN_tally){
	print QC $total_D_SYN/$total_D_SYN_tally, "\n"
}else{
	print QC "NA\n"
}
#######
$total_D_NONSYN = 0 ;
$total_D_NONSYN_tally = 0 ;
open OUT1, ">./rep${rep}_summaries/D.Dprime.vs.Quant.NONSYN_BetweenGene" ;
foreach my $quant (sort{$a <=> $b} keys %PairWise_D_NONSYN_BetweenGene){
	my $sum_D = 0 ;
	my $sum_Dpr = 0 ;
	my $sum_Rsq = 0 ;			
	foreach my $x ( @{$PairWise_D_NONSYN_BetweenGene{$quant}} ){
		$sum_D += $x ;
		$total_D_NONSYN += $x ;
		$total_D_NONSYN_tally++ ;
	}
	foreach my $x ( @{$PairWise_Dprime_NONSYN_BetweenGene{$quant}} ){
		$sum_Dpr += $x ;
	}
	#foreach my $x ( @{$PairWise_Rsq_vs_dist_NONSYN_BetweenGene{$quant}} ){
	#	$sum_Rsq += $x ;
	#}			
	print OUT1 $quant, "\t", scalar @{$PairWise_D_NONSYN_BetweenGene{$quant}}, "\t" ;
	print OUT1 $sum_D/(scalar @{$PairWise_D_NONSYN_BetweenGene{$quant}}), "\t", $sum_Dpr/(scalar @{$PairWise_Dprime_NONSYN_BetweenGene{$quant}}), "\n" ;
}
close OUT1 ;

print QC "AVG LD NONSYN BETWEEN GENE:\t" ;
if($total_D_NONSYN_tally){
	print QC $total_D_NONSYN/$total_D_NONSYN_tally, "\n"
}else{
	print QC "NA\n"
}

close QC ;


exit ;

