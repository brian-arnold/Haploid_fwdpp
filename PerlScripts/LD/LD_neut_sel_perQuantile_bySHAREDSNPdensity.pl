#!usr/bin/perl
use warnings ;
use strict ;

#############################################
# THIS SCRIPT LOOKS THROUGH A FILE OF MS-LIKE
# OUTPUT, ANALYZES A SINGLE REPLICATE
# CATEGORIZES DATA BY NONSYNONYMOUS SNP DENSITY
#############################################

my $results_file = $ARGV[0] ;
my $results2_file = $ARGV[1] ;
my $selpos_file = $ARGV[2] ;
my $rep = $ARGV[3] ;
my $gene_size = 100 ;
my $quantiles = 5 ;

unless(-e "rep${rep}_summaries"){
	system("mkdir rep${rep}_summaries") ;
}
open QC, ">./rep${rep}_summaries/QC_SharedSNPs.txt" ;

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
my %Segsites_tmphash ; # $Segsites{site} ++
my %Segsites2_tmphash ; # $Segsites{site} ++
my %Segsites ; # $Segsites{site} ++
my %Segsites2 ; # $Segsites{site} ++
my %Segsites_Total_hash ; # $Segsites{site} ++
my %Segsites_Total ; # $Segsites{site} ++

my %Segsites_Bi_PosInChromo; ; # $Segsites{index} = site ;
#my %Segsites_Bi_PosInGenostring ; # for mapping back to positin in string of 0's and 1's

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
				$Segsites2_tmphash{$line[$x]}++ ; # to check if sel positions exist
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
print QC "Num segsites N1\t", scalar keys %Segsites_tmphash, "\n" ;
print QC "Num segsites N2\t", scalar keys %Segsites2_tmphash, "\n" ;

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
	if( !exists $Segsites_tmphash{$pos} && !exists $Segsites2_tmphash{$pos}){
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
foreach my $pos ( sort{$a <=> $b} keys %Segsites_tmphash ){
	$Segsites_Total_hash{$pos}++ ;
}
foreach my $pos ( sort{$a <=> $b} keys %Segsites2_tmphash ){
	$Segsites_Total_hash{$pos}++ ;
}
my $cnt2=0 ;
foreach my $pos ( sort{$a <=> $b} keys %Segsites_Total_hash ){
	$Segsites_Total{$cnt2} = int($pos*$num_sites) ;
	$cnt2++ ;
}


print QC "Segsites_Total_hash b4 filter: ", scalar keys %Segsites_Total_hash, "\n" ;
print QC "Segsites_Total b4 filter: ", scalar keys %Segsites_Total, "\n" ;
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
print QC "Segsites_Total aft filter: ", scalar keys %Segsites_Total, "\n" ;



my $count = 0 ;
# since these are all positions in sample, first is pos 0 in genostring, second is pos 1, etc...
foreach my $pos ( sort{$a <=> $b} keys %Segsites_tmphash ){
	my $chromopos = int($pos*$num_sites) ;
	if(!exists $SitesToExclude{$chromopos} ){
		$Segsites{$chromopos}++ ;
    	$Segsites_Bi_PosInChromo{$count} = $chromopos ; # position in chromosome
	}
	# always increment counter, since this is position in genostring, and we need to account
	# for sites skipped b/c multiallelic
	$count++ ;
}


foreach my $pos ( sort{$a <=> $b} keys %Segsites2_tmphash ){
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

print QC "segsites b4 and after filtering \n" ;
print QC scalar keys %Segsites_tmphash, "\t", scalar keys %Segsites, "\n" ;
print QC scalar keys %Segsites2_tmphash, "\t", scalar keys %Segsites2, "\n" ;

#############################################
##	MAKE FRAGMENT INTO GENES
#############################################
my %Functional_Effect_BiAllelic ;
my %GenostringSite_2_Gene ;

foreach my $count ( sort{$a<=>$b} keys %Segsites_Bi_PosInChromo ){
	my $gene = int($Segsites_Bi_PosInChromo{$count}/$gene_size) ; # genes indexed from 0 to ...
	my $chromo_site = $Segsites_Bi_PosInChromo{$count} ;
	my $genostring_site = $count ;

	$GenostringSite_2_Gene{$genostring_site} = $gene ;
	if( exists $SelPos{$chromo_site} ){
		$Functional_Effect_BiAllelic{$gene}{$genostring_site} = "NONSYNONYMOUS" ;
	}else{
		$Functional_Effect_BiAllelic{$gene}{$genostring_site} = "SYNONYMOUS" ;
	}
}

#############################################
##	QUANTILE GENES BY SHARED SNP DENSITY
#############################################
my %SHARED_SNP_densities ; # $SNP_densities{$gene} = density
my %Genes_Quantiled ; # $Genes_Quantiled{gene} = quantile
my %Quant_means ;
foreach my $gene( sort {$a <=> $b} keys %Functional_Effect_BiAllelic ){
	$SHARED_SNP_densities{$gene} = 0 ;
	foreach my $genostring_site (sort {$a <=> $b} keys %{$Functional_Effect_BiAllelic{$gene}} ){
		
		my $chromopos = $Segsites_Bi_PosInChromo{$genostring_site} ;
		if(exists $Segsites2{$chromopos} ){
			$SHARED_SNP_densities{$gene} += 1/$gene_size ;
		}
	}
}

my $counter = 0 ;
my $numGenesPerQuantile = (scalar keys %SHARED_SNP_densities)/$quantiles ;
print QC "numGenesPerQuantile", "\t", $numGenesPerQuantile, "\n" ;
foreach my $gene ( sort{ $SHARED_SNP_densities{$a} <=> $SHARED_SNP_densities{$b} } keys %SHARED_SNP_densities){
	my $quant = int($counter/$numGenesPerQuantile) ;
	$Genes_Quantiled{$gene} = $quant ;
	$Quant_means{$quant} += $SHARED_SNP_densities{$gene}/$numGenesPerQuantile ;
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
foreach my $cnt ( sort{$a <=> $b} keys %Segsites_Bi_PosInChromo ){
	my $chromo_site = $Segsites_Bi_PosInChromo{$cnt} ;
	my $genostring_site = $cnt ;
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























=begin
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
=cut

exit ;


