#!usr/bin/perl
use warnings ;
use strict ;

my $file = $ARGV[0] ;

my @Thetas ;
my @Ds ;

open IN, "<${file}" ;
while(<IN>){
	chomp $_ ;
	my @line = split(/\t/, $_) ;
	if($_ =~ m/FiniteSites_Theta/){
		push @Thetas, $line[2] ;
	}elsif($_ =~ m/Average_D/){
		push @Ds, $line[2] ;
	}
}
close IN ;

@Ds = sort{$a <=> $b} @Ds ;
@Thetas = sort{$a <=> $b} @Thetas ;

my $quant1 = int(0.25*(scalar @Ds) - 1) ;
my $quant2 = int(0.5*(scalar @Ds) - 1) ;
my $quant3 = int(0.75*(scalar @Ds) - 1) ;

print "Thetas:\t" ;
print $Thetas[0], "\t" ;
print $Thetas[$quant1], "\t" ;
print $Thetas[$quant2], "\t" ;
print $Thetas[$quant3], "\t" ;
print $Thetas[$#Thetas], "\n" ;

print "Ds:\t" ;
print $Ds[0], "\t" ;
print $Ds[$quant1], "\t" ;
print $Ds[$quant2], "\t" ;
print $Ds[$quant3], "\t" ;
print $Ds[$#Ds], "\n" ;


exit ;



