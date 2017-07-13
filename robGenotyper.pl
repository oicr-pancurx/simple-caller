#!/usr/bin/perl

use strict;
use warnings;

my $bedFile = $ARGV[0];

my $l;
my @f;

my $homoDepth = 8;
my $homoFreq = 0.02;

my $hetDepth = 4;
my $hetBases = 2;
my $hetFreq = 0.2;

my %positionsOfInterest;
my @posOrder;

my $pos;
my @bases = qw(A C G T);
my $call;

open (BED, "<$bedFile") or die "Couldn't open $bedFile\n";
while ($l = <BED>)
{
	unless ($l =~ /^#/)
	{
		chomp $l;
		@f = split(/\t/, $l);

		$pos = "$f[0]\t$f[1]";

		if (exists $f[4])		# vcf, don't want indels
		{
			if (length($f[3]) == length($f[4]))
			{
				$positionsOfInterest{$pos}{label} = "$f[3]->$f[4]";
				$positionsOfInterest{$pos}{call} = "";
				push(@posOrder, $pos);
			}
		}
		else
		{
			$positionsOfInterest{$pos}{label} = $f[3];
			$positionsOfInterest{$pos}{call} = "";
			push(@posOrder, $pos);
		}
	}
}
close BED;




while ($l = <STDIN>)		# assuming no indels
{
	chomp $l;
	@f = split(/\t/, $l);
	
	$pos = "$f[0]\t$f[1]";
	$positionsOfInterest{$pos}{call} = "";

	if (exists $positionsOfInterest{$pos})
	{
		$positionsOfInterest{$pos}{depth} = $f[6];

		$positionsOfInterest{$pos}{A} = $f[7];
		$positionsOfInterest{$pos}{C} = $f[8];
		$positionsOfInterest{$pos}{G} = $f[9];
		$positionsOfInterest{$pos}{T} = $f[10];

		# test for homozygous
		if ($positionsOfInterest{$pos}{depth} >= $homoDepth)
		{
			for my $base (@bases)
			{
				if ($positionsOfInterest{$pos}{$base} >= 1 - $homoFreq)
				{
					$positionsOfInterest{$pos}{call} = "$base$base";
				}
			}
		}


		# test for heterozygous
		if ($positionsOfInterest{$pos}{call} eq "")
		{
			if ($positionsOfInterest{$pos}{depth} >= $hetDepth)
			{
				for my $base (@bases)
				{
					if (($positionsOfInterest{$pos}{$base} >= $hetFreq) and ($positionsOfInterest{$pos}{$base} * $positionsOfInterest{$pos}{depth} >= 2))
					{
						$positionsOfInterest{$pos}{call} .= "$base";
					}
				}
				if (length($positionsOfInterest{$pos}{call}) < 2)
				{
					$positionsOfInterest{$pos}{call} = "";		# a het with only one base is no good
				}
			}
			
			unless ($positionsOfInterest{$pos}{call} eq "")
			{
				$positionsOfInterest{$pos}{call} = join('', sort { $a cmp $b } split(//, $positionsOfInterest{$pos}{call}));
			}
		}
	}

}


for $pos (@posOrder)
{
	print "$pos\t$positionsOfInterest{$pos}{label}\t$positionsOfInterest{$pos}{call}\n";
}


