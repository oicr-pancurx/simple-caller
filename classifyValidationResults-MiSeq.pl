#!/usr/bin/perl

use strict;
use warnings;


my $reference_coverage_minimum = 15;
my $primary_coverage_minimum = 15;

my $reference_somatic_freq_max = 0.01;
my $primary_somatic_freq_min = 0.02;
my $primary_variant_reads_min = 3;

my $reference_germline_freq_min = 0.10;
my $primary_germline_freq_min = 0.10;

my $reference_wildtype_freq_max = 0.02;
my $primary_wildtype_freq_max = 0.02;

my $indel_homopolymer_freq_max = 0.10;
my $doHomopolymerFilter = 0;


my $normalFile = $ARGV[0];
my $tumourFile = $ARGV[1];


my %normalCalls;
my %tumourCalls;

my $l;
my @f;

my ($chr, $start, $end, $pos, $ref, $alt, $freq, $depth);

open (NORM, "<$normalFile") or die "Couldn't open [$normalFile]\n";

while ($l = <NORM>)
{
	chomp $l;
	unless ($l =~ /^#/)
	{
		@f = split(/\t/, $l);

		$chr = $f[0];
		$pos = $f[1];
		$ref = $f[3];
		$alt = $f[4];

		if ($f[7] =~ /AF=(.*?);DP=(.*?);/)
		{
			$freq = $1;
			$depth = $2;
		}
		else
		{
			die "Couldn't parse $f[7]\n";
		}

		$normalCalls{"$chr\t$pos"}{depth} = $depth;
		$normalCalls{"$chr\t$pos"}{"$ref\t$alt"}{freq} = $freq;
	}
}
close NORM;

open (TUMOUR, "<$tumourFile") or die "Couldn't open [$tumourFile]\n";
while ($l = <TUMOUR>)
{
	chomp $l;
	unless ($l =~ /^#/)
	{
		@f = split(/\t/, $l);

		$chr = $f[0];
		$pos = $f[1];
		$ref = $f[3];
		$alt = $f[4];

		if ($f[7] =~ /AF=(.*?);DP=(.*?);/)
		{
			$freq = $1;
			$depth = $2;
		}
		else
		{
			die "Couldn't parse $f[7]\n";
		}

		$tumourCalls{"$chr\t$pos"}{depth} = $depth;
		$tumourCalls{"$chr\t$pos"}{"$ref\t$alt"}{freq} = $freq;
	}
}
close TUMOUR;



my @alts;

my $total = 0;
my $snvs = 0;
my $indels = 0;

my $depthFailPri = 0;
my $depthFailRef = 0;
my $depthFail = 0;

my $somatic = 0;
my $germline = 0;
my $wildtype = 0;
my $homopolymer = 0;
my $unknown = 0;

my $chrPos;
my $change;
my $outcome;

my $refAlt;
my $insFreq;
my $delFreq;
my $indelRef;
my $indelAlt;
my @indelSplit;

my ($name,$pDepth,$pFreq,$rDepth,$rFreq);
my $normalAltDepth;
my $tumourAltDepth;
my $info;

while ($l = <STDIN>)
{
	chomp $l;

	if ($l =~ /^#/)
	{
		print "$l\n";
	}
	else
	{

		@f = split(/\t/, $l);

		$chr = $f[0];
		$pos = $f[1];
		$ref = $f[3];
		$info = $f[7];

		$rDepth = 0;
		$pDepth = 0;
		$rFreq = 0;
		$pFreq = 0;

		@alts = split(/,/, $f[4]);
		for $alt (@alts)
		{
			$total++;
	
			if (length($ref) == length($alt))
			{
				$snvs++;
			}
			else
			{
				$indels++;
			}
	
			$chrPos = "$chr\t$pos";
			$change = "$ref\t$alt";
	
			$outcome = "";
	
	
			if ((exists $tumourCalls{$chrPos}) and (exists $normalCalls{$chrPos}))
			{
				$pDepth = $tumourCalls{$chrPos}{depth};
				$rDepth = $normalCalls{$chrPos}{depth};
	
				if (exists $tumourCalls{$chrPos}{$change}{freq})
				{
					$pFreq = $tumourCalls{$chrPos}{$change}{freq};
				}
				else
				{
					$pFreq = 0;
				}
		
				if (exists $normalCalls{$chrPos}{$change}{freq})
				{
					$rFreq = $normalCalls{$chrPos}{$change}{freq};
				}
				else
				{
					$rFreq = 0;
				}
		
		
				if ($pDepth < $primary_coverage_minimum)
				{
					$depthFailPri++;
					$depthFail++;
					$outcome = "INSUFFICIENT_COVERAGE";
				}
				elsif ($rDepth < $reference_coverage_minimum)
				{
					$depthFailRef++;
					$depthFail++;
					$outcome = "INSUFFICIENT_COVERAGE";
				}
	
	
				# passed QC
				elsif (
					($rFreq < $reference_somatic_freq_max)
						and
					($pFreq >= $primary_somatic_freq_min)
						and
					(($pDepth * $pFreq) >= $primary_variant_reads_min)
				)
				{
					$somatic++;
					$outcome = "TRUE-SOMATIC";
				}
			
				elsif (
					($rFreq >= $reference_germline_freq_min)
						and
					($pFreq >= $primary_germline_freq_min)
				)
				{
					$germline++;
					$outcome = "FALSE-GERMLINE";
				}
			
				elsif (
					($rFreq < $reference_wildtype_freq_max)
						and
					($pFreq < $primary_wildtype_freq_max)
				)
				{
					$wildtype++;
					$outcome = "FALSE-WILDTYPE";
				}
			
				else		# did not match anything, so classify unknown
				{
					$unknown++;
					$outcome = "UNKNOWN";
				}
	
				if ($doHomopolymerFilter == 1)
				{
					# check indels for evidence of homopolymer noise
					unless (($outcome eq "TRUE-SOMATIC") or ($outcome eq "INSUFFICIENT_COVERAGE"))
					{
						if (length($ref) != length($alt))
						{
							$insFreq = 0;
							$delFreq = 0;
	
							for my $refAlt (keys %{ $normalCalls{$chrPos} })
							{
								unless ($refAlt eq "depth")
								{
									@indelSplit = split(/\t/, $refAlt);
									$indelRef = $indelSplit[0];
									$indelAlt = $indelSplit[1];
									
									if (length($indelRef) > length($indelAlt))
									{
										$delFreq += $normalCalls{$chrPos}{$refAlt}{freq};
									}
									elsif (length($indelRef) < length($indelAlt))
									{
										$insFreq += $normalCalls{$chrPos}{$refAlt}{freq};
									}
								
								}
							}
	
							if (($delFreq > $indel_homopolymer_freq_max) and ($insFreq > $indel_homopolymer_freq_max))
							{
								if ($outcome eq "FALSE-WILDTYPE")
								{
									$wildtype--;
								}
								elsif ($outcome eq "FALSE-GERMLINE")
								{
									$germline--;
								}
								elsif ($outcome eq "UNKNOWN")
								{
									$unknown--;
								}
	
								$homopolymer++;
								$outcome = "NOISY_INDEL";
							}
						}
					}
				}
			}
			else
			{
				# insufficient coverage
				$outcome = "INSUFFICIENT_COVERAGE";
				$depthFail++;

	
				unless (exists $tumourCalls{$chrPos})
				{
					$depthFailPri++;
				}
				unless (exists $normalCalls{$chrPos})
				{
					$depthFailRef++;
				}
			}

			$normalAltDepth = int($rFreq*$rDepth);
			$tumourAltDepth = int($pFreq*$pDepth);

			$f[4] = $alt;
			$f[7] = "$info;VERIFICATION=$outcome;VN_ADP,TDP=$normalAltDepth,$rDepth;VT_ADP,TDP=$tumourAltDepth,$pDepth";
	
			print $f[0];
			for (my $i = 1; $i < scalar(@f); $i++)
			{
				print "\t$f[$i]";
			}
			print "\n";
		}

	}
}


warn "\nTotal Variants:             $total\n";
warn "SNVs:                       $snvs\n";
warn "Indels:                     $indels\n";

warn "\nInsufficient tumour depth:  $depthFailPri\n";
warn "Insufficient normal depth:  $depthFailRef\n\n";

warn "Somatic:     $somatic\n";
warn "Germline:    $germline\n";
warn "Wildtype:    $wildtype\n\n";

warn "Homopolymer: $homopolymer\n";
warn "Unknown:     $unknown\n";
warn "Depth fail:  $depthFail\n\n";













