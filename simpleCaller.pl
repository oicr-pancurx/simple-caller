#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Std;
use vars qw/ %opt /;
use Bio::DB::Fasta;


my %cuts = (
	qual => 0,
	depth => 100,
	freq => 0.02,

	homoDepth => 8,
	homoFreq => 0.98,

	hetDepth => 4,
	hetBases => 2,
	hetFreq => 0.2,

	outputRef => "no"
);

my %fastaHandles = ("path" => "/oicr/data/reference/genomes/homo_sapiens_mc/UCSC/hg19_random/Genomic/references/fasta/");
my $maxGeneralizeLength = 2000;

my $dumpTime = 1000;
my $dumpToFileCounter = 0;

my %refHash = ("chunkSize" => 10000);

sub usage
{
	print "\nUsage is samtools view -hF 4 file.bam | simpleCaller.pl [options].\n";
	print "Options are as follows:\n";
	print "\t-q minimum quality for a base to be considered.  Default is $cuts{qual}.\n";
	print "\t-d minimum depth of coverage for reporting a variant.  Default is $cuts{depth}.\n";
	print "\t-f minumum minor allele frequency for reporting.  Default is $cuts{freq}.\n";
	print "\t-A flag to turn on outputting reference base depth.  Forces -f 0 and -d 0 (for verification of somatic results).\n";
	print "\t-F path to fasta files (one per chromosome, for now).  Default is $fastaHandles{\"path\"}.\n";
	print "\t-M maximum generalize length.  Default is $maxGeneralizeLength.\n";
	print "\t-h displays this helpful usage.\n";

	die "\n@_\n\n";
}


my $optString = "q:d:f:AF:h";
getopts ($optString, \%opt) or usage("Incorrect arguments.");

if (exists $opt{h})
{
	usage("Help requested.");
}

if (exists $opt{q})
{
	$cuts{qual} = $opt{q};
}

if (exists $opt{d})
{
	$cuts{depth} = $opt{d};
}
if (exists $opt{f})
{
	$cuts{freq} = $opt{f};
}

if (exists $opt{A})
{
	$cuts{outputRef} = "yes";
	$cuts{depth} = 0;
	$cuts{freq} = 0;
}

if (exists $opt{F})
{
	$fastaHandles{"path"} = $opt{F};
}
if (exists $opt{M})
{
	$maxGeneralizeLength = $opt{M};
}


## print VCF header


my ($l, @f);

my %snvHash;		# $snvHash{chr\tpos}{alt}{count} = count, {strand} = ++ or --, {cycle_sum} = ++  | cycle - readLength/2 | {qual_sum}
my %insHash;		# $insHash{chr\tpos}{alt}{count} = count, $insHash{chr\tpos}{alt}{isGeneralized} = bool
my %delHash;		# $delHash{chr\tpos}{alt}{count} = count, $delHash{chr\tpos}{alt}{isGeneralized} = bool
my %allHash;		# $depthHash{chr\tpos}{count}++;
my %printBuffer;

my ($flags, $chr, $startPos, $cigar, @bases, @quals, $base, $alt, $q, $qSum, $mapq, $pos, $anchorPos, $readLength);

my $isForward;
my @cigarOp;

my $doPrintHeader = 0;
my @chrOrder;
my @readGroup;

my %samples;
my %libraries;
my %lanes;
my $rgid;

while ($l = <STDIN>)
{
	chomp $l;
	if ($l =~ /^@/)		# header
	{
		if ($l =~ /^\@SQ\tSN:(.*?)\t/)
		{
			push(@chrOrder, $1);
		}
		elsif ($l =~ /^\@RG\t(.*)/)
		{
			@readGroup = split(/\t/, $1);

			for my $rg (@readGroup)
			{
				if ($rg =~ /^SM:(.*)/)
				{
					$samples{$1}++;
				}
				elsif ($rg =~ /^LB:(.*)/)
				{
					$libraries{$1}++;
				}
				elsif ($rg =~ /^ID:(.*)/)
				{
					$rgid = $1;
					if ($rgid =~ /SWID/)
					{
						for my $rg2 (@readGroup)
						{
							if ($rg2 =~ /^PU:(.*)/)
							{
								$lanes{"${rgid}_$1"}++;
							}
						}
					}
					else
					{
						$lanes{$rgid}++;
					}
				}
			}
		}
	}
	else
	{
		if ($doPrintHeader == 0)
		{
			$doPrintHeader = 1;
			printHeader(\@chrOrder, \%samples, \%libraries, \%lanes, \%cuts);
		}

		@f = split(/\t/, $l);

		$flags = $f[1];
		$chr = $f[2];
		$startPos = $f[3];
		$mapq = $f[4];
		$cigar = $f[5];
		@bases = split(//, $f[9]);
		@quals = split(//, $f[10]);

		$readLength = scalar(@bases);

		$isForward = procFlags($flags);
		@cigarOp = procCigar($cigar);
		$pos = $startPos;

		foreach my $c (@cigarOp)
		{
			if ($c =~ /(.*)M/)
			{
				for (my $i = 0; $i < $1; $i++)
				{
					$base = shift(@bases);
					$q = toPhred(shift(@quals));

					if ($q >= $cuts{qual})
					{
						$snvHash{$chr}{$pos}{$base}{count}++;
						$snvHash{$chr}{$pos}{$base}{qual_sum} += $q;
						$snvHash{$chr}{$pos}{$base}{mapq_sum} += $mapq;
						$snvHash{$chr}{$pos}{$base}{cycle_sum} += 1 - abs($readLength/2 - ($pos - $startPos)) / ($readLength/2);
						if ($isForward == 1)
						{
							$snvHash{$chr}{$pos}{$base}{strand}++;
						}
						else
						{
							$snvHash{$chr}{$pos}{$base}{strand}--;
						}
					}


					if ($isForward == 1)
					{
						$allHash{$chr}{$pos}{strand}++;
					}
					else
					{
						$allHash{$chr}{$pos}{strand}--;
					}
					$allHash{$chr}{$pos}{depth}++;

					$pos++;
				}
			}
			elsif ($c =~ /(.*)S/)
			{
				for (my $i = 0; $i < $1; $i++)
				{
					$base = shift(@bases);
					$q = shift(@quals);
				}
			}
			elsif ($c =~ /(.*)I/)
			{
				$anchorPos = $pos - 1;
				$alt = getBase($chr, $anchorPos, \%refHash, \%fastaHandles);
				$qSum = 0;


				for (my $i = 0; $i < $1; $i++)
				{
					$base = shift(@bases);
					$q = toPhred(shift(@quals));
					$alt .= $base;
					$qSum += $q;
				}


				$insHash{$chr}{$anchorPos}{$alt}{count}++;
				$insHash{$chr}{$anchorPos}{$alt}{qual_sum} += $qSum;
				$insHash{$chr}{$anchorPos}{$alt}{mapq_sum} += $mapq;
				$insHash{$chr}{$anchorPos}{$alt}{cycle_sum} += 1 - abs($readLength/2 - ($anchorPos + (length($alt) / 2) - $startPos)) / ($readLength/2);
				if ($isForward == 1)
				{
					$insHash{$chr}{$anchorPos}{$alt}{strand}++;
				}
				else
				{
					$insHash{$chr}{$anchorPos}{$alt}{strand}--;
				}

			}
			elsif ($c =~ /(.*)D/)
			{
				$anchorPos = $pos - 1;
				$alt = getBase($chr, $anchorPos, \%refHash, \%fastaHandles);
				$qSum = 0;

				for (my $i = 0; $i < $1; $i++)
				{
					$alt .= getBase($chr, $pos, \%refHash, \%fastaHandles);

					if ($isForward == 1)
					{
						$allHash{$chr}{$pos}{strand}++;
					}
					else
					{
						$allHash{$chr}{$pos}{strand}--;
					}
					$allHash{$chr}{$pos}{depth}++;

					$snvHash{$chr}{$pos}{del}{count}++;

					$pos++;
				}

				$delHash{$chr}{$anchorPos}{$alt}{count}++;
				$delHash{$chr}{$anchorPos}{$alt}{mapq_sum} += $mapq;
				$delHash{$chr}{$anchorPos}{$alt}{cycle_sum} += 1 - abs($readLength/2 - ($anchorPos + (length($alt) / 2) - $startPos)) / ($readLength/2);
				if ($isForward == 1)
				{
					$delHash{$chr}{$anchorPos}{$alt}{strand}++;
				}
				else
				{
					$delHash{$chr}{$anchorPos}{$alt}{strand}--;
				}

			}
			elsif ($c =~ /(.*)H/)
			{
			}
			else
			{
				die "Can't handle CIGAR operation: $c\n";
			}
		}


		$dumpToFileCounter++;
		if ($dumpToFileCounter >= $dumpTime)
		{
			processVariants($chr, $pos, $maxGeneralizeLength, \%snvHash, \%insHash, \%delHash, \%allHash, \%cuts, \%refHash, \%fastaHandles, \%printBuffer);
			sortAndPrintBuffer(\%printBuffer);
			freeMemChunks($chr, $pos, \%refHash);
			$dumpToFileCounter = 0;
		}
	}

}

# print remaining variants
processVariants("do the rest", $pos + 1, 0, \%snvHash, \%insHash, \%delHash, \%allHash, \%cuts, \%refHash, \%fastaHandles, \%printBuffer);
sortAndPrintBuffer(\%printBuffer);




sub printHeader
{
	my $chrOrder = $_[0];
	my $samples = $_[1];
	my $libraries = $_[2];
	my $lanes = $_[3];
	my $cuts = $_[4];

	print "##fileformat=VCFv4.1\n";
	print "##source=simpleCaller.pl\n";

	for my $cut (sort keys %{ $cuts })
	{
		print "##cutoff=<ID=$cut,VALUE=$cuts->{$cut}>\n";
	}

	for my $chr (@{ $chrOrder })
	{
		print "##contig=<ID=$chr>\n";
	}

	for my $lane (sort keys %{ $lanes })
	{
		print "##lane=<ID=$lane>\n";
	}

	for my $lib (sort keys %{ $libraries })
	{
		print "##library=<ID=$lib,LANES=$libraries->{$lib}>\n";
	}

	my $sampCount = 0;
	my $sample = "";
	for my $samp (sort keys %{ $samples })
	{
		print "##sample=<ID=$samp>\n";
		$sampCount++;
		$sample .= "$samp,";
	}
	if ($sampCount > 1)
	{
		warn "Multiple samples in the bam \@RG info!\n";
	}

	$sample =~ s/,$//;


	print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sample\n";
}


sub processVariants
{
	my $printChr = $_[0];
	my $pos = $_[1];
	my $maxGeneralizeLength = $_[2];
	my $snvHash = $_[3];
	my $insHash = $_[4];
	my $delHash = $_[5];
	my $allHash = $_[6];
	my $cuts = $_[7];
	my $refHash = $_[8];
	my $fastaHandles = $_[9];
	my $printBuffer = $_[10];

	my $printPos = $pos - $maxGeneralizeLength;
	my $generalizePos = $pos - ($maxGeneralizeLength / 2);

	my $snvDepth;
	my $newPos;

	# generalize first and adjust ref and del counts in the snv hash when deletions move
	for my $chr (keys %{ $insHash })
	{
		if ($chr ne $printChr)
		{
			for my $pos (keys %{ $insHash->{$chr} })
			{
				for my $alt (keys %{ $insHash->{$chr}{$pos} })
				{
					unless (exists $insHash->{$chr}{$pos}{$alt}{isGeneralized})
					{
						$newPos = generalizeIndel("+", $chr, $pos, $alt, $insHash, $refHash, $fastaHandles);
					}
				}
			}
		}
		else
		{
			for my $pos (keys %{ $insHash->{$chr} })
			{
				if ($pos < $generalizePos)
				{
					for my $alt (keys %{ $insHash->{$chr}{$pos} })
					{
						unless (exists $insHash->{$chr}{$pos}{$alt}{isGeneralized})
						{
							$newPos = generalizeIndel("+", $chr, $pos, $alt, $insHash, $refHash, $fastaHandles);
						}
					}
				}
			}
		}
	}



	my $delCount;
	for my $chr (keys %{ $delHash })
	{
		if ($chr ne $printChr)
		{
			for my $pos (keys %{ $delHash->{$chr} })
			{
				for my $alt (keys %{ $delHash->{$chr}{$pos} })
				{
					unless (exists $delHash->{$chr}{$pos}{$alt}{isGeneralized})
					{
						$delCount = $delHash->{$chr}{$pos}{$alt}{count};
						$snvHash->{$chr}{$pos}{del}{count} -= $delCount;

						$newPos = generalizeIndel("-", $chr, $pos, $alt, $delHash, $refHash, $fastaHandles);

						$snvHash->{$chr}{$newPos}{del}{count} += $delCount;
					}
				}
			}
		}
		else
		{
			for my $pos (keys %{ $delHash->{$chr} })
			{
				if ($pos < $generalizePos)
				{
					for my $alt (keys %{ $delHash->{$chr}{$pos} })
					{
						unless (exists $delHash->{$chr}{$pos}{$alt}{isGeneralized})
						{
							$delCount = $delHash->{$chr}{$pos}{$alt}{count};
							$snvHash->{$chr}{$pos}{del}{count} -= $delCount;

							$newPos = generalizeIndel("-", $chr, $pos, $alt, $delHash, $refHash, $fastaHandles);

							$snvHash->{$chr}{$newPos}{del}{count} += $delCount;
						}
					}
				}
			}
		}
	}




	for my $chr (keys %{ $snvHash })
	{
		if ($chr ne $printChr)
		{
			for my $pos (keys %{ $snvHash->{$chr} })
			{
				for my $alt (keys %{ $snvHash->{$chr}{$pos} })
				{
					if ((($alt ne getBase($chr, $pos, $refHash, $fastaHandles)) or $cuts->{outputRef} eq "yes") and ($alt ne "del") and ($alt ne "N"))
					{
						$snvDepth = 0;
						for my $base (qw/A C G T N del/)
						{
							unless (exists $snvHash->{$chr}{$pos}{$base})
							{
								$snvHash->{$chr}{$pos}{$base}{count} = 0;
							}
							$snvDepth += $snvHash->{$chr}{$pos}{$base}{count};
						}

						if (($snvDepth >= $cuts->{depth}) and ($snvDepth > 0))
						{
							if ($snvHash->{$chr}{$pos}{$alt}{count} / $snvDepth >= $cuts->{freq})
							{
								printSNV($chr, $pos, $alt, $snvHash, $refHash, $fastaHandles, $cuts, $printBuffer);
							}
						}
					}
				}
			}
			delete $snvHash->{$chr};
		}
		else
		{
			for my $pos (keys %{ $snvHash->{$chr} })
			{
				if ($pos < $printPos)
				{
					for my $alt (keys %{ $snvHash->{$chr}{$pos} })
					{
						if ((($alt ne getBase($chr, $pos, $refHash, $fastaHandles)) or $cuts->{outputRef} eq "yes") and ($alt ne "del") and ($alt ne "N"))
						{
							$snvDepth = 0;
							for my $base (qw/A C G T N del/)
							{
								unless (exists $snvHash->{$chr}{$pos}{$base})
								{
									$snvHash->{$chr}{$pos}{$base}{count} = 0;
								}
								$snvDepth += $snvHash->{$chr}{$pos}{$base}{count};
							}
							if (($snvDepth >= $cuts->{depth}) and ($snvDepth > 0))
							{
								if ($snvHash->{$chr}{$pos}{$alt}{count} / $snvDepth >= $cuts->{freq})
								{
									printSNV($chr, $pos, $alt, $snvHash, $refHash, $fastaHandles, $cuts, $printBuffer);
								}
							}
						}
					}
					delete $snvHash->{$chr}{$pos};
				}
			}
		}
	}


	for my $chr (keys %{ $insHash })
	{
		if ($chr ne $printChr)
		{
			for my $pos (keys %{ $insHash->{$chr} })
			{
				unless (exists $allHash{$chr}{$pos}{depth})
				{
					$allHash{$chr}{$pos}{depth} = 0;
				}
				if ($allHash{$chr}{$pos}{depth} >= $cuts->{depth})
				{
					for my $alt (keys %{ $insHash->{$chr}{$pos} })
					{
						unless ($allHash{$chr}{$pos}{depth} == 0)		# if depth is zero, insertion is at the start of its reads (which is dumb) so skip
						{
							if (($insHash->{$chr}{$pos}{$alt}{count} / $allHash{$chr}{$pos}{depth}) >= $cuts->{freq})
							{
								printIns($chr, $pos, $alt, $insHash, $allHash, $refHash, $fastaHandles, $printBuffer, $cuts);
							}
						}
					}
				}
			}
			delete $insHash->{$chr};
		}
		else
		{
			for my $pos (keys %{ $insHash->{$chr} })
			{
				if ($pos < $printPos)
				{
					unless (exists $allHash{$chr}{$pos}{depth})
					{
						$allHash{$chr}{$pos}{depth} = 0;
					}
					if ($allHash{$chr}{$pos}{depth} >= $cuts->{depth})
					{
						for my $alt (keys %{ $insHash->{$chr}{$pos} })
						{
							unless ($allHash{$chr}{$pos}{depth} == 0)		# if depth is zero, insertion is at the start of its reads (which is dumb) so skip
							{
								if (($insHash->{$chr}{$pos}{$alt}{count} / $allHash{$chr}{$pos}{depth}) >= $cuts->{freq})
								{
									printIns($chr, $pos, $alt, $insHash, $allHash, $refHash, $fastaHandles, $printBuffer, $cuts);
								}
							}
						}
					}
					delete $insHash->{$chr}{$pos};
				}
			}
		}
	}

	for my $chr (keys %{ $delHash })
	{
		if ($chr ne $printChr)
		{
			for my $pos (keys %{ $delHash->{$chr} })
			{
				if ($allHash{$chr}{$pos}{depth} >= $cuts->{depth})
				{
					for my $alt (keys %{ $delHash->{$chr}{$pos} })
					{
						if (($delHash->{$chr}{$pos}{$alt}{count} / $allHash{$chr}{$pos}{depth}) >= $cuts->{freq})
						{
							printDel($chr, $pos, $alt, $delHash, $allHash, $refHash, $fastaHandles, $printBuffer, $cuts);
						}
					}
				}
			}
			delete $delHash->{$chr};
		}
		else
		{
			for my $pos (keys %{ $delHash->{$chr} })
			{
				if ($pos < $printPos)
				{
					if (exists $allHash{$chr}{$pos}{depth})
					{
						if ($allHash{$chr}{$pos}{depth} >= $cuts->{depth})
						{
							for my $alt (keys %{ $delHash->{$chr}{$pos} })
							{
								if (($delHash->{$chr}{$pos}{$alt}{count} / $allHash{$chr}{$pos}{depth}) >= $cuts->{freq})
								{
									printDel($chr, $pos, $alt, $delHash, $allHash, $refHash, $fastaHandles, $printBuffer, $cuts);
								}
							}
						}
					}
					delete $delHash->{$chr}{$pos};
				}
			}
		}
	}

	for my $chr (keys %{ $allHash })
	{
		if ($chr ne $printChr)
		{
			delete $allHash->{$chr};
		}
		else
		{
			for my $pos (keys %{ $allHash->{$chr} })
			{
				if ($pos < $printPos)
				{
					delete $allHash->{$chr}{$pos};
				}
			}
		}
	}
}


sub printSNV
{
	my $chr = $_[0];
	my $pos = $_[1];
	my $alt = $_[2];
	my $snvHash = $_[3];
	my $refHash = $_[4];
	my $fastaHandles = $_[5];
	my $cuts = $_[6];
	my $printBuffer = $_[7];

	my $ref = getBase($chr, $pos, $refHash, $fastaHandles);

	my $totalDepth = 0;
	my $totalStrand = 0;
	for my $base (qw/A C G T N del/)
	{
		unless (exists $snvHash->{$chr}{$pos}{$base}{count})
		{
			$snvHash->{$chr}{$pos}{$base}{count} = 0;
		}
		unless (exists $snvHash->{$chr}{$pos}{$base}{strand})
		{
			$snvHash->{$chr}{$pos}{$base}{strand} = 0;
		}
		$totalDepth += $snvHash->{$chr}{$pos}{$base}{count};
		$totalStrand += $snvHash->{$chr}{$pos}{$base}{strand};
	}

	my $freq = sprintf("%.4f", $snvHash->{$chr}{$pos}{$alt}{count} / $totalDepth);

	my $altSB = sprintf("%.2f", $snvHash->{$chr}{$pos}{$alt}{strand} / $snvHash->{$chr}{$pos}{$alt}{count});
	my $refSB = 0;
	if ($snvHash->{$chr}{$pos}{$ref}{count} > 0)
	{
		$refSB = sprintf("%.2f", $snvHash->{$chr}{$pos}{$ref}{strand} / $snvHash->{$chr}{$pos}{$ref}{count});
	}
	my $totSB = sprintf("%.2f", $totalStrand / $totalDepth);


	my $info = "AF=$freq;DP=$totalDepth;ALT_SB=$altSB;REF_SB=$refSB;TOT_SB=$totSB";

	my $format = "GT:A:C:G:T:N:del";

	my $genotype;

	if ($ref ne $alt)
	{
		$genotype = getGenotype($totalDepth, $snvHash->{$chr}{$pos}{$ref}{count}, $snvHash->{$chr}{$pos}{$alt}{count}, $cuts);
	}
	else
	{
		$genotype = getGenotype($totalDepth, $snvHash->{$chr}{$pos}{$ref}{count}, 0, $cuts);
	}

	for my $base (qw/A C G T N del/)
	{
		$genotype .= ":$snvHash->{$chr}{$pos}{$base}{count}";
	}


	if ($ref ne $alt)
	{
		$printBuffer->{$chr}{$pos}{".\t$ref\t$alt"} = ".\t.\t$info\t$format\t$genotype";
	}
	else
	{
		$printBuffer->{$chr}{$pos}{".\t$ref\t."} = ".\t.\t$info\t$format\t$genotype";
	}
}


sub printIns
{
	my $chr = $_[0];
	my $pos = $_[1];
	my $alt = $_[2];
	my $insHash = $_[3];
	my $allHash = $_[4];
	my $refHash = $_[5];
	my $fastaHandles = $_[6];
	my $printBuffer = $_[7];
	my $cuts = $_[8];

	my $ref = getBase($chr, $pos, $refHash, $fastaHandles);

	my $freq = sprintf("%.4f", $insHash->{$chr}{$pos}{$alt}{count} / $allHash->{$chr}{$pos}{depth});

	my $altSB = sprintf("%.2f", $insHash->{$chr}{$pos}{$alt}{strand} / $insHash->{$chr}{$pos}{$alt}{count});
	my $totSB = sprintf("%.2f", $allHash->{$chr}{$pos}{strand} / $allHash->{$chr}{$pos}{depth});

	my $info = "AF=$freq;DP=$allHash->{$chr}{$pos}{depth};ALT_SB=$altSB;TOT_SB=$totSB";
	my $format = "GT:ins:other";
	my $genotype = getGenotype($allHash->{$chr}{$pos}{depth}, $allHash->{$chr}{$pos}{depth} - $insHash->{$chr}{$pos}{$alt}{count}, $insHash->{$chr}{$pos}{$alt}{count}, $cuts);
	$genotype .= ":$insHash->{$chr}{$pos}{$alt}{count}:";
	$genotype .= $allHash->{$chr}{$pos}{depth} - $insHash->{$chr}{$pos}{$alt}{count};

	$printBuffer->{$chr}{$pos}{".\t$ref\t$alt"} = ".\t.\t$info\t$format\t$genotype";

}



sub printDel
{
	my $chr = $_[0];
	my $pos = $_[1];
	my $ref = $_[2];
	my $delHash = $_[3];
	my $allHash = $_[4];
	my $refHash = $_[5];
	my $fastaHandles = $_[6];
	my $printBuffer = $_[7];
	my $cuts = $_[8];

	my $alt = getBase($chr, $pos, $refHash, $fastaHandles);

	my $freq = sprintf("%.4f", $delHash->{$chr}{$pos}{$ref}{count} / $allHash->{$chr}{$pos}{depth});

	my $altSB = sprintf("%.2f", $delHash->{$chr}{$pos}{$ref}{strand} / $delHash->{$chr}{$pos}{$ref}{count});
	my $totSB = sprintf("%.2f", $allHash->{$chr}{$pos}{strand} / $allHash->{$chr}{$pos}{depth});

	my $info = "AF=$freq;DP=$allHash->{$chr}{$pos}{depth};ALT_SB=$altSB;TOT_SB=$totSB";
	my $format = "GT:del:other";
	my $genotype = getGenotype($allHash->{$chr}{$pos}{depth}, $allHash->{$chr}{$pos}{depth} - $delHash->{$chr}{$pos}{$ref}{count}, $delHash->{$chr}{$pos}{$ref}{count}, $cuts);
	$genotype .= ":$delHash->{$chr}{$pos}{$ref}{count}:";
	$genotype .= $allHash->{$chr}{$pos}{depth} - $delHash->{$chr}{$pos}{$ref}{count};
	$printBuffer->{$chr}{$pos}{".\t$ref\t$alt"} = ".\t.\t$info\t$format\t$genotype";

}



sub sortAndPrintBuffer
{
	my $printBuffer = $_[0];

	# should improve this sort to always sort chrs in the correct order - pull the chr list from the bam header?

	for my $chr (sort keys %{ $printBuffer })
	{
		for my $pos (sort keys %{ $printBuffer->{$chr} })
		{
			for my $base (sort keys %{ $printBuffer->{$chr}{$pos} })
			{
				print "$chr\t$pos\t$base\t$printBuffer->{$chr}{$pos}{$base}\n";
			}
		}
	}
	%{ $printBuffer } = ();
}

# ("del", $chr, $pos, $alt, $delHash, $refHash, $fastaHandles)
sub generalizeIndel
{
	my $type = $_[0];
	my $chr = $_[1];
	my $oldPos = $_[2];
	my $oldAlt = $_[3];
	my $hash = $_[4];
	my $refHash = $_[5];
	my $fastaHandles = $_[6];

	my $indel = substr($oldAlt, 1);

	# find ambiguous alignments
	my $candidatePos = findCandidatePositions($type, $indel, $chr, $oldPos, $refHash, $fastaHandles);

	if (exists $candidatePos->{"deletionNotPossible"})
	{
		die "$type at $chr:$oldPos of $oldAlt? Deletion not possible omgwtfbbq!\n";
	}

	my @sortedPositions = sort { $a <=> $b } keys %$candidatePos;
	my $genPos = $sortedPositions[0];
	if ($genPos == $oldPos)	# indel was left aligned already
	{
		return $oldPos;
		$hash->{$chr}{$oldPos}{$oldAlt}{isGeneralized} = 1;
	}

	my $genAlt = getBase($chr, $genPos, $refHash, $fastaHandles);
	$genAlt .= $candidatePos->{$genPos};


	# add to generalized pos and delete old pos

	for my $metric (qw/count qual_sum mapq_sum cycle_sum strand/)
	{
		if (exists $hash->{$chr}{$genPos}{$genAlt}{$metric})
		{
			$hash->{$chr}{$genPos}{$genAlt}{$metric} += $hash->{$chr}{$oldPos}{$oldAlt}{$metric};
		}
		else
		{
			$hash->{$chr}{$genPos}{$genAlt}{$metric} = $hash->{$chr}{$oldPos}{$oldAlt}{$metric};
		}
	}

	#warn "generalized $chr:$oldPos:$oldAlt to $chr:$genPos:$genAlt\n";

	delete $hash->{$chr}{$oldPos}{$oldAlt};


	$hash->{$chr}{$genPos}{$genAlt}{isGeneralized} = 1;
	return $genPos;
}



# findCandidatePositions compares the indel sequence to the surrounding reference bases and builds a list of analogous indels
# input is the type and sequence of the indel, its chromosome and position, and references to the persistant reference hash and fasta handle hash
# output is a hash of candidate positions and resulting indel sequences
sub findCandidatePositions
{
	my $type = $_[0];
	my $indel = $_[1];
	my $chr = $_[2];
	my $pos = $_[3];
	my $reference = $_[4];
	my $fastaHandles = $_[5];

	my %candidatePos;
	my $candidateIndel;
	my $smallestInsertionSubunit;
	my $sortedIndel = join('',sort(split('',$indel)));
	my $indelLength = length($indel);

	if ($type eq "+")   # won't be able to find the insertion in reference seqeuence
	{
		$candidatePos{$pos} = $indel;
		$smallestInsertionSubunit = findSmallestInsertionSubunit($indel);
	}
	else
	{
		# check if deletion is possible
		
		$candidateIndel = getRange($chr, $pos + 1, $pos + 1 + $indelLength - 1, $reference, $fastaHandles);		# pos + 1 because pos is left of the deletion

		if ($indel eq $candidateIndel)
		{
			$candidatePos{$pos} = $indel;
		}
		else
		{
			warn "Deletion of $indel anchored at $chr:$pos is not possible (the reference is $candidateIndel), skipping generalization.\n";
			%candidatePos = ("deletionNotPossible" => 1);
			return \%candidatePos;
		}
	}

	
	# search left
	my $done = 0;
	my $searchPos = $pos;

	my $lastFoundInsertion = $pos;


	while ($done == 0)
	{
		if ($type eq "+")	# insertion
		{
			$searchPos -= length($smallestInsertionSubunit);
			$candidateIndel = getRange($chr, $searchPos + 1, $searchPos + 1 + length($smallestInsertionSubunit) - 1, $reference, $fastaHandles);  # check the bases to the right of the potential insertion for a subunit match

			if ($candidateIndel eq $smallestInsertionSubunit)
			{
				$candidatePos{$searchPos} = $indel;
			}
			else
			{
				$done = 1;
			}
		}
		else	# deletion
		{
			$candidateIndel = getRange($chr, $searchPos, $searchPos + $indelLength - 1, $reference, $fastaHandles);  # old code: substr($reference, $searchPos, $indelLength);
			if (join('',sort(split('',$candidateIndel))) eq $sortedIndel)
			{
				$candidatePos{$searchPos - 1} = $candidateIndel;	# searchPos - 1 so that we're returning the anchor base
			}
			else
			{
				$done = 1;
			}
			$searchPos--;
		}
	}


	# search right
	$searchPos = $pos;
	if ($type eq "+")
	{
		$searchPos -= length($smallestInsertionSubunit);
	}
	$done = 0;

	$lastFoundInsertion = $pos;

	while ($done == 0)
	{

		if ($type eq "+")
		{
			$searchPos += length($smallestInsertionSubunit);
			$candidateIndel = getRange($chr, $searchPos + 1, $searchPos + 1 + length($smallestInsertionSubunit) - 1, $reference, $fastaHandles);  # search the bases to the right of the position base for a subunit match

			if ($candidateIndel eq $smallestInsertionSubunit)
			{
				$candidatePos{$searchPos} = $indel;
				$candidatePos{$searchPos + length($smallestInsertionSubunit)} = $indel;     # can insert on either side of a repeated motif

			}
			else
			{
				$done = 1;
			}
		}
		else    # deletion
		{
			$searchPos++;
			$candidateIndel = getRange($chr, $searchPos, $searchPos + $indelLength - 1, $reference, $fastaHandles);  # old code: substr($reference, $searchPos, $indelLength);
			if (join('',sort(split('',$candidateIndel))) eq $sortedIndel)
			{
				$candidatePos{$searchPos - 1} = $candidateIndel; # searchPos - 1 so that we're returning the anchor base
			}
			else
			{
				$done = 1;
			}
		}
	}
	
	
	return \%candidatePos;
	
}

# findSmallestInsertionSubunit returns the simplist subunit inside an insertion (e.g. for an insertion of +AAA, A is the simplist subunit)
# input is the insertion string
# output is the subunit string
sub findSmallestInsertionSubunit
{
	my $insertion = $_[0];

	my $sim;
	my $simFound = 0;

	my @subunits;

	for (my $i = 1; $i * 2 <= length($insertion); $i++)
	{
		if ((length($insertion) % $i) == 0)
		{
			$sim = substr($insertion, 0, $i);
			@subunits = ();

			for (my $j = length($sim); $j < length($insertion); $j += length($sim))
			{
				push(@subunits, substr($insertion, $j, length($sim)));
			}

			$simFound = 1;	  # innocent until guilty
			for (my $j = 0; $j < scalar(@subunits); $j++)
			{
				unless ($sim eq $subunits[$j])
				{
					$simFound = 0;
				}
			}

			if ($simFound == 1)
			{
				return $sim;
			}
		}
	}

	return $insertion;
}






sub procCigar
{
    my $cigar = $_[0];
    my @cigarOp;

    while ($cigar =~ /^([0-9]+[MIDNSHPX=]).*$/)
    {
        push (@cigarOp, $1);
        $cigar =~ s/$1//;
    }

    return @cigarOp;
}

sub procFlags
{
    if ($_[0] & 16)
    {
        return 0;
    }
    else
    {
        return 1;
    }
}


sub toPhred
{
    my $char = $_[0];
    my $ascii = ord($char);
    my $offset = 33;
    return $ascii - $offset;
}


sub getGenotype
{
	my $depth = $_[0];
	my $refCount = $_[1];
	my $altCount = $_[2];
	my $cuts = $_[3];

	my $genotype = "./.";

	if ($depth >= $cuts->{homoDepth})
	{
		if (($altCount / $depth) >= $cuts->{homoFreq})
		{
			$genotype = "1/1";
		}
		elsif (($refCount / $depth) >= $cuts->{homoFreq})
		{
			$genotype = "0/0";
		}
	}
	if ($depth >= $cuts->{hetDepth})
	{
		if ((($altCount / $depth) >= $cuts->{hetFreq}) and (($refCount / $depth) >= $cuts->{hetFreq}) and ($altCount >= $cuts->{hetBases}) and ($refCount >= $cuts->{hetBases}))
		{
			$genotype = "0/1";
		}
	}

	return $genotype;
}



# getBase returns the base at a specific position in the reference.  It will initialize fasta handles and pull new chunks of reference if necessary
# input is chromosome, position, reference hash and fasta handles
# output is a single base (and the reference hash and fasta handles may be modified)
sub getBase
{
    my $chr = $_[0];
    my $pos = $_[1];
    my $reference = $_[2];
    my $fastaHandles = $_[3];

    my $chunkStart = int(($pos - 1) / $reference->{"chunkSize"}) * $reference->{"chunkSize"} + 1;       # +1 because the first base in the reference is 1, $pos - 1 so that multiples of chunk size resolve to the correct chunk
    my $chunkEnd = $chunkStart + $reference->{"chunkSize"} - 1;

    unless (exists $reference->{$chr}{$chunkStart}{$pos})       # if the position isn't in our hash, we need to get a new chunk from the reference
    {
        unless (exists ($fastaHandles->{$chr}))     # create a handle for the chromosome fasta, if it doesn't exist
        {
#           warn "Creating fasta handle for $chr\n";
            $fastaHandles->{$chr} = Bio::DB::Fasta->new("$fastaHandles{path}/$chr.fa");
        }

#       warn "Pulling $chr:$chunkStart-$chunkEnd from fasta\n";
        my $newChunk = uc($fastaHandles{$chr}->seq($chr, $chunkStart, $chunkEnd));
        my $i = $chunkStart;
        for my $base (split("", $newChunk))
        {
            $reference->{$chr}{$chunkStart}{$i} = $base;
            $i++;
        }
    }
#   warn "returning $reference->{$chr}{$chunkStart}{$pos}\n";
    return $reference->{$chr}{$chunkStart}{$pos};
}


# getRange returns a string of bases from the reference in the specified range by calling getBase
# input is chromosome, start pos, end pos, reference hash and fasta handles
# output is a string of bases
sub getRange
{
    my $chr = $_[0];
    my $start = $_[1];
    my $end = $_[2];
    my $reference = $_[3];
    my $fastaHandles = $_[4];

    my $seq = "";

    for (my $p = $start; $p <= $end; $p++)
    {
        $seq .= getBase($chr, $p, $reference, $fastaHandles);
#       warn "Got base: $chr:$p\t$seq\n";
    }

    return $seq;
}


# freeMemChunks tests if the next indel to be processed is on a different chromosome or more than a chunk away from the reference sequences currently in memory
#   if there exist chunks that we won't need again (assuming the input is sorted) the chunks will be removed from the reference hash
# input is the chromosome and position of the current indel, and the reference hash
# there is no output
sub freeMemChunks
{
    my $chr = $_[0];
    my $pos = $_[1];
    my $reference = $_[2];

    # delete chunks from non-current chromosomes
    for my $refChr (keys %$reference)
    {
        if (($refChr ne $chr) and ($refChr ne "chunkSize"))
        {
#           warn "deleting all chunks for $refChr.\n";
            delete $reference->{$refChr};
        }
    }

    # delete chunks if they are more than 1.5 chunks away from the current indel
    # 1.5 so that we are at least in the middle of the current chunk before dropping the previous one
    for my $chunkPos (keys %{ $reference->{$chr} })
    {
        if ($chunkPos < ($pos - (1.5 * $reference->{"chunkSize"})))
        {
#           warn "deleting $chr:$chunkPos chunk.\n";
            delete $reference->{$chr}{$chunkPos};
        }
    }

    return;
}
