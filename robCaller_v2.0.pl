#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Std;
use vars qw/ %opt /;

use Bio::DB::Fasta;


my %cuts = (
	qual => 30,
	depth => 100,
	freq => 0.02
);

my %fastaHandles = ("path" => "/oicr/data/reference/genomes/homo_sapiens_mc/UCSC/hg19_random/Genomic/references/fasta/");
my $dbSNPfile = "/oicr/data/reference/genomes/homo_sapiens_mc/dbSNP/hg19_random/Genomic/dbSNP132/dbSNP132.vcf";
my $useUnion = 0;
my $unionPosFile;

my $readLength = 200;

my $outputPos = 0;

my $ignoreQual = 0;
sub usage
{
	print "\nUsage is samtools view -F 4 file.bam | robCaller_v2.0.pl [options].\n";
	print "Options are as follows:\n";
	print "\t-q minimum quality for a base to be considered.  Default is $cuts{qual}.\n";
	print "\t-d minimum depth of coverage for reporting a variant.  Default is $cuts{depth}.\n";
	print "\t-f minumum minor allele frequency for reporting.  Default is $cuts{freq}.\n";
	#print "\t-a all?\n";
	print "\t-F path to fasta files (one per chromosome, for now).  Default is $fastaHandles{\"path\"}.\n";
	print "\t-D dbSNP file for annotation.  Default is $dbSNPfile.\n";
	print "\t-P toggles output in .pos format.\n";
	print "\t-u unionPos file (for row by column).\n";
	print "\t-h displays this helpful usage.\n";
	
	print "\t-Q ignore quality scores\n"; 

	die "\n@_\n\n";
}

my $opt_string = "q:d:f:F:D:Pu:Qh";
getopts ($opt_string, \%opt) or usage("Incorrect arguments.");

if (exists $opt{h})
{
	usage("Help requested.");
}
if (exists $opt{Q})
{
	$ignoreQual = 1;
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
if (exists $opt{D})
{
	$dbSNPfile = $opt{D};
}
if (exists $opt{F})
{
	$fastaHandles{"path"} = $opt{F};
}
if (exists $opt{P})
{
	$outputPos = 1;
}
if (exists $opt{u})
{
	$useUnion = 1;
	$unionPosFile = $opt{u};

	$cuts{depth} = 0;
	$cuts{freq} = 0;
}

my $line;
my @fields;

my %reference;			# $reference{"chr\tpos"} = base;

my %dbSNPhash;

#open (DBSNPFILE, "$dbSNPfile") or die "Couldn't open $dbSNPfile.\n";
#while ($line = <DBSNPFILE>)
#{
#	unless ($line =~ /^#/)
#	{
#		chomp $line;
#		@fields = split(/\t/, $line);
#
#		$dbSNPhash{"chr$fields[0]\t$fields[1]\t$fields[3]\t$fields[4]"} = $fields[2];
#	}
#}
#close DBSNPFILE;

my %unionPosHash;
$unionPosHash{"useMe"} = $useUnion;
if ($useUnion)
{
	open (UNIONFILE, $unionPosFile) or die "Couldn't open $unionPosFile .\n";
	while ($line = <UNIONFILE>)
	{
		chomp $line;
		$unionPosHash{"$line"} = 1;
	}
}


my %posHash;
# $posHash{position}{base}{qual}{start point} = count

my %indelHash;
# $indelHash{position}{indel[+A or -TCG or "no indel"]}{start point}{$isForward} = count

my $indel;
my $delCount;

my $chr;
my $startPos;
my $pPos;
my $pChr;
my $flags;
my $isForward;
my $cigar;
my @cigarOp;
my @bases;
my @quals;

my $base;
my $q;
my $refPos;
my $indelPos;

my $dumpToFileCounter = 0;
my $dumpTime = 5000;

#my @qualCuts = qw(25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40);

while (defined($line = <STDIN>))
{
	chomp $line;
	unless ($line =~ /^@/)
	{
		@fields = split(/\t/, $line);
	
		$flags = $fields[1];
		$chr = $fields[2];
		$startPos = $fields[3];
		$cigar = $fields[5];
		@bases = split(//, $fields[9]);
		unless ($ignoreQual == 1)
		{
			@quals = split(//, $fields[10]);
		}
		else
		{
			$q = 0;
		}
	
		$isForward = procFlags($flags);
	
		@cigarOp = procCigar($cigar);
	
		$refPos = $startPos;
	
		if ($cigarOp[0] =~ /(.*)S/)		# if first cigar operation is a soft clip, adjust start point
		{
			$startPos -= $1;
		}
	
		foreach my $c (@cigarOp)
		{
			if ($c =~ /(.*)M/)
			{
				for (my $i = 0; $i < $1; $i++)
				{
					$base = shift(@bases);
		unless ($ignoreQual == 1)
		{
					$q = shift(@quals);
		}	
					unless ($isForward == 1)
					{
						$base = lc($base);
					}
		
					unless (exists $posHash{"$chr\t$refPos"}{$base}{$q}{$startPos})
					{
						$posHash{"$chr\t$refPos"}{$base}{$q}{$startPos} = 1;
					}
					else
					{
						$posHash{"$chr\t$refPos"}{$base}{$q}{$startPos} += 1;
					}
	
					if (exists $indelHash{"$chr\t$refPos"}{"no indel"}{$startPos}{$isForward})
					{
						$indelHash{"$chr\t$refPos"}{"no indel"}{$startPos}{$isForward} += 1;
					}
					else
					{
						$indelHash{"$chr\t$refPos"}{"no indel"}{$startPos}{$isForward} = 1;
					}

					$refPos++;
				}
			}
			elsif ($c =~ /(.*)S/)
			{
				for (my $i = 0; $i < $1; $i++)
				{
					$base = shift(@bases);
		unless ($ignoreQual == 1)
		{
					$q = shift(@quals);
		}
				}
			}
			elsif ($c =~ /(.*)I/)
			{
				$indel = "+";
				$indelPos = $refPos - 1;
				for (my $i = 0; $i < $1; $i++)
				{
					$base = shift(@bases);
		unless ($ignoreQual == 1)
		{
					$q = shift(@quals);
		}
			$indel .= $base;
				}
				if (exists $indelHash{"$chr\t$indelPos"}{$indel}{$startPos}{$isForward})
				{
					$indelHash{"$chr\t$indelPos"}{$indel}{$startPos}{$isForward} += 1;
				}
				else
				{
					$indelHash{"$chr\t$indelPos"}{$indel}{$startPos}{$isForward} = 1;
				}
			}
			elsif ($c =~ /(.*)D/)
			{
				$delCount = 0;
				$indelPos = $refPos;
				for (my $i = 0; $i < $1; $i++)
				{
					if (exists $indelHash{"$chr\t$refPos"}{"no indel"}{$startPos}{$isForward})		# adding to "no indel" because this is used later as a raw read count - can't exactly use deletions because they are generalized and the counts would move to different positions
					{
						$indelHash{"$chr\t$refPos"}{"no indel"}{$startPos}{$isForward} += 1;
					}
					else
					{
						$indelHash{"$chr\t$refPos"}{"no indel"}{$startPos}{$isForward} = 1;
					}
					$refPos++;
					$delCount++
				}
				$indel = "-$delCount";
				if (exists $indelHash{"$chr\t$indelPos"}{$indel}{$startPos}{$isForward})
				{
					$indelHash{"$chr\t$indelPos"}{$indel}{$startPos}{$isForward} += 1;
				}
				else
				{
					$indelHash{"$chr\t$indelPos"}{$indel}{$startPos}{$isForward} = 1;
				}
			}
			elsif ($c =~ /(.*)H/)
			{
				# don't need to do anything?
			}
			else
			{
				die "Can't handle CIGAR operation: $c\n";
			}
		}
		
		$dumpToFileCounter++;
	
		if ($dumpToFileCounter >= $dumpTime)
		{
			# print file and clean up hash
	
			foreach my $p (sort keys %posHash)
			{
				@fields = split(/\t/, $p);
				$pChr = $fields[0];
				$pPos = $fields[1];
				
	
				if (($pPos < $startPos ) or ($pChr ne $chr)) 		# +1000 to leave space for indels
				{
					printPosition(\%posHash, $p, \%reference, \%fastaHandles, $outputPos, \%cuts, \%unionPosHash);
					delete $posHash{$p};
				}
			}
	
			$dumpToFileCounter = 0;
		}
	}
}


foreach my $p (sort keys %posHash)
{
	printPosition(\%posHash, $p, \%reference, \%fastaHandles, $outputPos, \%cuts, \%unionPosHash);
}

procAndPrintIndels(\%indelHash, \%reference, \%fastaHandles, $outputPos, \%cuts, \%unionPosHash);



sub printPosition
{
	my $posHashRef = $_[0];
	my $p = $_[1];
	my $refHashRef = $_[2];
	my $fastaHandlesRef = $_[3];
	my $outputPos = $_[4];
	my $cuts = $_[5];
	my $positionList = $_[6];

	my $refBase = uc(getBase($p, $refHashRef, $fastaHandlesRef));


if ($outputPos == 1)
{
	my @qualCuts = qw(25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40);




	# get total depth, read counts and start point histograms at current pos
	my $totalDepth = 0;
	my %readCounts;
	my %startPointHists;

	foreach my $qualCut (@qualCuts)
	{
		foreach my $base (qw(A a C c G g T t))
		{
			$readCounts{$qualCut}{$base} = 0;
		}
	}

	foreach my $base (qw(A a C c G g T t))
	{
		foreach my $qual (keys %{ $posHashRef->{$p}{$base} })
		{
			foreach my $stpnt (keys %{ $posHashRef->{$p}{$base}{$qual} })
			{
				$totalDepth += $posHashRef->{$p}{$base}{$qual}{$stpnt};

				foreach my $qualCut (@qualCuts)
				{
					if (toPhred($qual) >= $qualCut)
					{
						$readCounts{$qualCut}{$base} += $posHashRef->{$p}{$base}{$qual}{$stpnt};


						unless (exists $startPointHists{$qualCut}{$base}{$stpnt})
						{
							$startPointHists{$qualCut}{$base}{$stpnt} = $posHashRef->{$p}{$base}{$qual}{$stpnt};
						}
						else
						{
							$startPointHists{$qualCut}{$base}{$stpnt} += $posHashRef->{$p}{$base}{$qual}{$stpnt};
						}
					}

				}

				
			}
		}
	}


	# calculate quality depth, MAF, T/F/r, strand bias metric
	
	my %qualDepths;
	my %MAFs;
	my %TFrs;
	my %strandBias;

	my $forwardFreq;
	my $reverseFreq;

	foreach my $qualCut (keys %readCounts)
	{
		foreach my $base (qw(A a C c G g T t))
		{
			unless (exists $qualDepths{$qualCut})
			{
				$qualDepths{$qualCut} = $readCounts{$qualCut}{$base};
			}
			else
			{
				$qualDepths{$qualCut} += $readCounts{$qualCut}{$base};
			}
		}

		foreach my $base (qw(A C G T))
		{
			if ($qualDepths{$qualCut} > 0)
			{
				$MAFs{$qualCut}{$base} = ($readCounts{$qualCut}{uc($base)} + $readCounts{$qualCut}{lc($base)}) / $qualDepths{$qualCut};
			}
			else
			{
				$MAFs{$qualCut}{$base} = 0;
			}
				$TFrs{$qualCut}{$base} = ($readCounts{$qualCut}{uc($base)} + $readCounts{$qualCut}{lc($base)}) . "/" . $readCounts{$qualCut}{uc($base)} . "/" . $readCounts{$qualCut}{lc($base)};
			
			
			if (($readCounts{$qualCut}{uc($base)} + $readCounts{$qualCut}{lc($base)}) > 0)
			{
				$forwardFreq = ($readCounts{$qualCut}{uc($base)} / ($readCounts{$qualCut}{uc($base)} + $readCounts{$qualCut}{lc($base)}));
				$reverseFreq = ($readCounts{$qualCut}{lc($base)} / ($readCounts{$qualCut}{uc($base)} + $readCounts{$qualCut}{lc($base)}));
				$strandBias{$qualCut}{$base} = $forwardFreq * $reverseFreq;
			}
			else
			{
				$strandBias{$qualCut}{$base} = 0;
			}
		}
	}

	# calculate # start points and upper bin metric
	
	my %startPointsCollect;
	my %upperBinsCollect;	# wherein we do the collecting
	my %startPoints;
	my %upperBins;

	my $numStartPoints;

	foreach my $qualCut (@qualCuts)		#keys %startPointHists)		# $startPointHists{$qualCut}{$base}{$stpnt}
	{
		foreach my $base (qw(A a C c G g T t))
		{
			$startPointsCollect{$qualCut}{$base} = 0;
		}
		foreach my $base (qw(A C G T))
		{
			delete $upperBinsCollect{$qualCut}{$base};
		}

		if (exists $startPointHists{$qualCut})
		{
			foreach my $base (keys %{ $startPointHists{$qualCut} })
			{

				foreach my $stpnt (keys %{ $startPointHists{$qualCut}{$base} })
				{
					$startPointsCollect{$qualCut}{$base}++;
					unless (exists $upperBinsCollect{$qualCut}{uc($base)}{$stpnt})
					{
						$upperBinsCollect{$qualCut}{uc($base)}{$stpnt} = $startPointHists{$qualCut}{$base}{$stpnt};
					}
					else
					{
						$upperBinsCollect{$qualCut}{uc($base)}{$stpnt} += $startPointHists{$qualCut}{$base}{$stpnt};
					}
				}
			}
		}

		foreach my $base (qw(A C G T))
		{
			$startPoints{$qualCut}{$base} = ($startPointsCollect{$qualCut}{uc($base)} + $startPointsCollect{$qualCut}{lc($base)}) . "/" . $startPointsCollect{$qualCut}{uc($base)} . "/" . $startPointsCollect{$qualCut}{lc($base)};

			if (exists $upperBinsCollect{$qualCut}{$base} )
			{
				$upperBins{$qualCut}{$base} = upperBin(\%{ $upperBinsCollect{$qualCut}{$base} }, 10);
			}
			else
			{
				$upperBins{$qualCut}{$base} = 0;
			}
		}


	}


	# output results
	my $output;

	my $worthReporting;
	my $worthReportingMAF = $cuts->{freq};
	foreach my $base (qw(A C G T))
	{
		unless ($base eq uc($refHashRef->{$p}))
		{
			$worthReporting = 0;
			$output = "$p\t" . $refHashRef->{$p} . "\t$base\t";
	
			if (exists $dbSNPhash{"$p\t$refBase\t$base"})
			{
				$output = $output . $dbSNPhash{"$p\t$refBase\t$base"};
			}
			else
			{
				$output = $output . ".";
			}
	
			if ($positionList->{"useMe"} == 1)
			{
				if (exists($positionList->{$output}))
				{
					$worthReporting = 1;
				}
				else
				{
					$worthReporting = 0;
				}

			}

			$output = $output . "\t$totalDepth";
			foreach my $qualCut (@qualCuts) #qw(30 35))
			{
				unless ($qualDepths{$qualCut} < 1)
				{
					$output = $output . "\tq$qualCut";
					$output = $output . ",$qualDepths{$qualCut}";
					$output = $output . ",$MAFs{$qualCut}{$base}";
					$output = $output . ",$strandBias{$qualCut}{$base}";
					$output = $output . ",$upperBins{$qualCut}{$base}";
					$output = $output . ",$strandBias{$qualCut}{$refBase}";
					$output = $output . ",$upperBins{$qualCut}{$refBase}";
	
					$output = $output . ",$TFrs{$qualCut}{$base}";
					$output = $output . ",$startPoints{$qualCut}{$base}";
					$output = $output . ",$TFrs{$qualCut}{$refBase}";
					$output = $output . ",$startPoints{$qualCut}{$refBase}";
	
					if (($MAFs{$qualCut}{$base} > $worthReportingMAF) and ($positionList->{"useMe"} == 0))
					{
						$worthReporting = 1;
					}
				}
			}
			if ($worthReporting == 1)
			{
				print $output . "\n";
			}
		}
	}
}
else		# more simple output
{
        # get total depth and base counts
        my $totalDepth = 0;
        my %baseCounts;

        foreach my $base (qw(A C G T))
        {
                $baseCounts{$base} = 0;
        }

        foreach my $base (qw(A a C c G g T t))
        {
                foreach my $qual (keys %{ $posHashRef->{$p}{$base} })
                {
                        foreach my $stpnt (keys %{ $posHashRef->{$p}{$base}{$qual} })
                        {
                                if (toPhred($qual) >= $cuts->{qual})
                                {
                                	$totalDepth += $posHashRef->{$p}{$base}{$qual}{$stpnt};

									$baseCounts{uc($base)} += $posHashRef->{$p}{$base}{$qual}{$stpnt};
								}
                        }
                }
        }


		my $MAF;

		foreach my $base (qw(A C G T))
		{
			unless ($base eq $refBase)
			{
				if ($totalDepth >= $cuts->{depth})
				{
					unless ($totalDepth == 0)
					{
						$MAF = $baseCounts{$base} / $totalDepth;
					}
					else
					{
						$MAF = 0;
					}

					if ($MAF >= $cuts->{freq})
					{
						print "$p\t$refBase\t$base\t.\t$MAF\t$totalDepth";

						foreach my $moreBase (qw(A C G T))
						{
							unless ($totalDepth == 0)
							{
								$MAF = $baseCounts{$moreBase} / $totalDepth;
								print "\t$MAF";
							}
							else
							{
								$MAF = 0;
								print "\t$MAF";
							}
						}

						print "\n";
					}
				}
			}
		}


}



}

sub procAndPrintIndels
{
	my $indelHash = $_[0];
	my $reference = $_[1];
	my $fastaHandles = $_[2];

	my $outputPos = $_[3];
	my $cuts = $_[4];

	my $positionList = $_[5];

	my @qualCuts = qw(25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40);

	my $delCount;
	my $deletion;

	my $chr;
	my $pos;
	my $before;
	my $after;

	# use reference to populate deletion bases
	for my $p (keys %{$indelHash})
	{
		for my $indel (keys %{ $indelHash->{$p} })
		{
			if ($indel =~ /^-(.*)/)		# is a deletion
			{
				$delCount = $1;
				for my $sp (keys %{ $indelHash->{$p}{$indel} })
				{
					for my $isF (keys %{ $indelHash->{$p}{$indel}{$sp} })
					{
						$deletion = "-";
						for (my $i = 0; $i < $delCount; $i++)
						{
							($chr, $pos) = split(/\t/, $p);
							$pos += $i;

							$deletion .= uc(getBase("$chr\t$pos", $reference, $fastaHandles));
						}

						$indelHash->{$p}{$deletion}{$sp}{$isF} = $indelHash->{$p}{$indel}{$sp}{$isF};
						delete $indelHash->{$p}{$indel}{$sp}{$isF};
					}
				}
			}
		}
	}

	my $generalIndels = {};
	my $indel;

	# generalize indels
	for my $p (keys %{$indelHash})
	{
		for my $i (keys %{ $indelHash->{$p} })
		{
			for my $sp (keys %{ $indelHash->{$p}{$i} })
			{
				for my $isF (keys %{ $indelHash->{$p}{$i}{$sp} })
				{
					unless ($i eq "no indel")
					{
						($chr, $pos, $indel, $before, $after) = generalizeIndel($i, $p, $reference, $fastaHandles);		# would be better to not do this for each start point

						unless ($chr eq "skip")
						{
							$generalIndels->{"$chr\t$pos"}{"$before\t$after ($i)"}{$sp}{$isF} += $indelHash{$p}{$i}{$sp}{$isF};
						}
						delete $indelHash{$p}{$i}{$sp}{$isF};
					}
					else
					{
						$generalIndels->{$p}{"$i\t"}{$sp}{$isF} += $indelHash{$p}{$i}{$sp}{$isF};
						delete $indelHash{$p}{$i}{$sp}{$isF};
					}
				}
			}
		}
	}


	my $worthReporting;
	my $totalDepth;
	my $indelDepth;

	my $indelFrequency;
	my $notIndelFrequency;

	my $output;

	# process and print each indel
	for my $p (sort keys %{$generalIndels})
	{
		for my $i (keys %{ $generalIndels->{$p} })
		{
			$worthReporting = 1;
			$totalDepth = 0;
			$indelDepth = 0;

			# calculate total depth			# ideally would only do this once per position
			for my $j (keys %{$generalIndels->{$p} })
			{
				if ($j eq "no indel\t") 		# only count reads!
				{
					for my $sp (keys %{ $generalIndels->{$p}{$j} })
					{
						for my $isF (keys %{ $generalIndels->{$p}{$j}{$sp} })
						{
							$totalDepth += $generalIndels->{$p}{$j}{$sp}{$isF};
						}
					}
				}
			}

			# calculate indel depth
			for my $sp (keys %{ $generalIndels->{$p}{$i} })
			{
				for my $isF (keys %{ $generalIndels->{$p}{$i}{$sp} })
				{
					$indelDepth += $generalIndels->{$p}{$i}{$sp}{$isF};
				}
			}

			if ($totalDepth > 0)
			{
				$indelFrequency = $indelDepth / $totalDepth;
			}
			else
			{
				$indelFrequency = 1;
			}
			$notIndelFrequency = 1 - $indelFrequency;



			if ($outputPos == 1)
			{
				if (($totalDepth > $cuts->{depth}) and ($indelFrequency >= $cuts->{freq}) and ($i ne "no indel\t"))
				{
					# calculate start point bias and whatnot
					
					
					$output = "$p\t" . "$i\t";
					$output .= ".\t";	# dbSNP
					$output .= "$totalDepth\t";
	
					$output .= "qI,";		# qual cut
					$output .= "$totalDepth,";
				$output .= "$indelFrequency,";

					$output .= "1,";		# strand bias minor
					$output .= "1,";		# sp bias minor
					$output .= "1,";		# strand bias ref
					$output .= "1,";		# sp bias ref

					$output .= "0/0/0,0/0/0,0/0/0,0/0/0";		# forward/reverse base and sp counts for ref and minor

					if ($positionList->{"useMe"} == 1)
					{
						if (exists $positionList->{"$p\t$i\t."})		# . for dbSNP
						{
							print $output . "\n";
						}
					}
					else
					{
						print $output . "\n";
					}
				}

			}
			else
			{
				if (($totalDepth > $cuts->{depth}) and ($indelFrequency >= $cuts->{freq}) and ($i ne "no indel\t"))
				{
					# print
					$output = "$p\t" . "$i\t";

					$output .= ".\t";	# dbSNP
					$output .= "$indelFrequency\t";
					$output .= "$totalDepth\n";


					print $output;
				}
			}
		}
	}

}


sub getBase
{
	my $chrPos = $_[0];
	my $reference = $_[1];
	my $fastaHandles = $_[2];

	if (exists $reference->{$chrPos})
	{
		return $reference->{$chrPos};
	}

	my ($chr, $pos) = split(/\t/, $chrPos);

	unless (exists ($fastaHandles->{$chr}))
	{
		$fastaHandles->{$chr} = Bio::DB::Fasta->new("$fastaHandles{path}/$chr.fa");
	}
	unless (exists ($reference->{"$chrPos"}))
	{
		$reference{$chrPos} = $fastaHandles{$chr}->seq($chr, $pos, $pos);
	}
	return $reference->{$chrPos};
}


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
                $seq .= uc(getBase("$chr\t$p", $reference, $fastaHandles));
        }

        return $seq;
}

sub generalizeIndel
{
	my $type;
	my $indel;


	if ($_[0] =~ /([+-])([ACGT]+)/)      # in the format +A or -TAG
	{
	        $type = $1;
	        $indel = $2;
	}
	else
	{
	        warn "Couldn't generalize indel [$_[1]: $_[0]], should be in the style of +A or -TAG\n";
			return ("skip", "skip", "skip", "skip", "skip");
	}
	

	my $indelLength = length($indel);
	
	my $chrPos = $_[1];
	
	my ($chr, $pos) = split(/\t/, $chrPos);
	
	my $referenceHash = $_[2];
	my $fastaHandles = $_[3];
	
	my $before;
	my $after;
	
	my $sortedIndel = join('',sort(split('',$indel)));
	
	my %candidatePos;
	my $candidateIndel;
	
	my $smallestInsertionSubunit;
	
	if ($type eq "+")       # won't be able to find the insertion in reference seqeuence
	{
	        $candidatePos{$pos} = $type . $indel;
	        $smallestInsertionSubunit = findSmallestInsertionSubunit($indel);
	}
	else
	{
	        # check if deletion is possible
	        $candidateIndel = getRange($chr, $pos, $pos + $indelLength - 1, $referenceHash, $fastaHandles);
	
	        if ($indel eq $candidateIndel)
	        {
	                $candidatePos{$pos} = $type . $indel;
	        }
	        else
	        {
	                die "Deletion of $indel at $chr $pos ($candidateIndel) is not possible\n";
	        }
	}
	
	# search left
	my $done = 0;
	my $searchPos = $pos;
	
	my $lastFoundInsertion = $pos;
	
	while ($done == 0)
	{
	
	
	        if ($type eq "+")
	        {
	                $searchPos -= length($smallestInsertionSubunit);
	                $candidateIndel = getRange($chr, $searchPos + 1, $searchPos + 1 + length($smallestInsertionSubunit) - 1, $referenceHash, $fastaHandles);  # check the bases to the right of the potential insertion for a subunit match
	
	                if ($candidateIndel eq $smallestInsertionSubunit)
	                {
	                        $candidatePos{$searchPos} = $type . $indel;
	
	                }
	                else
	                {
	                                $done = 1;
	                }
	        }
	        else    # deletion
	        {
	                $searchPos--;
	                $candidateIndel = getRange($chr, $searchPos, $searchPos + $indelLength - 1, $referenceHash, $fastaHandles);  # old code: substr($reference, $searchPos, $indelLength);

	                if (join('',sort(split('',$candidateIndel))) eq $sortedIndel)
	                {
	                        $candidatePos{$searchPos} = $type . $candidateIndel;
	                }
	                else
	                {
	                        $done = 1;
	                }
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
	                $candidateIndel = getRange($chr, $searchPos + 1, $searchPos + 1 + length($smallestInsertionSubunit) - 1, $referenceHash, $fastaHandles);  # search the bases to the right of the position base for a subunit match
	
	                if ($candidateIndel eq $smallestInsertionSubunit)
	                {
	                        $candidatePos{$searchPos} = $type . $indel;
	                        $candidatePos{$searchPos + length($smallestInsertionSubunit)} = $type . $indel;         # can insert on either side of a repeated motif
	
	                }
	                else
	                {
	                                $done = 1;
	                }
	        }
	        else    # deletion
	        {
	                $searchPos++;
	                $candidateIndel = getRange($chr, $searchPos, $searchPos + $indelLength - 1, $referenceHash, $fastaHandles);  # old code: substr($reference, $searchPos, $indelLength);

	                if (join('',sort(split('',$candidateIndel))) eq $sortedIndel)
	                {
	                        $candidatePos{$searchPos} = $type . $candidateIndel;
	                }
	                else
	                {
	                        $done = 1;
	                }
	        }
	}
	
	# return results
	my @sortedPositions = sort { $a <=> $b } keys %candidatePos;
	
	if ($type eq "+")
	{
	        $before = getRange($chr, $sortedPositions[0] + 1, $sortedPositions[$#sortedPositions], $referenceHash, $fastaHandles);        # old code: substr($reference, $sortedPositions[0], $sortedPositions[$#sortedPositions] - $sortedPositions[0]);
	        $after = "$before$indel";
	        if ($before eq "")
	        {
	                $before = "-";
	        }
	}
	elsif ($type eq "-")
	{
	        $before = getRange($chr, $sortedPositions[0], $sortedPositions[$#sortedPositions] + $indelLength - 1, $referenceHash, $fastaHandles); # old code: substr($reference, $sortedPositions[0], $sortedPositions[$#sortedPositions] - ($sortedPositions[0] - $indelLength));
	        $after = getRange($chr, $sortedPositions[0], $sortedPositions[$#sortedPositions] - 1, $referenceHash, $fastaHandles); # old code: substr($reference, $sortedPositions[0], $sortedPositions[$#sortedPositions] - $sortedPositions[0]);
	        if ($after eq "")
	        {
	                $after = "-";
	        }
	}
	
	# return left format for now, I guess (better to be consistent with dbSNP?)
	$pos = $sortedPositions[0];
	
	return ($chr, $pos, $candidatePos{$pos}, $before, $after);

}


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

                        $simFound = 1;          # innocent until guilty
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

sub mean
{
	my %hash = %{ $_[0] };
	my $sum = 0;
	my $count = 0;

	foreach my $h (keys %hash)
	{
		$sum += $hash{$h};
		$count++;
	}

	return ($sum / $count);
}

sub median
{
	my %hash = %{ $_[0] };
	my $count = 0;

	my @array;

	foreach my $h (keys %hash)
	{
		push(@array, $hash{$h});
		$count++;
	}
	my @sortArray = sort { $a <=> $b } @array;

	return $sortArray[$count / 2];

}

sub count
{
	my %hash = %{ $_[0] };
	my $count = 0;
	foreach my $h (keys %hash)
	{
		$count++;
	}
	return $count;
}

sub min
{
	my %hash = %{ $_[0] };
	my $min = 999999999999999;
	foreach my $h (keys %hash)
	{
		if ($hash{$h} < $min)
		{
			$min = $hash{$h};
		}
	}
	return $min;
}

sub max
{
	my %hash = %{ $_[0] };
	my $max = 0;
	foreach my $h (keys %hash)
	{
		if ($hash{$h} > $max)
		{
			$max = $hash{$h};
		}
	}
	return $max;
}

sub upperBin
{
	my %hash = %{ $_[0] };
	my $numBins = $_[1];

	my $hashSize = keys %hash;
	if ($hashSize < 1)
	{
		return 0;
	}


	my @sortArray = sort {$hash{$b} <=> $hash{$a}} keys %hash;

	my $total = 0;
	foreach my $h (keys %hash)
	{
		$total += $hash{$h};
	}
	my $binSize = $total / $numBins;

	my $iInBin = 0;
	my $sum = 0;

	foreach my $i (@sortArray)
	{
		if (($hash{$i} + $sum) > $binSize)
		{
			$iInBin += ($binSize - $sum) / $hash{$i};
			return $iInBin;
		}
		else
		{
			$iInBin++;
			$sum += $hash{$i};
		}
	}
	
}

sub toPhred
{
	my $char = $_[0];
	my $ascii = ord($char);
	my $offset = 33;
	return $ascii - $offset;
}

