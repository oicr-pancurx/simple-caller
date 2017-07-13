#!/usr/bin/perl

use strict;
use warnings;

my @fileList = @ARGV;

my $somatotypeWarningLevel = 0.50;
my $minAssayedCount = 10;

my %laneHash;
my %calls;


my $lane;
my $sample;

my $type1;
my $type2;

my $l;
my ($chr, $pos, $id, $call);

for my $file (@fileList)
{
	if ($file =~ /\.\.\/\.\.\/(.*?)\/.*\/(.*?_.*?_.*?_.*?_.*?_.*?_.*?)_(.*?)_(.)_(.*?_.*?_.*?_.*?).somatotype/)
	{
		$lane = "$2_$5_$4_$3";
		$sample = $1;

		$laneHash{$sample}{$lane} = $file;
		open (FILE, "<$file") or die "Couldn't open $file\n";
		while ($l = <FILE>)
		{
			chomp $l;
			($chr, $pos, $id, $call) = split(/\t/, $l);
			$calls{$sample}{"$chr,$pos,$id"}{$file} = $call;
		}
		close FILE;
	}
	else
	{
		warn "couldn't parse $file , skipping\n";
	}
}

my $assayedCount;
my $matchedCount;

my $doWarning;

for $sample (keys %laneHash)
{
	$doWarning = 0;
	open (OUT, ">$sample.somatotype.csv") or die "Couldn't open >$sample.somatotype.csv\n";

	for $lane (sort keys %{ $laneHash{$sample} })
	{
		print OUT ",$lane";
	}
	print OUT "\n";

	for $lane (sort keys %{ $laneHash{$sample} })
	{
		$lane =~ /^.*?_.*?_.*?_(.*?)_/;
		$type1 = $1;

		print OUT "$lane";
		for my $lane2 (sort keys %{ $laneHash{$sample} })
		{
			$lane2 =~ /^.*?_.*?_.*?_(.*?)_/;
			$type2 = $1;

			$assayedCount = 0;
			$matchedCount = 0;

			for $pos (keys %{ $calls{$sample} })
			{
				if ((exists $calls{$sample}{$pos}{$laneHash{$sample}{$lane}}) and (exists $calls{$sample}{$pos}{$laneHash{$sample}{$lane2}}))
				{
					if (($calls{$sample}{$pos}{$laneHash{$sample}{$lane}} ne "") and ($calls{$sample}{$pos}{$laneHash{$sample}{$lane2}} ne ""))
					{
						$assayedCount++;
						if ($calls{$sample}{$pos}{$laneHash{$sample}{$lane}} eq $calls{$sample}{$pos}{$laneHash{$sample}{$lane2}})
						{
							$matchedCount++;
						}
					}
				}
			}
			unless ($assayedCount == 0)
			{
				printf OUT (",%.2f%% (%d/%d)", ($matchedCount / $assayedCount) * 100, $matchedCount, $assayedCount);
				if ($type1 eq $type2)
				{
					if ((($matchedCount / $assayedCount) < $somatotypeWarningLevel) and ($assayedCount > $minAssayedCount))
					{
						$doWarning = 1;
					}
				}
			}
			else
			{
				print OUT ",NA (0/0)";
			}
		}
		print OUT "\n";
	}

	if ($doWarning == 1)
	{
		print "$sample has a normal lane more than $somatotypeWarningLevel somatic, or a non-normal lane with less than $somatotypeWarningLevel somatic!\n";
		warn "$sample has a normal lane more than $somatotypeWarningLevel somatic, or a non-normal lane with less than $somatotypeWarningLevel somatic!\n";
	}

	print OUT "\n\n";

	print OUT "Chromosome,Position,Id";
	for $lane (sort keys %{ $laneHash{$sample} })
	{
		print OUT ",$lane";
	}
	print OUT "\n";

	for $pos (sort keys %{ $calls{$sample} })
	{
		print OUT "$pos";
		for $lane (sort keys %{ $laneHash{$sample} })
		{
			unless (exists $calls{$sample}{$pos}{$laneHash{$sample}{$lane}})
			{
				$calls{$sample}{$pos}{$laneHash{$sample}{$lane}} = "";
			}
			print OUT ",$calls{$sample}{$pos}{$laneHash{$sample}{$lane}}";
		}
		print OUT"\n";
	}

	close OUT;
}


open (OUT, ">all.somatotype.csv") or die "Couldn't open >all.somatotype.csv\n";

for $sample (sort keys %laneHash)
{
	for $lane (sort keys %{ $laneHash{$sample} })
	{
		print OUT ",$lane";
	}
}
print OUT "\n";

for $sample (sort keys %laneHash)
{
	for $lane (sort keys %{ $laneHash{$sample} })
	{
		print OUT "$lane";
		for my $sample2 (sort keys %laneHash)
		{
			for my $lane2 (sort keys %{ $laneHash{$sample2} })
			{
				$assayedCount = 0;
				$matchedCount = 0;
				for $pos (keys %{ $calls{$sample} })
				{
					if ((exists $calls{$sample}{$pos}{$laneHash{$sample}{$lane}}) and (exists $calls{$sample2}{$pos}{$laneHash{$sample2}{$lane2}}))
					{
						if (($calls{$sample}{$pos}{$laneHash{$sample}{$lane}} ne "") and ($calls{$sample2}{$pos}{$laneHash{$sample2}{$lane2}} ne ""))
						{
							$assayedCount++;
							if ($calls{$sample}{$pos}{$laneHash{$sample}{$lane}} eq $calls{$sample2}{$pos}{$laneHash{$sample2}{$lane2}})
							{
								$matchedCount++;
							}
						}
					}
				}
				unless ($assayedCount == 0)
				{
					printf OUT (",%.2f%% (%d/%d)", ($matchedCount / $assayedCount) * 100, $matchedCount, $assayedCount);
				}
				else
				{
					print OUT ",NA (0/0)";
				}
			}
		}
		print OUT "\n";
	}
}



