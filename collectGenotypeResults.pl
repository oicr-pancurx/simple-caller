#!/usr/bin/perl

use strict;
use warnings;

my @fileList = @ARGV;

my $genotypeWarningLevel = 0.85;
my $minAssayedCount = 10;

my %laneHash;
my %calls;


my $lane;
my $sample;

my $l;
my ($chr, $pos, $id, $call);

for my $file (@fileList)
{
	# SWID_497274_CPCG_0234_Pr_P_PE_332_WG_130228_SN801_0099_AC1FMVACXX_NoIndex_L008_R1_001.fastq.gz.annotated.bam
	if ($file =~ /SWID_.*?_(.*?_.*?_.*?_.*?_.*?_.*?_.*?)_(.*?_.*?_.*?_.*?)_(.*?)_L00(.)_.*/)
	{
		$lane = "$1_$2_$4_$3";
		$sample = "$1";
		$sample =~ s/^(.*?)_(.*?)_.*$/$1$2/;

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
	#PCSI_0019_Pa_P_PE_579_WG_NoIndex_5_130904_SN1080_0147_AC2D64ACXX.bam
	elsif ($file =~ /(.*?_.*?_.*?_.*?_.*?_.*?_.*?)_(.*?)_(.)_(.*?_.*?_.*?_.*?).bam/)
	{
		$lane = "$1_$4_$3_$2";
		$sample = $1;
		$sample =~ s/^(.*?)_(.*?)_.*$/$1$2/;

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
	elsif ($file =~ /\.\.\/\.\.\/(.*?)\/.*\/(.*?_.*?_.*?_.*?_.*?_.*?_.*?)_(.*?)_(.)_(.*?_.*?_.*?_.*?).genotype/)
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
	elsif ($file =~ /.*\/(.*?_.*?_.*?_.*?_.*?_.*?_.*?)_(.*?)_(.)_(.*?_.*?_.*?_.*?).somatotype/)
	{
		$lane = "$1_$4_$3_$2";
		$sample = $1;
		$sample =~ s/^(.*?)_(.*?)_.*$/$1$2/;

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
	open (OUT, ">$sample.genotype.csv") or die "Couldn't open >$sample.genotype.csv\n";

	for $lane (sort keys %{ $laneHash{$sample} })
	{
		print OUT ",$lane";
	}
	print OUT "\n";

	for $lane (sort keys %{ $laneHash{$sample} })
	{
		print OUT "$lane";
		for my $lane2 (sort keys %{ $laneHash{$sample} })
		{
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
				if ((($matchedCount / $assayedCount) < $genotypeWarningLevel) and ($assayedCount > $minAssayedCount))
				{
					$doWarning = 1;
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
		print "$sample has a lane with concordance less than $genotypeWarningLevel!\n";
		warn "$sample has a lane with concordance less than $genotypeWarningLevel!\n";
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


open (OUT, ">all.genotype.csv") or die "Couldn't open >all.genotype.csv\n";
open (OUTVIS, ">all-vis.genotype.csv") or die "Couldn't open >all.genotype.csv\n";

for $sample (sort keys %laneHash)
{
	for $lane (sort keys %{ $laneHash{$sample} })
	{
		print OUT ",${sample}_$lane";
		print OUTVIS ",${sample}_$lane";
	}
}
print OUT "\n";
print OUTVIS "\n";

for $sample (sort keys %laneHash)
{
	for $lane (sort keys %{ $laneHash{$sample} })
	{
		print OUT "${sample}_$lane";
		print OUTVIS "${sample}_$lane";
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
					printf OUTVIS (",%.2f", ($matchedCount / $assayedCount) * 100, $matchedCount, $assayedCount);
				}
				else
				{
					print OUT ",NA (0/0)";
					print OUTVIS ",NA";
				}

				if ($sample ne $sample2)
				{
					if (($assayedCount > 100) and (($matchedCount / $assayedCount) > 0.7))
					{
						warn "${sample}_$lane matches ${sample2}_$lane2 in $matchedCount / $assayedCount SNPs\n";
					}
				}
			}
		}
		print OUT "\n";
		print OUTVIS "\n";
	}
}


close OUT;
close OUTVIS;



