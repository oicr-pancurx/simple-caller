#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

my $vcfTarget = shift(@ARGV);
my @bams = @ARGV;

my $simpleCaller = "/.mounts/labs/PCSI/users/rdenroche/simpleCaller/simpleCaller.pl";
my $bedtools = "/oicr/local/analysis/sw/bedtools/2.24.0/bedtools2/bin/bedtools";
my $samtools = "/oicr/local/analysis/sw//samtools/samtools-0.1.18/samtools";
my $outDir = "sc_out";
mkdir $outDir;

my %vars;
my %posHash;

my $l;
my ($chr, $pos, $id, $ref, $alts, $alt, $qual, $filter, $info, $format, $gt, $gt2, $varKey);

# parse vcf
open (VCF, $vcfTarget) or die "Couldn't open $vcfTarget\n";

while ($l = <VCF>)
{
	chomp $l;
	unless ($l =~ /^#/)
	{
		($chr, $pos, $id, $ref, $alts, $qual, $filter, $info, $format, $gt, $gt2) = split(/\t/, $l);

		for $alt (split(/,/, $alts))
		{
			$varKey = "$chr\t$pos\t$ref\t$alt";

			$vars{$varKey}{chr} = $chr;
			$vars{$varKey}{pos} = $pos;
			$vars{$varKey}{ref} = $ref;
			$vars{$varKey}{alt} = $alt;
			$vars{$varKey}{info} = $info;

			$posHash{"$chr\t$pos"}++;
		}
	}

}
close VCF;


# build bed file from vcf
my $bedFile = "./$vcfTarget.bed";
$bedFile =~ s/^.*\//.\/$outDir\//;

my $mergedBed = "./$vcfTarget.merged.bed";
$mergedBed =~ s/^.*\//.\/$outDir\//;

my ($start, $end);

open (BED, ">$bedFile") or die "Couldn't open $bedFile for writing\n";

for my $v (sort keys %vars)
{
	$start = $vars{$v}{pos};
	$end = $vars{$v}{pos} + length($vars{$v}{ref});
	print BED "$vars{$v}{chr}\t$start\t$end\n";
}
close BED;

`$bedtools sort -i $bedFile | $bedtools merge > $mergedBed`;


# run simpleCaller to generate single sample vcfs
my $outVcf;
my @outVcfList;
for my $bam (@bams)
{
	$outVcf = "./$bam.sc.vcf";
	$outVcf =~ s/^.*\//.\/$outDir\//;

	push (@outVcfList, $outVcf);

	`$samtools view -hF4 -L $mergedBed $bam | $simpleCaller -A > $outVcf`;
}


# create multi-sample vcf
my $sample;
my @samples;
my %jointHash;
for my $bamVcf (@outVcfList)
{
	open (VCF, $bamVcf) or die "Couldn't open $bamVcf\n";
	$sample = ".";
	while ($l = <VCF>)
	{
		chomp $l;
		if ($l =~ /^#/)
		{
			if ($l =~ /##sample=<ID=(.*?)>/)
			{
				$sample = $1;
				push(@samples, $sample);
			}
		}
		else
		{
			($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $gt) = split(/\t/, $l);
			if (exists $posHash{"$chr\t$pos"})
			{
				if ($format eq "GT:A:C:G:T:N:del")
				{
					if ($info =~ /DP=(.*?);/)
					{
						$jointHash{"$chr\t$pos"}{$sample}{depth} = $1;
					}
					else
					{
						warn "Couldn't parse $info\n";
					}

					if ($gt =~ /^.*?:(.*?):(.*?):(.*?):(.*?):/)
					{
						$jointHash{"$chr\t$pos"}{$sample}{A} = $1;
						$jointHash{"$chr\t$pos"}{$sample}{C} = $2;
						$jointHash{"$chr\t$pos"}{$sample}{G} = $3;
						$jointHash{"$chr\t$pos"}{$sample}{T} = $4;
					}
					else
					{
						warn "Couldn't parse $gt\n";
					}
				}
				elsif ($format eq "GT:del:other")
				{
					if ($info =~ /DP=(.*?);/)
					{
						$jointHash{"$chr\t$pos"}{$sample}{depth} = $1;
					}
					if ($gt =~ /^.*?:(.*?):(.*?)$/)
					{
						$jointHash{"$chr\t$pos"}{$sample}{"$ref>$alt"} = $1;
						$jointHash{"$chr\t$pos"}{$sample}{"!$ref>$alt"} = $2;
					}
					
				}
				elsif ($format eq "GT:ins:other")
				{
					if ($info =~ /DP=(.*?);/)
					{
						$jointHash{"$chr\t$pos"}{$sample}{depth} = $1;
					}
					if ($gt =~ /^.*?:(.*?):(.*?)$/)
					{
						$jointHash{"$chr\t$pos"}{$sample}{"$ref>$alt"} = $1;
						$jointHash{"$chr\t$pos"}{$sample}{"!$ref>$alt"} = $2;
					}
				}
				else
				{
					warn "Didn't match $format\n";
				}
			}
		}
	}
	close VCF;
}

print "##fileformat=VCFv4.1\n";
print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
for $sample (sort @samples)
{
	print "\t$sample";
}
print "\n";

my $freq;

for my $v (sort keys %vars)
{
	$chr = $vars{$v}{chr};
	$pos = $vars{$v}{pos};
	$ref = $vars{$v}{ref};
	$alt = $vars{$v}{alt};
	$info = $vars{$v}{info};

	print "$chr\t$pos\t.\t$ref\t$alt\t.\t.\t$info\talt_freq:alt_depth:total_depth";

	for $sample (sort @samples)
	{
		if (exists $jointHash{"$chr\t$pos"}{$sample})
		{
			if (length($ref) == length($alt))
			{
				$freq = $jointHash{"$chr\t$pos"}{$sample}{$alt} / $jointHash{"$chr\t$pos"}{$sample}{depth};
				print "\t$freq:$jointHash{\"$chr\t$pos\"}{$sample}{$alt}:$jointHash{\"$chr\t$pos\"}{$sample}{depth}";
			}
			elsif (exists $jointHash{"$chr\t$pos"}{$sample}{"$ref>$alt"})
			{
				$freq = $jointHash{"$chr\t$pos"}{$sample}{"$ref>$alt"} / $jointHash{"$chr\t$pos"}{$sample}{depth};
				print "\t$freq:$jointHash{\"$chr\t$pos\"}{$sample}{\"$ref>$alt\"}:$jointHash{\"$chr\t$pos\"}{$sample}{depth}";
			}
			else
			{
				print "\t0:0:$jointHash{\"$chr\t$pos\"}{$sample}{depth}";
			}
		}
		else
		{
			print "\t.:.:.";
		}
	}
		print "\n";
}






