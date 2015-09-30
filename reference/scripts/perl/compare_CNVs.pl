#!/usr/bin/perl -w

use strict;
use Math::BigInt;
use Math::BigFloat;

use CNV_utils;

my $truthFile = shift || undef;

my $testFile = shift || undef;
my $testFamFile = shift || undef;

my $targetFile = shift || undef;

my $OUT_FILE = shift || undef;


my $REQUIRE_SAME_CNV_TYPE = 1;

if (!$truthFile || !$testFile || !$testFamFile || !$targetFile) {
  die "Usage: $0 truthFile testFile testFamFile targetFile [OUT_FILE]";
}

my $truthCNVs = readCNVfile($truthFile);
my $testCNVs = readCNVfile($testFile);

my $testSamples = readFamFile($testFamFile);

my $targets = readTargetFile($targetFile);

if (defined($OUT_FILE)) {
  open(STDOUT, ">$OUT_FILE") or die "can't open $OUT_FILE: $!";
}

print STDOUT "#SAMPLE= Sample ID of the 'truth' CNV\n";
print STDOUT "#CHR= Chromosome locus of the 'truth' CNV\n";
print STDOUT "#TYPE= Type of the 'truth' CNV\n";
print STDOUT "#TRUTH_CNV_START= Start position of the 'truth' CNV\n";
print STDOUT "#TRUTH_CNV_STOP= Stop position of the 'truth' CNV\n";
print STDOUT "#DIST= Distance between closest 'test' CNV(s) and the 'truth' CNV, in SAMPLE\n";
print STDOUT "#NUM_CLOSEST_CNVS= Number of 'test' CNV(s) in SAMPLE at distance DIST from the 'truth' CNV\n";
print STDOUT "#CLOSEST_CNVS= Loci of 'test' CNV(s) in SAMPLE at distance DIST from the 'truth' CNV\n";
print STDOUT "#OVERLAP= Cumulative overlap (in bases) between closest 'test' CNV(s) and the 'truth' CNV\n";
print STDOUT "#OVERLAP_FRAC= Fraction of the 'truth' CNV cumulatively overlapped (in bases) between closest 'test' CNV(s)\n";
print STDOUT "#CLOSEST_TARGET_DIST= Distance between closest exome target(s) and the 'truth' CNV\n";
print STDOUT "#NUM_CLOSEST_TARGETS= Number of exome target(s) at distance CLOSEST_TARGET_DIST from the 'truth' CNV\n";
print STDOUT "#CLOSEST_TARGETS= Loci of exome target(s) at distance CLOSEST_TARGET_DIST from the 'truth' CNV\n";
print STDOUT "#NORMALIZED_DIST= DIST - CLOSEST_TARGET_DIST\n";
print STDOUT "#NUM_CLOSEST_TARGETS_OVERLAPPED= How many of the 'truth' CNV-overlapping targets (CLOSEST_TARGETS) are overlapped by some 'test' CNV(s) that is closest to the 'truth' CNV (CLOSEST_CNVS)?  [i.e., target-level sensitivity to 'truth']\n";
print STDOUT "#OVERLAP_CLOSEST_TARGETS_FRAC= What fraction of the span of the targets closest to the 'truth' CNV are overlapped by some 'test' CNV(s)?\n";
print STDOUT "#NUM_TARGETS_OVERLAPPED= How many of the exome target(s) are overlapped by some 'test' CNV(s) that are closest to the 'truth' CNV [CLOSEST_CNVS]?  [i.e., target-level specificity to 'truth']\n";
print STDOUT "#TARGETS_OVERLAPPED= The exome target(s) overlapped the 'test' CNV(s) that are closest to the 'truth' CNV [CLOSEST_CNVS]\n";

print STDOUT "SAMPLE\tCHR\tTYPE\tTRUTH_CNV_START\tTRUTH_CNV_STOP\tDIST\tNUM_CLOSEST_CNVS\tCLOSEST_CNVS\tOVERLAP\tOVERLAP_FRAC\tCLOSEST_TARGET_DIST\tNUM_CLOSEST_TARGETS\tCLOSEST_TARGETS\tNORMALIZED_DIST\tNUM_CLOSEST_TARGETS_OVERLAPPED\tOVERLAP_CLOSEST_TARGETS_FRAC\tNUM_TARGETS_OVERLAPPED\tTARGETS_OVERLAPPED\n";

foreach my $sample (sort(keys(%{$truthCNVs}))) {
  # Don't compare with samples that could *never* be called:
  next if (!exists($testSamples->{$sample}));

  my %sampleCNVs = %{$truthCNVs->{$sample}};

  foreach my $chr (sort byChromosome (keys(%sampleCNVs))) {
    my %sampleChrCNVs = %{$sampleCNVs{$chr}};

    foreach my $type (keys(%sampleChrCNVs)) {
      my @TRUTH_sampleChrTypeCNVs = @{$sampleChrCNVs{$type}};

      foreach my $truthCNV (@TRUTH_sampleChrTypeCNVs) {
	my @TEST_sampleChrCNVs;
	if (exists($testCNVs->{$sample}) && exists($testCNVs->{$sample}{$chr})) { # otherwise, nothing to compare to:
	  if ($REQUIRE_SAME_CNV_TYPE) {
	    if (exists($testCNVs->{$sample}{$chr}{$type})) {
	      @TEST_sampleChrCNVs = @{$testCNVs->{$sample}{$chr}{$type}};
	    }
	  }
	  else {
	    foreach my $innerType (keys(%{$testCNVs->{$sample}{$chr}})) {
	      push(@TEST_sampleChrCNVs, @{$testCNVs->{$sample}{$chr}{$innerType}});
	    }
	  }
	}

	# Which test CNVs are closest to the truth CNV ?
	my ($dist, $closestCNVs) = calcNearestIntervals($truthCNV, \@TEST_sampleChrCNVs);
	my $overlap = Math::BigInt->bzero();
	for my $CNV (@{$closestCNVs}) {
	  $overlap->badd(calcOverlap($truthCNV, $CNV));
	}

	my $truthCNVlength = calcLength($truthCNV);
	my $overlapFrac = sprintf("%0.2f", Math::BigFloat->new($overlap)->bdiv($truthCNVlength));
	
	# Which are the targets closest to the truth CNV ?
	my ($closestTargetDist, $closestTargets) = (Math::BigInt->binf(), []);
	($closestTargetDist, $closestTargets) = calcNearestIntervals($truthCNV, \@{$targets->{$chr}}) if (exists($targets->{$chr}));

	# How many of the truth CNV-overlapping targets are overlapped by some test CNV that is closest to the truth CNV ?
	my $closestTargetsOverlap = Math::BigInt->bzero();
	my $closestTargetsLength = Math::BigInt->bzero();

	my $numClosestTargetsOverlapped = 0;
	foreach my $targInterval (@{$closestTargets}) {
	  $closestTargetsLength->badd(calcLength($targInterval));

	  foreach my $CNV (@{$closestCNVs}) {
	    if (calcMinDistance($CNV, $targInterval) == 0) {
	      $numClosestTargetsOverlapped++;
	      $closestTargetsOverlap->badd(calcOverlap($CNV, $targInterval));
	      last;
	    }
	  }
	}

	my $overlapClosestTargetsFrac = 0.00;
	$overlapClosestTargetsFrac = sprintf("%0.2f", Math::BigFloat->new($closestTargetsOverlap)->bdiv($closestTargetsLength)) if ($closestTargetsLength > 0);
	
	# How many of the targets are overlapped by some test CNV that is closest to the truth CNV ?
	my @overlappedTargets;
	foreach my $targInterval (@{$targets->{$chr}}) {
	  foreach my $CNV (@{$closestCNVs}) {
	    if (calcMinDistance($CNV, $targInterval) == 0) {
	      push(@overlappedTargets, $targInterval);
	      last;
	    }
	  }
	}

	my $normalizedDist = $dist->copy();
	$normalizedDist->bsub($closestTargetDist);

	my $output = "$sample\t$chr\t$type\t".$truthCNV->[0]."\t".$truthCNV->[1];

	$output .= "\t$dist\t".scalar(@{$closestCNVs})."\t";
	push(@{$closestCNVs}, ["NULL", "NULL"]) if (scalar(@{$closestCNVs}) == 0);
	foreach my $CNV (@{$closestCNVs}) {
	  $output .= $CNV->[0]."-".$CNV->[1].";";
	}
	$output .= "\t$overlap";
	$output .= "\t$overlapFrac";

	$output .= "\t$closestTargetDist\t".scalar(@{$closestTargets})."\t";
	push(@{$closestTargets}, ["NULL", "NULL"]) if (scalar(@{$closestTargets}) == 0);
	foreach my $target (@{$closestTargets}) {
	  $output .= $target->[0]."-".$target->[1].";";
	}

	$output .= "\t$normalizedDist";
	$output .= "\t$numClosestTargetsOverlapped";

	$output .= "\t$overlapClosestTargetsFrac";

	$output .= "\t".scalar(@overlappedTargets)."\t";
	push(@overlappedTargets, ["NULL", "NULL"]) if (scalar(@overlappedTargets) == 0);
	foreach my $CNV (@overlappedTargets) {
	  $output .= $CNV->[0]."-".$CNV->[1].";";
	}

	print STDOUT "$output\n";
      }
    }
  }
}



sub byChromosome {
  my $aIsNumber = isNumeric($a);
  my $bIsNumber = isNumeric($b);

  if ($aIsNumber && $bIsNumber) {
    return ($a <=> $b);
  }
  else {
    if ($aIsNumber) {
      return -1;
    }
    elsif ($bIsNumber) {
      return 1;
    }
  }

  return $a cmp $b;
}
