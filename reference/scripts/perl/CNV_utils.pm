package CNV_utils;

use strict;
use base 'Exporter';

use Math::BigInt;

our @EXPORT = qw(readCNVfileFullInfo readCNVfileHeader readCNVfile readFamFile readTargetFile parseIntervalString getIntervalString readGeneListFile readTranscriptsFile readTriosFile readTrioParentsFile readTabDelimitedFile calcMinDistance calcLength calcOverlap getOverlappingCNVsinSample min max intervalLength dist calcNearestIntervals isNumeric);
our @EXPORT_OK = qw();

#Variable constants to export:


#General subs:
sub readCNVfileFullInfo {
  my ($file) = @_;

  open FILE, $file or die "$!: $file";

  # Read the header fields:
  my %HEADER_FIELDS;

  my $line = <FILE>;
  chomp($line);
  # Split on white-space:
  my @headerFields = split(' ', $line);
  for (my $i = 0; $i < @headerFields; $i++) {
    $HEADER_FIELDS{$headerFields[$i]} = $i;
  }

  my @CNVlines;
  while ($line = <FILE>) {
    chomp($line);
    # Split on white-space:
    my @fields = split(' ', $line);

    my %dataLine;
    for (my $i = 0; $i < scalar(@headerFields); $i++) {
      my $col = $headerFields[$i];
      $dataLine{$col} = $fields[$i];
    }
    push(@CNVlines, \%dataLine);
  }

  close(FILE);

  return (\@CNVlines, \@headerFields);
}

sub readCNVfileHeader {
  my ($file) = @_;

  open FILE, $file or die "$!: $file";

  my $line = <FILE>;
  chomp($line);
  # Split on white-space:
  my @headerFields = split(' ', $line);

  close(FILE);

  return (\@headerFields);
}

sub readCNVfile {
  my ($file) = @_;

  my ($CNVlines, $headerFields) = readCNVfileFullInfo($file);

  my %CNVs;
  foreach my $line (@{$CNVlines}) {
    my $sample = $line->{"FID"};
    my $chr = $line->{"CHR"};
    my $start = $line->{"BP1"};
    my $stop = $line->{"BP2"};

    my $CNVtype = "";
    my $copy = $line->{"TYPE"};
    if ($copy < 2) {
      $CNVtype = "del";
    }
    elsif ($copy > 2) {
      $CNVtype = "dup";
    }

    push(@{$CNVs{$sample}{$chr}{$CNVtype}}, [$start, $stop, $line]) if ($CNVtype ne '');
  }

  return \%CNVs;
}

sub readFamFile {
  my ($file) = @_;

  my %samples;

  open FILE, $file or die "$!: $file";

  my $line;
  while ($line = <FILE>) {
    chomp($line);
    # Split on white-space:
    my @fields = split(' ', $line);

    my $sample = $fields[0];
    $samples{$sample} = 1;
  }

  close(FILE);

  return \%samples;
}

sub readTargetFile {
  my ($file) = @_;

  my %targets;

  open FILE, $file or die "$!: $file";

  my $line;
  while ($line = <FILE>) {
    chomp($line);

    my ($chr, $interval) = parseIntervalString($line);
    push(@{$targets{$chr}}, $interval);
  }

  close(FILE);

  return \%targets;
}

sub parseIntervalString {
  my ($str) = @_;

  my ($chr, $interval);
  if ($str =~ m/^(.+):(.+)-(.+)$/) {
    $chr = $1;
    my $start = $2;
    my $stop = $3;

    $interval = [$start, $stop];
  }
  elsif ($str =~ m/^(.+):(.+)$/) {
    $chr = $1;
    my $start = $2;

    $interval = [$start, $start];
  }
  else {
    die("Invalid interval: $str");
  }

  return ($chr, $interval);
}

sub getIntervalString {
  my ($chr, $interval) = @_;

  return $chr.":".$interval->[0]."-".$interval->[1];
}

sub readGeneListFile {
  my ($file) = @_;

  my %genes;

  open FILE, $file or die "$!: $file";

  my $line;
  while ($line = <FILE>) {
    chomp($line);
    # Split on white-space:
    my @fields = split(' ', $line);

    if (scalar(@fields) == 4) {
      my $chr = $fields[0];
      push(@{$genes{$chr}}, [$fields[1], $fields[2], $fields[3]]);
    }
    else {
      die("Invalid line in gene file $file: $line");
    }
  }

  close(FILE);

  return \%genes;
}

sub readTranscriptsFile {
  my ($file) = @_;

  my %transcripts;

  open FILE, $file or die "$!: $file";

  my $line;
  while ($line = <FILE>) {
    chomp($line);
    # Split on white-space:
    my @fields = split(' ', $line);

    if (scalar(@fields) == 6) {
      my $t = $fields[0];

      $transcripts{$t}{'gene'} = $fields[1];
      $transcripts{$t}{'chr'} = $fields[2];
      $transcripts{$t}{'span'} = [$fields[3], $fields[4]];
      $transcripts{$t}{'strand'} = $fields[5];
    }
    else {
      die("Invalid line in gene file $file: $line");
    }
  }

  close(FILE);

  return \%transcripts;
}

sub readTriosFile {
  my ($file) = @_;

  my %trios;

  open FILE, $file or die "$!: $file";

  my $line;
  while ($line = <FILE>) {
    chomp($line);
    # Split on white-space:
    my @fields = split(' ', $line);

    if (scalar(@fields) == 3) {
      my $proband = $fields[0];
      $trios{$proband} = {'father' => $fields[1], 'mother' => $fields[2]};
    }
    else {
      die("Invalid line in trios file $file: $line");
    }
  }

  close(FILE);

  return \%trios;
}

sub readTrioParentsFile {
  my ($file) = @_;

  my %trioParents;

  open FILE, $file or die "$!: $file";

  my $line;
  while ($line = <FILE>) {
    chomp($line);
    # Split on white-space:
    my @fields = split(' ', $line);

    if (scalar(@fields) == 3) {
      my $proband = $fields[0];
      my $father = $fields[1];
      my $mother = $fields[2];

      $trioParents{$father} = ['father', $proband];
      $trioParents{$mother} = ['mother', $proband];
    }
    else {
      die("Invalid line in trioParents file $file: $line");
    }
  }

  close(FILE);

  return \%trioParents;
}

sub readTabDelimitedFile {
  my ($file) = @_;

  open FILE, $file or die $!;

  my @lines;
  while (my $line = <FILE>) {
    chomp($line);
    # Split on white-space:
    my @fields = split(' ', $line);

    push(@lines, \@fields);
  }

  close(FILE);

  return \@lines;
}

sub calcMinDistance {
  my ($interval1, $interval2) = @_;

  my $start1 = $interval1->[0];
  my $stop1 = $interval1->[1];

  my $start2 = $interval2->[0];
  my $stop2 = $interval2->[1];

  if ($stop1 < $start2) { # 1 is before 2
    return $start2 - $stop1;
  }
  elsif ($stop2 < $start1) { # 2 is before 1
    return $start1 - $stop2;
  }

  # The segments overlap:
  return 0;
}

sub calcLength {
  my ($interval) = @_;

  return ($interval->[1] - $interval->[0] + 1);
}

sub calcOverlap {
  my ($interval1, $interval2) = @_;

  if (calcMinDistance($interval1, $interval2) > 0) {
    return 0; # no overlap if distance >= 1
  }

  my $start1 = $interval1->[0];
  my $stop1 = $interval1->[1];

  my $start2 = $interval2->[0];
  my $stop2 = $interval2->[1];

  return dist(max($start1, $start2), min($stop1, $stop2));
}

sub getOverlappingCNVsinSample {
  my ($sample, $type, $chr, $cnv, $CNVS, $REQUIRE_SAME_CNV_TYPE) = @_;

  my @checkSampleChrCNVs;
  if (exists($CNVS->{$sample}) && exists($CNVS->{$sample}{$chr})) { # otherwise, nothing to compare to:
    if ($REQUIRE_SAME_CNV_TYPE) {
      if (exists($CNVS->{$sample}{$chr}{$type})) {
	@checkSampleChrCNVs = @{$CNVS->{$sample}{$chr}{$type}};
      }
    }
    else {
      foreach my $innerType (keys(%{$CNVS->{$sample}{$chr}})) {
	push(@checkSampleChrCNVs, @{$CNVS->{$sample}{$chr}{$innerType}});
      }
    }
  }

  my ($dist, $closestCNVs) = calcNearestIntervals($cnv, \@checkSampleChrCNVs);
  if ($dist == 0) {
    return $closestCNVs;
  }
  else {
    return [];
  }
}

sub min {
  my ($a, $b) = @_;

  if ($a < $b) {
    return $a;
  }

  return $b;
}

sub max {
  my ($a, $b) = @_;

  if ($a > $b) {
    return $a;
  }

  return $b;
}

sub intervalLength {
  my ($interval) = @_;

  return dist($interval->[0], $interval->[1]);
}

sub dist {
  my ($start, $stop) = @_;

  return ($stop - $start + 1);
}

sub calcNearestIntervals {
  my ($interval, $intervalsArray) = @_;

  my $dist = Math::BigInt->binf();

  my @closestIntervals;
  foreach my $testInterval (@{$intervalsArray}) {
    my $minDist = Math::BigInt->new(calcMinDistance($interval, $testInterval));
    my $diff = $minDist->bcmp($dist);

    if ($diff < 0) {
      $dist = $minDist;
      @closestIntervals = ();
    }
    if ($diff <= 0) {
      push(@closestIntervals, $testInterval);
    }
  }

  return ($dist, \@closestIntervals);
}

sub isNumeric {
  my ($val) = @_;

  return ($val =~ m/^(\d+\.?\d*|\.\d+)$/);
}






1;
__END__
