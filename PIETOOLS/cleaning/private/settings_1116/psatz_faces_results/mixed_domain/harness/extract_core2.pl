#!/usr/bin/perl
# extract_core2.pl PROPOSAL CORE CORE_B: splits the shared construction text
# of the proposal into CORE (from "\section{One spatial variable}" up to
# "\section{The present tiers in the row}") and CORE_B (from there up to
# "\section{Requirements from the mixed-domain measurements"), so the
# report can place its proposed tiers between the settings and the present
# tiers (maintainer, 11/17: the "what" before the "why"). Replaces
# extract_core.pl. SCRATCH, MMP, 11/17/2026.
use strict; use warnings;
my ($in, $out, $outb) = @ARGV;
open(my $fh, '<', $in) or die $!;  my $s = do { local $/; <$fh> };  close $fh;
my $a = index($s, '\section{One spatial variable}');
my $m = index($s, '\section{The present tiers in the row}');
my $b = index($s, '\section{Requirements from the mixed-domain measurements');
die "markers not found" if $a < 0 || $m < 0 || $b < 0 || !($a < $m && $m < $b);
open($fh, '>', $out) or die $!;   print $fh substr($s, $a, $m - $a);  close $fh;
open($fh, '>', $outb) or die $!;  print $fh substr($s, $m, $b - $m);  close $fh;
print "core: ", $m - $a, " chars; core_b: ", $b - $m, " chars\n";
