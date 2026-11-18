#!/usr/bin/perl
# extract_core.pl PROPOSAL OUT: writes to OUT the body of the settings
# proposal from "\section{One spatial variable}" up to (not including)
# "\section{Requirements from the mixed-domain measurements", i.e. the
# construction (1-D, conventions, N variables) and the settings row and
# reader, for \input into the report. SCRATCH, MMP, 11/16/2026.
use strict; use warnings;
my ($in, $out) = @ARGV;
open(my $fh, '<', $in) or die $!;  my $s = do { local $/; <$fh> };  close $fh;
my $a = index($s, '\section{One spatial variable}');
my $b = index($s, '\section{Requirements from the mixed-domain measurements');
die "markers not found" if $a < 0 || $b < 0 || $b <= $a;
my $core = substr($s, $a, $b - $a);
open($fh, '>', $out) or die $!;  print $fh $core;  close $fh;
print "core: ", length($core), " chars\n";
