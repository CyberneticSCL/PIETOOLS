#!/usr/bin/perl
# fix_tiers_ref.pl PROPOSAL CORE: re-extracts CORE from PROPOSAL (as
# extract_core.pl), then replaces the marker (TIERS-REF) by a reference to
# the report's tier table in CORE and by the report's file name in
# PROPOSAL, and shortens the verbatim comments in CORE. SCRATCH, MMP,
# 11/16/2026.
use strict; use warnings;
my ($prop, $core) = @ARGV;
open(my $fh, '<', $prop) or die $!;  my $s = do { local $/; <$fh> };  close $fh;
my $a = index($s, '\section{One spatial variable}');
my $b = index($s, '\section{Requirements from the mixed-domain measurements');
die "markers not found" if $a < 0 || $b < 0 || $b <= $a;
my $c = substr($s, $a, $b - $a);
$c =~ s/\(TIERS-REF\)/(Table~\\ref{tab:tiers})/g;
open($fh, '>', $core) or die $!;  print $fh $c;  close $fh;
$s =~ s/\(TIERS-REF\)/(the report \\texttt{settings\\_report\\_2026-11-16.pdf}, Section~5)/g;
open($fh, '>', $prop) or die $!;  print $fh $s;  close $fh;
print "core: ", length($c), " chars; refs fixed\n";
