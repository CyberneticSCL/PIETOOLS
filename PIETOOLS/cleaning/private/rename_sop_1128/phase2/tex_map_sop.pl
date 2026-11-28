#!/usr/bin/perl
# TEX_MAP_SOP.PL TABLE FILE... - 11/28/2026, Phase 2 of the _sop rename (../PLAN.md): the map of
# Phase 1 on LaTeX notes, in place (adapted from ../../rename_copvar_opvar_1125/phase2/tex_map.pl).
# Each renamed row of TABLE (../rename_table.tsv: functions, tests, folders) as a whole name, with
# its underscores written "\_" (in \code{}) or "_" (in \file{}); the new name takes the same
# spelling. Kept: lpis_sopvar, the folder merged into lpi_programming_sopvar on 10/21/2026.
# Prints each changed line, old and new.
use strict; use warnings;
my $tab = shift or die "usage: TABLE FILE...\n";
open(my $T,'<',$tab) or die "$tab: $!"; <$T>; my @P;
while (<$T>) { s/\r?\n$//; my @f = split /\t/; push @P, [$f[0],$f[3]] if @f >= 4 && $f[3] ne $f[0] }
close $T;
@P = sort { length($b->[0]) <=> length($a->[0]) } @P;
for my $f (@ARGV) {
  open(my $h,'<:raw',$f) or die "$f: $!"; local $/; my $s = <$h>; close $h;
  my @L = split /(?<=\n)/, $s; my $n = 0;
  for my $i (0..$#L) {
    my $o = $L[$i]; my $x = $o;
    for my $p (@P) { my ($old,$new) = @$p;
      my $re = join '(\\\\?_)', map { quotemeta } split /_/, $old, -1;
      $x =~ s{(?<![A-Za-z0-9_\\])$re(?![A-Za-z0-9]|\\?_[A-Za-z0-9])}{
        my $m = $&; my $esc = $m =~ /\\_/; my $r = $new; $r =~ s/_/\\_/g if $esc; $r }ge }
    if ($x ne $o) { $n++; print "$f:", $i+1, ":\n  - $o  + $x"; $L[$i] = $x }
  }
  open($h,'>:raw',$f) or die "$f: $!"; print $h join('',@L); close $h;
  print STDERR "$f: $n lines changed\n";
}
