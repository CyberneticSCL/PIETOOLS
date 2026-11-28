#!/usr/bin/perl
# MKQ_AFTER.PL TABLE < q_new.tsv > q_after.tsv - 11/28/2026: the Phase 3 queue of the container
# rename (rows at abfb0e83, the code of 91ea628e) for the tree after the _sop rename. Columns 3 to
# 5 (folder, command, extra folders) take the names of ../rename_table.tsv; the row ids stay, so
# that p3_cmp.pl pairs the rows; the export path of the P rows becomes C:\wt1128sop\after.
use strict; use warnings;
my $tab = shift or die "usage: TABLE < q_new.tsv\n";
open(my $T,'<',$tab) or die "$tab: $!"; <$T>; my %m;
while (<$T>) { s/\r?\n$//; my @f = split /\t/; $m{$f[0]} = $f[3] if @f >= 4 && $f[3] ne $f[0] } close $T;
my $re = join '|', map { quotemeta } sort { length($b) <=> length($a) } keys %m;
my $old = quotemeta('C:\wt1127\new\PIETOOLS'); my $new = 'C:\wt1128sop\after\PIETOOLS';
while (<STDIN>) { s/\r?\n$//; my @f = split /\t/, $_, -1;
  for my $k (2..$#f) { $f[$k] =~ s/(?<![A-Za-z0-9_])($re)(?![A-Za-z0-9_])/$m{$1}/g; $f[$k] =~ s/$old/$new/g }
  print join("\t", @f), "\n" }
