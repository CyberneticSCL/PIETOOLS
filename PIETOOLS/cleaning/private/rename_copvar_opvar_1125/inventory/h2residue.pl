#!/usr/bin/perl
# Labels each line of an h2scan.pl output by the reason the H2 pass keeps it (11/26/2026):
#   log       part of a dated header change-log entry (header logs stay as written);
#   was       a (was) line: quoted old code;
#   dated     a comment line with a trailing date marker, or inside a dated BEGIN ... END note;
#   quoted    commented-out old code (the line declares or fills a 2025 object);
#   qualified names the class as retired, stock, legacy or 2025;
#   baseline  a dated baseline folder under testing_codes/private/baselines;
#   tests_1d  the 2025-API low-rank tests (open case C8);
#   name2d    names only opvar2d/dopvar2d, which no container takes;
#   generic   lpivar as the declaration of a free LPI variable;
#   retired   a retired function name that holds a container token (X, Y);
#   design    transl_mod lines written for the opvar of after the rename (C1);
#   code      code (strings kept) naming a 2025 class: H1, decided in REVIEW_PHASE0.md.
# Unlabelled lines are printed with the label '-'.
# Usage: perl h2residue.pl <PIETOOLS folder> < h2_lines.tsv
use strict; use warnings;
my $root = shift; my (%hdr,%blk);
sub header_info { my $f=shift; return $hdr{$f} if $hdr{$f};
  open(my $h,'<',"$root/$f") or return {}; my @L=map { s/\r?\n$//r } <$h>; close $h;
  my ($in,$k,%log,%dblk)=(0,0);
  # first code line: the function line and its continuations do not count, as many files
  # put the help and the log after the function line
  my $first = 0; for my $i (0..$#L) { next if $L[$i] =~ /^\s*function\b/ || ($i && $L[$i-1] =~ /\.\.\.\s*$/);
    if ($L[$i] !~ /^\s*%/ && $L[$i] =~ /\S/) { $first=$i+1; last } }
  for my $i (0..$#L) { my $l=$L[$i];
    if ($first && $i+1 < $first) {
      if ($l =~ /^%\s*(?:Initial coding|[A-Z]{2,3}(?:,? ?[A-Z]{2,3})*(?: \(CC\))?,? ?[ -] ?\d{1,2}\/\d{1,2}\/\d{2,4}|[A-Z]{2,3}(?:, [A-Z]{2,3})*, \d{1,2}\/\d{1,2}\/\d{2,4}|\d{1,2}\/\d{1,2}\/\d{2,4}\s*-)/) { $in=1; $log{$i+1}=1; next }
      if ($in && $l =~ /^%\s{2,}\S/) { $log{$i+1}=1; next }
      $in=0 }
    $in=0 }
  my $d=0; for my $i (0..$#L) { my $l=$L[$i];
    $d=1 if $l =~ /^\s*%.*(?:\bBEGIN\b.*\d{1,2}\/\d{1,2}\/\d{4}|\d{1,2}\/\d{1,2}\/\d{4}.*\bBEGIN\b)/ || $l =~ /Start of changes/;
    $dblk{$i+1}=1 if $d;
    $d=0 if $l =~ /^\s*%.*(?:\bEND\b.*\d{1,2}\/\d{1,2}\/\d{4}|\d{1,2}\/\d{1,2}\/\d{4}.*\bEND\b)/ || $l =~ /End of changes/ || ($d && $l !~ /^\s*%/ && $l =~ /\S/) }
  return $hdr{$f} = {log=>\%log, dblk=>\%dblk} }
while (my $r=<STDIN>) { chomp $r; my ($area,$kind,$pat,$loc,$t)=split /\t/,$r,5; my ($f,$n)=$loc=~/^(.*):(\d+)$/;
  my $hi = header_info($f); my $lab='-';
  my $prev = ''; if (open(my $h,'<',"$root/$f")) { my $k=0; while (my $l=<$h>) { $k++; if ($k==$n-1) { $prev=$l; last } } close $h }
  my %reviewed = map { $_=>1 } ('PIESIM/PIESIM_routines/PIESIM_Setup/transl_mod.m:18',  # C1: the opvar of after the rename
    'PIETOOLS_demos/cuadmm/cuadmm_path.m:26',   # stays true: a shadow named poslpivar shadows the renamed one too
    'PIETOOLS_demos/cuadmm/cuadmm_path.m:81',   # a measurement dated 2026-09-25 in the same sentence
    'PIETOOLS_demos/cuadmm/cuadmm_path.m:98',   # the body note "CC, 11/07/2026: ..." (dated) of line 96
    'PIETOOLS_demos/cuadmm/cuadmm_path.m:99');  # commented-out old code under that note
  if    ($reviewed{$loc})                                           { $lab='reviewed' }
  elsif ($f =~ m#^testing_codes/private/baselines/20#)            { $lab='baseline' }
  elsif ($f =~ m#/tests_1d/#)                                       { $lab='tests_1d' }
  elsif ($hi->{log}{$n})                                            { $lab='log' }
  elsif ($t =~ /\(was\)\s*$/)                                       { $lab='was' }
  elsif ($f =~ /transl_mod\.m$/ && $n >= 74)                        { $lab='design' }
  elsif ($kind eq 'code')                                           { $lab='code' }
  elsif ($t =~ /%\s*[A-Z]{2,3}(?:, ?[A-Z]{2,3})*(?: \(CC\))?, \d{2}\/\d{2}\/\d{4}\s*$/ || $hi->{dblk}{$n}) { $lab='dated' }
  elsif ($t =~ /^%+\s*(?:%\s*)?\s*(?:opvar|dopvar)\s+\w/i)          { $lab='quoted' }
  elsif ($t =~ /\b(?:retired|stock|legacy|2025)\b/i || $prev =~ /\b(?:retired|stock|legacy|2025)\s*$/i) { $lab='qualified' }
  elsif ($pat !~ /[XYZ]/ && $t !~ /\b(?:opvar|dopvar)\b/i)          { $lab='name2d' }
  elsif ($pat =~ /^[XY]+$/)                                         { $lab='retired' }
  elsif ($pat eq 'Z' && $t =~ /\blpivar\b/ && $t !~ /\bposlpivar\b|\bmat2opvar\b/) { $lab='generic' }
  print "$lab\t$area\t$kind\t$pat\t$loc\t$t\n" }
