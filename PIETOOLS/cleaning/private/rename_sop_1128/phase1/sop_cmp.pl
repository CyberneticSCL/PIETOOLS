#!/usr/bin/perl
# SOP_CMP.PL BASE NEW OUT DATE TABLE - checks 3c and 3a of Phase 1 (called by sop_check.sh),
# 11/28/2026; adapted from ../../rename_copvar_opvar_1125/phase1/ren_cmp.pl (no article rule).
# 3c: for each pair of OUT/pairs.tsv, NEW equals the sed-mapped BASE (OUT/expect) line by line,
#     except: the two D1 lines; a record line (header log, "(was)"), which equals BASE; runs of
#     blanks.
# 3a: in every live .m file of NEW, no old name of TABLE (lower, upper or capitalized) outside the
#     records.
use strict; use warnings; use File::Find;
my ($B,$N,$OUT,$DATE,$TAB) = @ARGV;
sub rd { open(my $h,'<',shift) or die $!; my @L = map { s/\r?\n$//r } <$h>; close $h; @L }
sub records { my @L=@_; my %r; my ($first,$seen)=(0,0);
  for my $i (0..$#L) { $r{$i}=1 if $L[$i] =~ /\(was\)/ }
  for my $i (0..$#L) { my $l=$L[$i];
    next if $l =~ /^\s*(?:function|classdef)\b/ || (!$seen && $i && $L[$i-1] =~ /\.\.\.\s*$/);
    if ($l =~ /^\s*%/) { $seen=1; next } next unless $l =~ /\S/; next if !$seen && $i<3;
    $first=$i+1; last }
  my $in=0;
  for my $i (0..($first ? $first-2 : -1)) { my $l=$L[$i];
    if ($l =~ /^%\s*(?:Initial coding|[A-Z]{2,3}(?:,? ?[A-Z]{2,3})*(?: \(CC\))?,? ?[ -] ?\d{1,2}\/\d{1,2}\/\d{2,4}|[A-Z]{2,3}(?:, [A-Z]{2,3})*, \d{1,2}\/\d{1,2}\/\d{2,4}|\d{1,2}\/\d{1,2}\/\d{2,4}\s*-)/) { $in=1; $r{$i}=1; next }
    if ($in && $l =~ /^%\s{2,}\S/) { $r{$i}=1; next } $in=0 }
  return \%r }
sub ws { my $s=shift; $s =~ s/[ \t]+/ /g; $s =~ s/ $//; $s }
my ($np,$bad,$nd1,$nrec,$nchg) = (0,0,0,0,0);
open(my $P,'<',"$OUT/pairs.tsv") or die $!;
while (my $p=<$P>) { chomp $p; my ($o,$n)=split /\t/,$p; $np++;
  my @O=rd("$B/$o"); my @S=rd("$OUT/expect/$n"); my @M=rd("$N/$n"); my @X;
  for (my $i=0; $i<@M; $i++) {
    if ($M[$i] eq "% MMP, $DATE: _sop dropped from the names of the executives, LPI" && $i+1<@M
        && $M[$i+1] =~ /^%\s+and settings functions \(2025 names where they exist\); names only\.$/) { $nd1++; $i++; next }
    push @X,$M[$i] }
  if (@X != @O || @S != @O) { print "3c LINES $n: base ".@O.", new ".@X.", sed ".@S."\n"; $bad++; next }
  my $r = records(@O);
  for my $i (0..$#O) { next if $X[$i] eq $O[$i] && $S[$i] eq $O[$i];
    if ($r->{$i}) { $nrec++ if $S[$i] ne $O[$i];
      if ($X[$i] ne $O[$i]) { print "3c RECORD $n:".($i+1)."\n"; $bad++ } next }
    $nchg++ if $X[$i] ne $O[$i];
    next if ws($X[$i]) eq ws($S[$i]);
    print "3c DIFF $n:".($i+1)."\n  new: $X[$i]\n  sed: $S[$i]\n"; $bad++ } }
close $P;
print "3c. comment identity: $np files, $nchg changed lines, $bad defects; D1 blocks $nd1; record lines kept $nrec\n";
# 3a: the old names of the table
open(my $T,'<',$TAB) or die $!; <$T>; my @old;
while (<$T>) { s/\r?\n$//; my @f = split /\t/; push @old, $f[0] if @f >= 4 && $f[3] ne $f[0] } close $T;
my %old = map { ($_ => 1, uc($_) => 1, ucfirst($_) => 1) } @old;
my ($nres,$krec)=(0,0);
find({ no_chdir=>1, wanted=>sub { my $f=$File::Find::name; return unless -f $f && $f =~ /\.m$/;
  (my $rel=$f) =~ s#^\Q$N\E/##; return if $rel =~ m#^(cleaning|orphaned_file|SOSTOOLS400)/#;
  my @L=rd($f); return unless grep { /_sop/i } @L;
  my $r=records(@L);
  for my $i (0..$#L) { next unless $L[$i] =~ /_sop/i;
    my @hit = grep { $old{$_} } ($L[$i] =~ /(\w+)/g); next unless @hit;
    if ($r->{$i}) { $krec++; next }
    $nres += @hit; print "3a RESIDUE $rel:".($i+1)." @hit\n" } }}, $N);
printf "3a. old names outside the records: %d; record lines that keep one: %d\n", $nres, $krec;
