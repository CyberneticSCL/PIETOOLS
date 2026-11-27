#!/usr/bin/perl
# REN_CMP.PL BASE NEW OUT DATE - checks 3c and 3a of Phase 1 (called by ren_check.sh), 11/26/2026.
# 3c: for each pair of OUT/pairs.tsv, NEW equals the sed-mapped BASE (OUT/expect) line by line,
#     except: the two D1 lines; a record line (header log, "(was)"), which equals BASE; runs of
#     blanks; an "an" that ends a comment line before a mapped cop... token.
# 3a: in every live .m file of NEW, a token naming a container outside the records is a retired
#     name.
use strict; use warnings; use File::Find;
my ($B,$N,$OUT,$DATE) = @ARGV;
my @RET = qw(op2copvar_sop pie2copvar opvar2copvar copvar2opvar opvar2d2copvar cg_copvar_pd
             test_cg_copvar_nd test_copvar_kyp test_poscopvar_vs_poslpivar_2d);
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
my ($np,$bad,$nd1,$nrec,$nart,$nchg) = (0,0,0,0,0,0);
open(my $P,'<',"$OUT/pairs.tsv") or die $!;
while (my $p=<$P>) { chomp $p; my ($o,$n)=split /\t/,$p; $np++;
  my @O=rd("$B/$o"); my @S=rd("$OUT/expect/$n"); my @M=rd("$N/$n"); my @X;
  for (my $i=0; $i<@M; $i++) {
    if ($M[$i] eq "% MMP, $DATE: copvar, cdopvar renamed opvar, dopvar, with the" && $i+1<@M
        && $M[$i+1] =~ /^%\s+functions of the rename table; names only\.$/) { $nd1++; $i++; next }
    push @X,$M[$i] }
  if (@X != @O || @S != @O) { print "3c LINES $n: base ".@O.", new ".@X.", sed ".@S."\n"; $bad++; next }
  my $r = records(@O);
  for my $i (0..$#O) { next if $X[$i] eq $O[$i] && $S[$i] eq $O[$i];
    if ($r->{$i}) { $nrec++ if $S[$i] ne $O[$i];
      if ($X[$i] ne $O[$i]) { print "3c RECORD $n:".($i+1)."\n"; $bad++ } next }
    $nchg++ if $X[$i] ne $O[$i];
    next if ws($X[$i]) eq ws($S[$i]);
    (my $y=$X[$i]) =~ s/\b([Aa])n(\s*)$/$1$2/;
    if (ws($y) eq ws($S[$i]) && $i<$#O && $O[$i+1] =~ /^\s*%+\s*['"`(]?(?:copvar|copquadvar|copvar_space_list|copvar2tensopmat)\b/i) { $nart++; next }
    print "3c DIFF $n:".($i+1)."\n  new: $X[$i]\n  sed: $S[$i]\n"; $bad++ } }
close $P;
print "3c. comment identity: $np files, $nchg changed lines, $bad defects; D1 blocks $nd1; record lines kept $nrec; article at a line end $nart\n";
my ($nres,%kept)=(0);
find({ no_chdir=>1, wanted=>sub { my $f=$File::Find::name; return unless -f $f && $f =~ /\.m$/;
  (my $rel=$f) =~ s#^\Q$N\E/##; return if $rel =~ m#^(cleaning|orphaned_file|SOSTOOLS400)/#;
  my @L=rd($f); return unless grep { /c(?:d)?opvar|copquadvar|apply_cop|poslpivar_sop/i } @L;
  my $r=records(@L);
  for my $i (0..$#L) { next unless $L[$i] =~ /c(?:d)?opvar|copquadvar|apply_cop|poslpivar_sop/i;
    if ($r->{$i}) { $kept{record}++; next }
    for my $t ($L[$i] =~ /(\w+)/g) { next unless $t =~ /c(?:d)?opvar|copquadvar|apply_cop|poslpivar_sop/i;
      if (grep { index(lc $t,$_) >= 0 } @RET) { $kept{retired}++ } else { $nres++; print "3a RESIDUE $rel:".($i+1)." $t\n" } } } }}, $N);
printf "3a. residue outside the records: %d; kept: %d record lines, %d retired-name tokens\n", $nres, $kept{record}//0, $kept{retired}//0;
