#!/usr/bin/perl
# H2 scan, 11/26/2026: lines of live .m files whose text the rename makes ambiguous.
#   W  names a retired 2025 class as a word: opvar, dopvar, opvar2d, dopvar2d (any case);
#   X  holds a token naming both a container and a 2025 class (opvar2copvar, COPVAR2OPVAR, ...);
#   Y  names poscopvar and poslpivar on one line;
#   Z  names a 2025 function whose name a container function takes (poslpivar, lpivar,
#      mat2opvar, rand_opvar) without the container name on the line.
# Kind: code (match in code, strings kept), trail (match in the comment after code),
# comment (comment-only line). Quote rule as h1scan.pl.
# Usage: perl h2scan.pl <PIETOOLS folder>  > h2_lines.tsv
use strict; use warnings;
my $root = shift;
sub split_line { my $l=shift; my ($o,$i,$n,$p)=('',0,length $l,'');
  while ($i<$n) { my $c=substr($l,$i,1);
    return ($o, substr($l,$i)) if $c eq '%' || substr($l,$i,3) eq '...';
    if ($c eq '"' || ($c eq "'" && $p !~ /[\w\)\]\}\.']/)) { my $q=$c; my $j=$i+1;
      while ($j<$n) { my $d=substr($l,$j,1); if ($d eq $q) { if ($j+1<$n && substr($l,$j+1,1) eq $q) { $j+=2; next } last } $j++ }
      $o.=substr($l,$i,$j-$i+1); $i=$j+1; $p=$q; next }
    $o.=$c; $p=$c; $i++ }
  return ($o, '') }
sub pats { my $t=shift; my @p;
  push @p,'W' if $t =~ /\b(?:opvar|dopvar|opvar2d|dopvar2d)\b/i;
  my $x=0; while ($t =~ /(\w+)/g) { my $w=$1;
    if ($w =~ /c(?:d)?opvar/i) { (my $r=$w) =~ s/c(?:d)?opvar//ig; $x=1 if $r =~ /opvar/i } }
  push @p,'X' if $x;
  push @p,'Y' if $t =~ /poscopvar/i && $t =~ /poslpivar/i;
  push @p,'Z' if $t =~ /\b(?:poslpivar|lpivar|mat2opvar|rand_opvar)\b/i && $t !~ /c(?:d)?opvar/i;
  return join('',@p) }
my @files = split /\n/, `git -C "$root" ls-files -- "*.m"`;
for my $f (@files) { next if $f =~ m#^(cleaning|orphaned_file|SOSTOOLS400)/#;
  open(my $h,'<',"$root/$f") or next; my $k=0;
  while (my $l=<$h>) { $k++; $l=~s/\r?\n$//; my ($code,$cmt)=split_line($l);
    my $pc=pats($code); my $pm=pats($cmt); next unless $pc || $pm;
    my $kind = $pc ? 'code' : ($code =~ /\S/ ? 'trail' : 'comment');
    my %u; my $pp=join('',grep { !$u{$_}++ } split //, $pc.$pm);
    (my $t=$l) =~ s/^\s+//; $t = substr($t,0,170);
    my $area = $f =~ m#^testing_codes/# ? 'test' : 'toolbox';
    print "$area\t$kind\t$pp\t$f:$k\t$t\n" } close $h }
