#!/usr/bin/perl
# Names visible from the toolbox that are defined ONLY in cleaning/, orphaned_file/ or testing_codes/,
# and every toolbox code token (comments stripped, strings kept) that names one.
use strict; use warnings;
my ($root, $all, $out) = @ARGV;
my (%T, %Z, %Zwhere);
open(my $fa, '<', $all) or die; my @files = map { chomp; $_ } <$fa>; close $fa;
for my $f (@files) {
  next unless $f =~ /\.(m|p|mexw32|mexw64|mlx|mlapp|mdl|slx)$/;
  my $isZ = $f =~ m{^PIETOOLS/(cleaning|orphaned_file|testing_codes)/};
  my ($base) = $f =~ m{([^/]+)\.\w+$};
  my @names = ($base);
  if ($f =~ m{/@(\w+)/}) { push @names, $1; }            # class name visible
  if ($f =~ m{/\+(\w+)/}) { push @names, "$1.$base"; }
  if ($isZ) {
    next if $f =~ m{/private/} && $f !~ m{^PIETOOLS/(cleaning|orphaned_file|testing_codes)/private/[^/]+$} && 0;
    for my $n (@names) { $Z{$n} = 1; push @{$Zwhere{$n}}, $f; }
  } else { $T{$_} = 1 for @names; }
}
# Z names: a file under any private/ folder is visible only to its parent folder, which lies in Z too;
# keep them anyway (a toolbox call to such a name would fail, which is also a dependency worth seeing).
my @zonly = sort grep { !$T{$_} } keys %Z;
open(my $fz, '>', "$out/zonly.txt") or die; print $fz "$_\t", join(' ', @{$Zwhere{$_}}[0..($#{$Zwhere{$_}} < 1 ? $#{$Zwhere{$_}} : 1)]), "\n" for @zonly; close $fz;
my %zo = map { $_ => 1 } @zonly;
open(my $fh, '>', "$out/code_hits.tsv") or die;
for my $f (@files) {
  next if $f =~ m{^PIETOOLS/(cleaning|orphaned_file|testing_codes)/};
  next unless $f =~ /\.m$/;
  open(my $in, '<', "$root/$f") or next; my $ln = 0; my $blk = 0;
  while (my $l = <$in>) { $ln++; $l =~ s/\r?\n$//;
    if ($l =~ /^\s*%\{\s*$/) { $blk++; next; } if ($blk) { $blk-- if $l =~ /^\s*%\}\s*$/; next; }
    # strip a trailing comment: first % not inside a '...' string (transpose ' after identifier/bracket is not a quote)
    my ($code, $cmt) = ('', ''); my $inq = 0; my $prev = '';
    for (my $i = 0; $i < length $l; $i++) { my $ch = substr($l, $i, 1);
      if ($inq) { $code .= $ch; if ($ch eq "'") { if (substr($l,$i+1,1) eq "'") { $code .= "'"; $i++; } else { $inq = 0; } } next; }
      if ($ch eq '"') { my $j = index($l, '"', $i+1); $j = length($l)-1 if $j < 0; $code .= substr($l, $i, $j-$i+1); $i = $j; next; }
      if ($ch eq "'") { if ($prev =~ /[\w\)\]\}\.']/) { $code .= $ch; } else { $inq = 1; $code .= $ch; } $prev = $ch; next; }
      if ($ch eq '%') { $cmt = substr($l, $i); last; }
      $code .= $ch; $prev = $ch if $ch !~ /\s/; }
    for my $part (['code', $code], ['comment', $cmt]) {
      my ($kind, $txt) = @$part;
      while ($txt =~ /(?<![\w.])([A-Za-z]\w*(?:\.[A-Za-z]\w*)?)/g) { my $t = $1; my $t1 = $t; $t1 =~ s/\..*//;
        for my $cand ($t, $t1) { if ($zo{$cand}) { print $fh join("\t", $kind, $f, $ln, $cand, substr($l =~ s/^\s+//r, 0, 150)), "\n"; last; } } } }
  }
  close $in;
}
close $fh;
printf "Z names %d, toolbox names %d, Z-only names %d\n", scalar(keys %Z), scalar(keys %T), scalar(@zonly);
