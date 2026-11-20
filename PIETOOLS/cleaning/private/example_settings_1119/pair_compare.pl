#!/usr/bin/perl
# pair_compare.pl HARNESS CONTROL TEST
# Compares two tier passes of the study harness row by row
# (HARNESS/results_<CONTROL>.csv against results_<TEST>.csv, keyed by
# library, id and tier): rows that pass in both, with an equal result (a
# threshold or value within 1e-3 relative; a verdict counts as equal), rows
# that pass in one only (listed), the ratio of decision variables
# control/test over the rows passing in both (geometric mean, median), and
# the summed variables of 'test first, control on failure' against the
# control alone. A pass is certified, converged or optimal. Every verdict
# helper returns one scalar (a regex match in list context returns an empty
# list on failure, which shifted the pairs of the first count of 11/19).
# MMP, 11/19/2026.
use strict; use warnings;
my ($H, $C, $T) = @ARGV;
sub rd { my ($f) = @_; my %r; open(my $h, '<', "$H/results_$f.csv") or die "no results_$f.csv"; <$h>;
  while (my $l = <$h>) { my @x = split /,/, $l; next if @x < 9; $r{"$x[0]:$x[1]:$x[2]"} = { st => $x[3], res => $x[4], ndv => $x[8] }; }
  close $h; return \%r; }
sub ok  { my $s = shift; return ($s =~ /^(?:certified|converged|optimal)$/) ? 1 : 0 }
sub num { my $v = shift; return (defined $v && $v =~ /^-?[\d.]+(?:e[-+]?\d+)?$/i) ? 1 : 0 }
my $A = rd($C); my $B = rd($T);
my ($n, $both, $eq, $c1, $c2) = (0, 0, 0, 0, 0); my (@onlyC, @onlyT, @ratio, @diff);
for my $k (sort keys %$A) {
  next unless $B->{$k}; my ($x, $y) = ($A->{$k}, $B->{$k}); $n++;
  my $ox = ok($x->{st}); my $oy = ok($y->{st});
  my $nx = $x->{ndv} =~ /^\d+$/ ? $x->{ndv} : 0; my $ny = $y->{ndv} =~ /^\d+$/ ? $y->{ndv} : 0;
  $c1 += $nx; $c2 += $oy ? $ny : $ny + $nx;
  if ($ox && !$oy) { push @onlyC, $k; next }
  if (!$ox && $oy) { push @onlyT, $k; next }
  next unless $ox && $oy; $both++;
  my ($rx, $ry) = ($x->{res}, $y->{res});
  if (num($rx) && num($ry) && $rx != 1 && abs($rx - $ry) > 1e-3 * abs($rx)) { push @diff, "$k $rx -> $ry" } else { $eq++ }
  push @ratio, $nx / $ny if $nx && $ny; }
my $g = 0; $g += log($_) for @ratio; $g = @ratio ? exp($g / @ratio) : 0; my @s = sort { $a <=> $b } @ratio;
printf "%s vs %s: rows %d | pass both %d (result equal %d) | %s only %d | %s only %d | variables %s/%s geo %.2f median %.2f | %s first, %s on failure: %.2fx the variables of %s\n",
  $C, $T, $n, $both, $eq, $C, scalar(@onlyC), $T, scalar(@onlyT), $C, $T, $g, @s ? $s[int(@s/2)] : 0, $T, $C, $c1 ? $c2 / $c1 : 0, $C;
print "  $C only: @onlyC\n" if @onlyC;
print "  $T only: @onlyT\n" if @onlyT;
print "  results differ: ", join('; ', @diff), "\n" if @diff;
