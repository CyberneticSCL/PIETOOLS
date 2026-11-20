#!/usr/bin/perl
# fallback.pl HARNESS_DIR DIR VERIFY_CSV
# For each example with a row that fails its verification (VERIFY_CSV of
# verify_check.pl, ok = 0), the cheapest setting of the tier passes and the
# 11/18-11/19 studies that reached the best result on every row of the
# example in the same pass (full bisections, values, verdicts: measured, not
# inferred). Prints the example, the pass, its settings and form, and the
# largest decision-variable count over its rows. MMP, 11/19/2026.
use strict; use warnings;
my ($H, $D, $V) = @ARGV;  my $RT = 1e-3;
my %TAG = (   # pass => [label, st changes from lpisettings_sop(tier)]
  M1L => ['light', ''], B2dD0 => ['light, dD 0', 'st.lpi.ineq.dD = 0;'], H1ctl => ['heavy', ''],
  B3q522 => ['heavy, deg_Q{1} [5 2 2]', 'st.sop.deg_Q{1} = [5 2 2];'], B3q622 => ['heavy, deg_Q{1} [6 2 2]', 'st.sop.deg_Q{1} = [6 2 2];'],
  C1ext => ['extreme', ''], C2ext => ['extreme, stripped 1-D entries', 'stripped1d'],
  C3ext => ['extreme, stripped 1-D entries, dD 0', 'stripped1d; st.lpi.ineq.dD = 0;'],
  D1P => ['light', ''], D1Ph => ['heavy', ''], D1Ppde => ['light', ''], D1Ppdeh => ['heavy', ''],
  GallP => ['light', ''], GallPh => ['heavy', ''], GallQ => ['light', ''], GallQh => ['heavy', ''],
  XSP => ['stripped', ''], XSPd0 => ['stripped, dD 0', 'st.lpi.ineq.dD = 0;'], XLPd0 => ['light, dD 0', 'st.lpi.ineq.dD = 0;'],
  XSQ => ['stripped', ''], XSQd0 => ['stripped, dD 0', 'st.lpi.ineq.dD = 0;'], XVH => ['veryheavy', ''], B4ctl => ['veryheavy', ''] );
sub isnum { my $v = shift; return defined $v && $v =~ /^-?[0-9.]+(?:[eE][-+]?\d+)?$/ }
sub csv { my ($f) = @_; my @r; open(my $h, '<', $f) or return (); my $hd = <$h>; $hd =~ s/\r?\n$//; my @c = split /,/, $hd;
  while (my $l = <$h>) { $l =~ s/\r?\n$//; my @x = ($l =~ /("(?:[^"]|"")*"|[^,]*)(?:,|$)/g); s/^"(.*)"$/$1/ for @x; my %e; @e{@c} = @x; push @r, \%e } return @r }
my %ch = map { ("$_->{library}:$_->{id}") => $_ } csv("$D/phaseB_choices.csv");
my %bad = map { $_->{key} => 1 } grep { $_->{ok} eq '0' } csv($V);
my @plan = csv("$D/phaseC_plan.csv");
my %res;
for my $t (keys %TAG) { open(my $h, '<', "$H/results_$t.csv") or next; <$h>;
  while (my $l = <$h>) { $l =~ s/\r?\n$//; my @x = split /,/, $l, 10; next if @x < 10;
    $res{$t}{"$x[0]:$x[1]"} = { tier => $x[2], st => $x[3], res => $x[4], ndv => $x[8], note => $x[9] }; } close $h; }
sub reach { my ($k, $e) = @_; my $c = $ch{$k} or return 0; my $kd = $c->{kind}; my $b = $c->{best};
  my ($ex) = $e->{note} =~ /exst=(\w+)/; my $eg = !defined $ex || $ex =~ /optimal|inaccurate/;
  if ($kd eq 'bisect') { return $e->{st} =~ /^(converged|max_solves|hi_pass)$/ && isnum($e->{res}) && $e->{res} >= $b - $RT * abs($b) }
  if ($kd eq 'value')  { return $e->{st} eq 'certified' && $eg && isnum($e->{res}) && $e->{res} <= $b + $RT * abs($b) }
  return $e->{st} eq 'certified' && $eg; }
for my $p (@plan) { my @rows = split / /, $p->{rows}; next unless grep { $bad{$_} } @rows;
  my @cand;
  for my $t (sort keys %TAG) { my $all = 1; my ($mx, %fm) = (0);
    for my $k (@rows) { my $e = $res{$t}{$k}; unless ($e && reach($k, $e)) { $all = 0; last }
      $mx = $e->{ndv} if isnum($e->{ndv}) && $e->{ndv} > $mx; my ($f) = $e->{note} =~ m{_sop/(\w+)}; $fm{$f // '-'}++; }
    push @cand, [$t, $mx, join('/', sort keys %fm), $res{$t}{$rows[0]}{tier}] if $all; }
  my ($best) = sort { $a->[1] <=> $b->[1] } @cand;
  printf "%-44s rows %-48s -> %s\n", $p->{group}, join(' ', @rows),
    $best ? sprintf("%s (%s; tier %s; changes '%s'; forms %s), max %d variables", $best->[0], $TAG{$best->[0]}[0], $best->[3], $TAG{$best->[0]}[1], $best->[2], $best->[1]) : 'no common pass';
}
