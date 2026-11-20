#!/usr/bin/perl
# tiers_vs_stock.pl HARNESS_DIR BASE_1028
# The present settings tiers (lpisettings_sop, container executives, current
# tree: the 11/18-11/19 control passes) against the stock tiers (stock
# executives and settings, the baseline of 10/28/2026, stock_results.csv),
# row by row on the 138 rows the baseline holds. Verdicts by the harness
# gate on both sides (|fr-1| <= 0.3, no numerr, pinf, dinf; bisections by
# the harness gate per step); thresholds (larger is better) and values
# (gains, norms: smaller is better) compared at 1e-3 relative; size by the
# decision variables of the last solve. MMP, 11/19/2026.
use strict; use warnings;
my ($H, $B) = @ARGV;
sub rd { my ($f, $tier) = @_; my %r; open(my $h, '<', $f) or return \%r; <$h>;
  while (my $l = <$h>) { $l =~ s/\r?\n$//; my @x = split /,/, $l, 10; next if @x < 10;
    next if defined $tier && $x[2] ne $tier;
    $r{"$x[0]:$x[1]"} = { st=>$x[3], res=>$x[4], n=>$x[6], sec=>$x[7], ndv=>$x[8], note=>$x[9] }; } close $h; return \%r; }
my $stock = { map { my $t = $_; ($t => rd("$B/stock_results.csv", $t)) } qw(stripped light heavy veryheavy extreme) };
# present tiers, current tree: 1-D and 2-D sources per tier
my %src = (
  light     => [rd("$H/results_M1L.csv"), rd("$H/results_B1c0.csv")],
  heavy     => [rd("$H/results_H1ctl.csv"), rd("$H/results_H2ctl.csv", 'heavy')],
  stripped  => [rd("$H/results_XSQ.csv"), rd("$H/results_X2SQ.csv"), rd("$H/results_E1str.csv")],
  extreme   => [rd("$H/results_C1ext.csv"), rd("$H/results_C1ext2.csv")],
  veryheavy => [rd("$H/results_B4ctl.csv"), rd("$H/results_E3vh.csv"), rd("$H/results_all_1118b.csv", 'veryheavy')] );
sub isnum { my $v = shift; return defined $v && $v =~ /^-?[0-9.]+(?:[eE][-+]?\d+)?$/ }
sub bis { $_[0] =~ /^(converged|max_solves|lo_fail|hi_pass)$/ }
sub cert { my ($e) = @_; return 0 unless $e;
  return (bis($e->{st}) ? ($e->{st} ne 'lo_fail' && isnum($e->{res})) : $e->{st} eq 'certified') ? 1 : 0; }
sub isval { my $e = shift; return $e && $e->{note} =~ /Hinf_|H2_/ && !bis($e->{st}) }
printf "%-9s %4s | %5s %6s %6s %6s | %-24s | %-24s | %s\n", 'tier', 'rows', 'both', 'stock', 'new', 'neith',
  'thresholds new vs stock', 'values new vs stock', 'variables new/stock (median, sum)';
for my $t (qw(stripped light heavy veryheavy extreme)) {
  my %c; my (@ratio, $s0, $s1) = (); $s0 = $s1 = 0; my (%th, %va, @lostk, @gaink);
  for my $k (sort keys %{$stock->{$t}}) {
    my $s = $stock->{$t}{$k}; my ($n) = grep { defined } map { $_->{$k} } @{$src{$t}};
    next unless $n;      # a row the present pass did not run (time projection, size gate)
    $c{rows}++; my ($cs, $cn) = (cert($s), cert($n));
    if ($cs && $cn) { $c{both}++ } elsif ($cs) { $c{stock}++; push @lostk, $k } elsif ($cn) { $c{new}++; push @gaink, $k } else { $c{neither}++ }
    if ($cs && $cn && bis($s->{st}) && isnum($s->{res}) && isnum($n->{res})) {
      my $d = ($n->{res} - $s->{res}) / (abs($s->{res}) || 1);
      $th{$d > 1e-3 ? 'higher' : $d < -1e-3 ? 'lower' : 'same'}++; }
    if ($cs && $cn && isval($s) && isnum($s->{res}) && isnum($n->{res})) {
      my $d = ($n->{res} - $s->{res}) / (abs($s->{res}) || 1);
      $va{$d < -1e-3 ? 'lower' : $d > 1e-3 ? 'higher' : 'same'}++; }
    if (isnum($s->{ndv}) && isnum($n->{ndv}) && $s->{ndv} > 0) { push @ratio, $n->{ndv} / $s->{ndv}; $s0 += $s->{ndv}; $s1 += $n->{ndv}; }
  }
  my @r = sort { $a <=> $b } @ratio; my $med = @r ? $r[int(@r/2)] : 0;
  printf "%-9s %4d | %5d %6d %6d %6d | %-24s | %-24s | %.2f, %d -> %d\n", $t, $c{rows} // 0, $c{both} // 0, $c{stock} // 0,
    $c{new} // 0, $c{neither} // 0, join(' ', map { "$_ " . ($th{$_} // 0) } qw(higher same lower)),
    join(' ', map { "$_ " . ($va{$_} // 0) } qw(lower same higher)), $med, $s0, $s1;
  print "            stock only: @lostk\n" if @lostk;
  print "            new only:   @gaink\n" if @gaink;
}
